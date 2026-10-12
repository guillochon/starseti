"""Raw ASAS-SN / ZTF photometry on disk, so a re-run costs no network.

TESS has had a light-curve cache since the beginning, but nobody wrote one:
MAST serves *files*, `lightkurve` downloads them to
``data/lc_cache/mastDownload/TESS/`` because that is the only way to read the
archive, and the cache is a side effect of that. ASAS-SN and ZTF are
API-shaped -- Sky Patrol hands back an in-memory frame, IRSA hands back a CSV
body -- so nothing was ever persisted, and every stage-1 pass over the 29,746
sparse stars in the store re-downloaded all of them. That is also what keeps
the injection grids on TESS: ``inject.cache_index`` draws hosts by globbing
FITS off disk, and for two thirds of the sample there is no disk to glob.

**Raw rows, not light curves.** What is stored is the frame the API returned,
before any cut -- the same contract the TESS FITS have, where every quality
mask, band selection and outlier clip is applied on *read*. So changing a
cleaning rule, a band or ``include_poor`` costs nothing, and a cached blob for
one band serves every other. The cleaning code in ``stage1_pulse`` is
unchanged and still runs on every load; this layer only decides whether the
bytes come from IRSA or from disk.

**Only successful, non-empty queries are cached.** ``_ztf_get`` collapses "the
survey has nothing here" and "three timeouts in a row" into the same ``None``,
and persisting the second as the first is precisely the bug that left the
whole ZTF arm reading "no data" for months (cdb7b9d). Negative caching is
worth having -- it needs that distinction made first, and is not made here.

**Entries expire.** TESS sectors are immutable once released; ASAS-SN and ZTF
keep appending epochs to the same object, so a cached blob is a snapshot and
not an archive. Every entry carries the time it was fetched and is refused
past ``DEFAULT_MAX_AGE_DAYS``, which ``--lc-cache-age`` overrides and
``--lc-cache-refresh`` ignores entirely.

Layout, under ``layout.LC_CACHE_DIR`` alongside ``mastDownload/``::

    asassn/3f/asassn_498216382399.3f9c1a2b4d5e.parquet   rows
    asassn/3f/asassn_498216382399.3f9c1a2b4d5e.json      target, params, fetched_at

The JSON sidecar is what `index()` reads, so listing 30k entries costs no
parquet footers and a cache directory can be understood with `cat`. The
two-hex shard keeps any one directory near a thousand files.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

import layout

# Bumped when the stored shape changes meaning. It is part of the key, so a
# bump orphans every old entry rather than reading it as if it were current;
# `prune(orphans=True)` is what removes them.
CACHE_VERSION = 1

# ASAS-SN and ZTF both append to an existing object, so an entry is a snapshot.
# 30 days is a compromise: a survey adds a handful of epochs to a star in that
# window, which cannot move a period fit built on years of baseline, and it
# keeps a month-long reprocessing campaign off the network entirely.
DEFAULT_MAX_AGE_DAYS = 30.0

# Which surveys have a cache here at all. TESS is deliberately absent: its
# cache is the FITS tree, and a second copy of it would be 7.3 GB of nothing.
SURVEYS = ("asassn", "ztf")

_SLUG_RE = re.compile(r"[^a-z0-9]+")
_SLUG_MAX = 48

# Module-level off switch, for a run that must see the live survey.
# `STARSETI_NO_LC_CACHE=1` sets it at import; `--lc-cache-refresh` passes
# `refresh` to the read instead, which still *writes* what it fetched.
ENABLED = os.environ.get("STARSETI_NO_LC_CACHE", "") not in ("1", "true", "yes")


@dataclass(frozen=True)
class CachePolicy:
    """How one run treats the cache. An argument, never module state.

    Windows spawns its workers, so a setting made on the module in the parent
    is simply absent in the child -- that is the `--prewhiten` bug (3.4), and
    this is the shape that fixed it: the policy travels with the work item.

    `refresh` re-fetches and *rewrites* rather than bypassing: a run told to
    ignore what is stored still leaves the cache better than it found it.
    """

    enabled: bool = True
    max_age_days: float | None = DEFAULT_MAX_AGE_DAYS
    refresh: bool = False


DEFAULT_POLICY = CachePolicy()
OFF = CachePolicy(enabled=False)


def root() -> Path:
    """The cache root, read through `layout` so a test can point it at tmp_path."""
    return Path(layout.LC_CACHE_DIR)


def _slug(target: str) -> str:
    """A filesystem-safe stem for a target name, for humans only.

    Target names carry spaces and dots ('ASASSN 498216382399', 'RADEC 10.0
    20.0', 'delta Cep'). The slug is lossy and never parsed back -- the key
    hash is what identifies an entry, and the JSON sidecar holds the real name.
    """
    s = _SLUG_RE.sub("_", target.strip().lower()).strip("_")
    return (s[:_SLUG_MAX] or "target").rstrip("_")


def key_for(survey: str, target: str, params: dict[str, Any] | None = None) -> str:
    """Stable 12-hex identity for one query.

    Everything that changes *what the API is asked for* belongs in `params`;
    everything applied to the rows afterwards must not, or a band cut would
    fork the cache into copies of the same download.
    """
    payload = json.dumps(
        {
            "v": CACHE_VERSION,
            "survey": survey,
            "target": target,
            "params": dict(sorted((params or {}).items())),
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]


def path_for(survey: str, target: str, params: dict[str, Any] | None = None) -> Path:
    key = key_for(survey, target, params)
    return root() / survey / key[:2] / f"{_slug(target)}.{key}.parquet"


def _meta_path(data_path: Path) -> Path:
    return data_path.with_suffix(".json")


def has(
    survey: str,
    target: str,
    params: dict[str, Any] | None = None,
    *,
    policy: CachePolicy | None = None,
) -> bool:
    """Whether `load` would return a frame, without reading the frame.

    Same answer as `load(...) is not None` for every reason that depends on
    the sidecar -- caching off, nothing stored, aged out -- and it stops
    there. `load` goes on to parse the parquet, which is what a caller
    deciding *what to download* pays for and never uses: filling the cache for
    a 62,000-star target list meant 21,702 full parquet parses before the
    first fetch. A corrupt or empty parquet still reports True here and is
    re-fetched when something actually reads it.
    """
    pol = policy or DEFAULT_POLICY
    if not ENABLED or not pol.enabled or pol.refresh:
        return False
    path = path_for(survey, target, params)
    meta_path = _meta_path(path)
    if not path.exists() or not meta_path.exists():
        return False
    if pol.max_age_days is None:
        return True
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return (time.time() - float(meta.get("fetched_at", 0.0))) / 86400.0 <= pol.max_age_days


def load(
    survey: str,
    target: str,
    params: dict[str, Any] | None = None,
    *,
    policy: CachePolicy | None = None,
) -> pd.DataFrame | None:
    """The cached raw frame, or None to go to the network.

    None for every reason a caller would rather not distinguish: caching off,
    nothing stored, the entry aged out, or the file is unreadable. A corrupt
    entry is a cache miss, never an exception -- a half-written parquet from
    an interrupted run must cost a download, not a crash.
    """
    pol = policy or DEFAULT_POLICY
    if not ENABLED or not pol.enabled or pol.refresh:
        return None
    path = path_for(survey, target, params)
    meta_path = _meta_path(path)
    if not path.exists() or not meta_path.exists():
        return None
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if pol.max_age_days is not None:
        age_days = (time.time() - float(meta.get("fetched_at", 0.0))) / 86400.0
        if age_days > pol.max_age_days:
            return None
    try:
        df = pd.read_parquet(path)
    except Exception:
        return None
    return df if len(df) else None


def save(
    survey: str,
    target: str,
    df: pd.DataFrame,
    params: dict[str, Any] | None = None,
    *,
    policy: CachePolicy | None = None,
) -> Path | None:
    """Persist one successful, non-empty query. Returns the path, or None.

    Never raises: a cache that cannot be written is a performance problem, and
    taking down a twelve-hour stage-1 run over a full disk would be a
    correctness one. Writes to a temporary name and renames, so a reader never
    sees a partial file, and writes the data before the sidecar, so an
    interrupted write reads as a miss rather than as empty rows.
    """
    pol = policy or DEFAULT_POLICY
    if not ENABLED or not pol.enabled or df is None or not len(df):
        return None
    path = path_for(survey, target, params)
    meta_path = _meta_path(path)
    tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    meta_tmp = meta_path.with_name(f"{meta_path.name}.{os.getpid()}.tmp")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(tmp, index=False)
        os.replace(tmp, path)
        meta = {
            "cache_version": CACHE_VERSION,
            "survey": survey,
            "target": target,
            "params": dict(sorted((params or {}).items())),
            "fetched_at": time.time(),
            "n_rows": len(df),
            "columns": [str(c) for c in df.columns],
        }
        meta_tmp.write_text(json.dumps(meta, indent=1), encoding="utf-8")
        os.replace(meta_tmp, meta_path)
    except Exception:
        for leftover in (tmp, meta_tmp):
            try:
                leftover.unlink(missing_ok=True)
            except OSError:
                pass
        return None
    return path


def entries(
    survey: str | None = None, *, base_dir: str | Path | None = None
) -> Iterator[tuple[Path, dict[str, Any]]]:
    """(parquet path, metadata) for every readable entry.

    Driven by the sidecars: a parquet with no sidecar is an interrupted write
    and a sidecar with no parquet is a deleted one, and neither is an entry.
    """
    top = Path(base_dir) if base_dir is not None else root()
    for name in (survey,) if survey else SURVEYS:
        base = top / name
        if not base.is_dir():
            continue
        for meta_path in base.rglob("*.json"):
            data_path = meta_path.with_suffix(".parquet")
            if not data_path.exists():
                continue
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            yield data_path, meta


def index(
    survey: str | None = None,
    *,
    max_age_days: float | None = None,
    base_dir: str | Path | None = None,
) -> dict[str, list[str]]:
    """{target: [paths]} over the cache, the shape `inject.cache_index` returns.

    A target with more than one entry has been fetched under different query
    parameters (ZTF with and without field merging); all of them are listed,
    newest first, so a caller that wants one can take the first.
    """
    now = time.time()
    out: dict[str, list[tuple[float, str]]] = {}
    for path, meta in entries(survey, base_dir=base_dir):
        target = str(meta.get("target") or "")
        if not target:
            continue
        fetched = float(meta.get("fetched_at", 0.0))
        if max_age_days is not None and (now - fetched) / 86400.0 > max_age_days:
            continue
        out.setdefault(target, []).append((fetched, str(path)))
    return {t: [p for _, p in sorted(v, reverse=True)] for t, v in sorted(out.items())}


def stats(survey: str | None = None) -> dict[str, dict[str, float]]:
    """Per-survey entry count, target count, rows, bytes and fetch times."""
    out: dict[str, dict[str, float]] = {}
    for path, meta in entries(survey):
        name = str(meta.get("survey") or "?")
        row = out.setdefault(
            name,
            {"entries": 0, "targets": 0, "rows": 0, "bytes": 0.0, "oldest": 0.0, "newest": 0.0},
        )
        row["entries"] += 1
        row["rows"] += int(meta.get("n_rows", 0))
        try:
            row["bytes"] += path.stat().st_size
        except OSError:
            pass
        fetched = float(meta.get("fetched_at", 0.0))
        row["oldest"] = fetched if not row["oldest"] else min(row["oldest"], fetched)
        row["newest"] = max(row["newest"], fetched)
    for name in out:
        out[name]["targets"] = len(index(name))
    return out


def prune(
    survey: str | None = None,
    *,
    max_age_days: float | None = None,
    orphans: bool = True,
    dry_run: bool = False,
) -> list[Path]:
    """Remove aged-out entries, and files that are not part of one.

    Orphans are the residue of interrupted writes and of a `CACHE_VERSION`
    bump: a parquet whose sidecar is missing, a sidecar whose parquet is, and
    any `.tmp` left by a process that died mid-write.
    """
    now = time.time()
    removed: list[Path] = []
    for name in (survey,) if survey else SURVEYS:
        base = root() / name
        if not base.is_dir():
            continue
        if max_age_days is not None:
            for path, meta in list(entries(name)):
                fetched = float(meta.get("fetched_at", 0.0))
                if (now - fetched) / 86400.0 > max_age_days:
                    removed += [path, _meta_path(path)]
        if orphans:
            for f in base.rglob("*"):
                if f.is_dir() or f in removed:
                    continue
                if f.name.endswith(".tmp"):
                    removed.append(f)
                elif f.suffix == ".parquet" and not _meta_path(f).exists():
                    removed.append(f)
                elif f.suffix == ".json" and not f.with_suffix(".parquet").exists():
                    removed.append(f)
    if not dry_run:
        for f in removed:
            try:
                f.unlink(missing_ok=True)
            except OSError:
                pass
    return removed
