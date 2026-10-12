#!/usr/bin/env python3
"""SQLite store used by every stage.

Canonical on-disk file is ``data/cepheid.sqlite``. Connections use WAL so a
second process can read (and, after a short wait, write) while stage 1
commits one star at a time. Callers may still write ``cep.table`` in SQL;
that prefix is stripped for SQLite.
"""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
from array import array
from collections.abc import Callable, Iterable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pandas as pd

import layout
from layout import DEFAULT_DB
from score import _json_obj, score_from_payloads
from sequence import SEQ_CODE_VERSION, is_square_match, score_oc_frame

BUSY_TIMEOUT_MS = 60_000

STAGE0_COLS = [
    "survey",
    "star_class",
    "target",
    "source_id",
    "type_best_classification",
    "type2_best_sub_classification",
    "mode_best_classification",
    "multi_mode_best_classification",
    "class_score",
    "period_days",
    "g_mag",
    "peak_to_peak_g",
    "ra",
    "dec",
    "l",
    "b",
    "host_galaxy",
    "tic",
    "tmag",
    "xmatch_sep_arcsec",
    "asas_sn_id",
    "asassn_tic_id",
    "asassn_match",
    "ztf_oid",
    "ztf_nobs",
    "ztf_nobs_total",
    "ztf_n_oids",
    "ztf_medianmag",
    "ztf_oids_all",
    "bp_rp",
]

OC_COLS = [
    "target",
    "survey",
    "cycle_index",
    "t_center",
    "shift_days",
    "shift_err_days",
    "npts",
    "chi2_red",
    "amp_ratio",
    "oc_detrended_days",
]

# The per-cycle series of a sparse survey (stage1_pulse.cycle_timings): one
# timing per *observed* cycle, from the epochs that fell in it, beside the
# block O-C in `stage1_oc`. It is what the epoch fold reads, and storing it
# is what lets `--recompute-routes` re-read that fold without a light curve.
EPOCH_OC_COLS = [
    "target",
    "survey",
    "cycle_index",
    "t_center",
    "shift_days",
    "shift_err_days",
    "npts",
    "chi2_red",
    "oc_detrended_days",
]


def add_db_arg(parser: argparse.ArgumentParser, default: Path = DEFAULT_DB) -> None:
    parser.add_argument(
        "--db", default=str(default), help=f"SQLite store; default: {layout.show(default)}"
    )


def ensure_sequence_columns(con: sqlite3.Connection) -> None:
    """ADD COLUMN for sequence fields introduced after the table already existed."""
    rows = con.execute("PRAGMA table_info(stage1_sequence)").fetchall()
    if not rows:
        return
    have = {r[1].lower() for r in rows}
    for name, decl in (
        ("seq_p", "DOUBLE"),
        ("seq_bit_len", "INTEGER"),
        ("seq_squareness", "DOUBLE"),
        ("seq_off_mask", "TEXT"),
        ("seq_code_version", "INTEGER"),
        ("seq_code_json", "TEXT"),
    ):
        if name not in have:
            con.execute(f"ALTER TABLE stage1_sequence ADD COLUMN {name} {decl}")


# The columns of a stage1_sequence row, in the order `seq_persist_tuple`
# writes them. Four places insert one of these rows; they share this list so a
# new column cannot reach three of them and be forgotten by the fourth.
SEQ_COLUMNS = (
    "target",
    "survey",
    "seq_matched",
    "seq_name",
    "seq_label",
    "seq_digit_count",
    "seq_sort_n",
    "oc_spark",
    "seq_p",
    "seq_bit_len",
    "seq_squareness",
    "seq_off_mask",
    "seq_code_version",
    "seq_code_json",
)


def seq_insert_sql(table: str = "cep.stage1_sequence", *, replace: bool = False) -> str:
    verb = "INSERT OR REPLACE" if replace else "INSERT"
    cols = ", ".join(SEQ_COLUMNS)
    marks = ", ".join("?" * len(SEQ_COLUMNS))
    return f"{verb} INTO {table} ({cols}) VALUES ({marks})"


def seq_persist_tuple(target: str, survey: str, seq: dict) -> tuple:
    # The decoded word goes in beside the fields read off it. It is what makes
    # a viewer request cheap: without it every look at a star re-runs the
    # decode, which on a few thousand cycles is seconds of description-length
    # search. The version stamp says which decoder wrote it (SEQ_CODE_VERSION).
    code = seq.get("code")
    return (
        target,
        survey,
        int(bool(seq.get("seq_matched"))),
        seq.get("seq_name") or "",
        seq.get("seq_label") or "unrecognized",
        int(seq.get("seq_digit_count") or 0),
        float(seq.get("seq_sort_n") or 0.0),
        seq.get("oc_spark") or "",
        float(seq.get("seq_p") if seq.get("seq_p") is not None else 1.0),
        int(seq.get("seq_bit_len") or 0),
        _float_or_none(seq.get("seq_squareness")),
        str(seq.get("seq_off_mask") or ""),
        int(SEQ_CODE_VERSION),
        None if not code else dumps_payload(code),
    )


# ---- phase fold -----------------------------------------------------------
# The fold is a display artifact: ~2500 points plus a 400-point template, and
# it was 96% of every stage-1 payload (990 MB of a 1.5 GB store). It lives in
# its own table as float32, so `payload` is small enough to read freely and the
# fold is only touched for the one star being looked at.


def pack_f32(values: Iterable[float | None] | None) -> bytes | None:
    """Little-endian float32 blob, or None for an empty series."""
    vals = [0.0 if v is None else float(v) for v in (values or [])]
    if not vals:
        return None
    buf = array("f", vals)
    if sys.byteorder != "little":
        buf.byteswap()
    return buf.tobytes()


def unpack_f32(blob: bytes | memoryview | None) -> list[float]:
    if not blob:
        return []
    buf = array("f")
    buf.frombytes(bytes(blob))
    if sys.byteorder != "little":
        buf.byteswap()
    return buf.tolist()


def fold_blobs(fold: dict | None) -> tuple | None:
    """(phase, flux, template_flux) blobs for one fold, or None.

    template_phase is linspace(0, 1, n) by construction (see
    stage1_pulse.pack_phase_fold), so it is rebuilt on read rather than stored.
    """
    if not fold:
        return None
    phase = pack_f32(fold.get("phase"))
    flux = pack_f32(fold.get("flux"))
    if phase is None or flux is None:
        return None
    return (phase, flux, pack_f32(fold.get("template_flux")))


def fold_from_blobs(
    phase: bytes | None, flux: bytes | None, template_flux: bytes | None
) -> dict | None:
    """Rebuild the dict shape the viewer expects."""
    ph = unpack_f32(phase)
    if not ph:
        return None
    tf = unpack_f32(template_flux)
    n = len(tf)
    return {
        "phase": ph,
        "flux": unpack_f32(flux),
        "template_phase": [i / (n - 1) for i in range(n)] if n > 1 else [],
        "template_flux": tf,
    }


def json_safe(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items() if not str(k).startswith("_")}
    if isinstance(obj, (list, tuple)):
        return [json_safe(v) for v in obj]
    if hasattr(obj, "tolist"):
        return json_safe(obj.tolist())
    if isinstance(obj, (float,)):
        return None if not math.isfinite(obj) else obj
    try:
        import numpy as np

        if isinstance(obj, (np.floating,)):
            return None if not math.isfinite(float(obj)) else float(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return json_safe(obj.tolist())
    except ImportError:
        pass
    return obj


def dumps_payload(obj: Any) -> str:
    return json.dumps(json_safe(obj), allow_nan=False, default=str)


# Fields that ride in a stage-1 record as far as the store and then go to a
# table of their own rather than into the summary payload: the phase fold (18
# kB of a 19 kB payload), the sequence row `analyse_oc` hands over so that no
# writer has to decode the word a second time, and the per-cycle series of a
# sparse survey (`epoch_oc`, a frame of up to a few thousand rows). Every
# writer filters through `payload_fields`, so a new one cannot quietly put
# them back.
FIELDS_WITH_OWN_TABLE = ("phase_fold", "sequence_row", "epoch_oc")


def payload_fields(rec: dict) -> dict:
    """`rec` without the fields the store keeps in their own tables."""
    return {k: v for k, v in rec.items() if k not in FIELDS_WITH_OWN_TABLE}


def _frame_rows(df: pd.DataFrame) -> list[tuple]:
    """Turn a DataFrame into sqlite-safe tuples (None for NA / non-finite)."""
    rows = []
    for rec in df.itertuples(index=False, name=None):
        out = []
        for v in rec:
            if v is None:
                out.append(None)
                continue
            try:
                if pd.isna(v):
                    out.append(None)
                    continue
            except (TypeError, ValueError):
                pass
            if isinstance(v, float) and not math.isfinite(v):
                out.append(None)
            elif hasattr(v, "item"):
                try:
                    item = v.item()
                except (ValueError, OverflowError):
                    item = None
                if isinstance(item, float) and not math.isfinite(item):
                    out.append(None)
                else:
                    out.append(item)
            else:
                out.append(v)
        rows.append(tuple(out))
    return rows


def connect_sqlite(path: str | Path, *, read_only: bool = False) -> sqlite3.Connection:
    """Open the store with WAL (writers) or a read-only URI (readers)."""
    path = Path(path)
    timeout_s = BUSY_TIMEOUT_MS / 1000
    if read_only:
        uri = f"file:{path.resolve().as_posix()}?mode=ro"
        con = sqlite3.connect(uri, uri=True, timeout=timeout_s)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        con = sqlite3.connect(str(path), timeout=timeout_s)
        con.execute("PRAGMA journal_mode=WAL").fetchone()
        con.execute("PRAGMA synchronous=NORMAL")
        con.execute("PRAGMA wal_autocheckpoint=1000")
    con.execute(f"PRAGMA busy_timeout={BUSY_TIMEOUT_MS}")
    con.execute("PRAGMA temp_store=MEMORY")
    con.isolation_level = None
    return con


def _sqlite_sql(sql: str) -> str:
    return sql.replace("cep.main.", "").replace("cep.", "")


class Store:
    """Native sqlite3 session. WAL + one transaction per upsert."""

    # WAL readers see a commit immediately; the explicit passive checkpoint
    # only folds pages back into the main file, and doing that after every
    # single star is pure extra I/O on a multi-gigabyte store.
    CHECKPOINT_EVERY = 64

    def __init__(self, path: str | Path = DEFAULT_DB) -> None:
        self.sqlite_path = Path(path)
        self.con: sqlite3.Connection | None = None
        self._commits = 0

    def __enter__(self) -> Store:
        self.con = connect_sqlite(self.sqlite_path)
        self._init_schema()
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        if self.con is not None:
            try:
                self.con.execute("PRAGMA wal_checkpoint(PASSIVE)").fetchone()
            except sqlite3.Error:
                pass
            self.con.close()
            self.con = None

    @property
    def path(self) -> Path:
        return self.sqlite_path

    def _init_schema(self) -> None:
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage0_targets (
                survey VARCHAR,
                star_class VARCHAR,
                target VARCHAR,
                source_id BIGINT,
                type_best_classification VARCHAR,
                type2_best_sub_classification VARCHAR,
                mode_best_classification VARCHAR,
                multi_mode_best_classification VARCHAR,
                class_score DOUBLE,
                period_days DOUBLE,
                g_mag DOUBLE,
                peak_to_peak_g DOUBLE,
                ra DOUBLE,
                dec DOUBLE,
                l DOUBLE,
                b DOUBLE,
                host_galaxy VARCHAR,
                tic BIGINT,
                tmag DOUBLE,
                xmatch_sep_arcsec DOUBLE,
                asas_sn_id BIGINT,
                asassn_tic_id BIGINT,
                asassn_match VARCHAR,
                ztf_oid BIGINT,
                ztf_nobs BIGINT,
                ztf_nobs_total BIGINT,
                ztf_n_oids BIGINT,
                ztf_medianmag DOUBLE,
                ztf_oids_all VARCHAR,
                bp_rp DOUBLE,
                PRIMARY KEY (survey, star_class, target)
            )
        """)
        self._ensure_stage0_columns()
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage1_summary (
                target VARCHAR,
                survey VARCHAR,
                star_class VARCHAR,
                mode VARCHAR,
                period_days DOUBLE,
                n_series INTEGER,
                oc_rms_seconds DOUBLE,
                p_ljung_box DOUBLE,
                verdict VARCHAR,
                error VARCHAR,
                in_triage BOOLEAN,
                payload VARCHAR,
                PRIMARY KEY (target, survey)
            )
        """)
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage1_oc (
                target VARCHAR,
                survey VARCHAR,
                cycle_index INTEGER,
                t_center DOUBLE,
                shift_days DOUBLE,
                shift_err_days DOUBLE,
                npts INTEGER,
                chi2_red DOUBLE,
                amp_ratio DOUBLE,
                oc_detrended_days DOUBLE
            )
        """)
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage1_epoch_oc (
                target VARCHAR,
                survey VARCHAR,
                cycle_index INTEGER,
                t_center DOUBLE,
                shift_days DOUBLE,
                shift_err_days DOUBLE,
                npts INTEGER,
                chi2_red DOUBLE,
                oc_detrended_days DOUBLE
            )
        """)
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage1_sequence (
                target VARCHAR,
                survey VARCHAR,
                seq_matched INTEGER,
                seq_name VARCHAR,
                seq_label VARCHAR,
                seq_digit_count INTEGER,
                seq_sort_n DOUBLE,
                oc_spark VARCHAR,
                seq_p DOUBLE,
                seq_bit_len INTEGER,
                seq_squareness DOUBLE,
                seq_off_mask VARCHAR,
                seq_code_version INTEGER,
                seq_code_json VARCHAR,
                PRIMARY KEY (target, survey)
            )
        """)
        self._ensure_sequence_columns()
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage1_fold (
                target VARCHAR,
                survey VARCHAR,
                phase BLOB,
                flux BLOB,
                template_flux BLOB,
                PRIMARY KEY (target, survey)
            )
        """)
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage1_score (
                target VARCHAR,
                survey VARCHAR,
                score DOUBLE,
                score_pct DOUBLE,
                classification VARCHAR,
                p_min DOUBLE,
                parts VARCHAR,
                stage2 INTEGER,
                PRIMARY KEY (target, survey)
            )
        """)
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.stage2_vetting (
                target VARCHAR,
                survey VARCHAR,
                star_class VARCHAR,
                classification VARCHAR,
                period_days DOUBLE,
                payload VARCHAR,
                PRIMARY KEY (target, survey)
            )
        """)
        # What `inject --find-threshold --publish` measured: the dt/P each star
        # type needs before the pipeline recovers it at some rate. It lives in
        # the main store rather than in injection.sqlite because the viewer
        # reads it -- the power a transmitter needs is a property of *this*
        # search's sensitivity, so a star's Kardashev class cannot be computed
        # without it. One row per (survey, class, word, cycles-per-bit, gate).
        self.execute("""
            CREATE TABLE IF NOT EXISTS cep.sensitivity (
                survey VARCHAR,
                star_class VARCHAR,
                word VARCHAR,
                cycles_per_bit INTEGER,
                gate VARCHAR,
                target DOUBLE,
                frac DOUBLE,
                is_bound INTEGER,
                n_hosts INTEGER,
                n_injections INTEGER,
                searched_hi DOUBLE,
                best_rate DOUBLE,
                run_id VARCHAR,
                measured_at DOUBLE,
                model VARCHAR,
                PRIMARY KEY (survey, star_class, word, cycles_per_bit, gate)
            )
        """)
        # The waveform a threshold was measured with. Rows published before it
        # was recorded were square-wave grids and arrive NULL; the viewer reads
        # NULL as `square`, so a stale figure is never passed off as the
        # physical model's.
        have = {r[1] for r in self.execute("PRAGMA cep.table_info(sensitivity)").fetchall()}
        if "model" not in have:
            self.execute("ALTER TABLE cep.sensitivity ADD COLUMN model VARCHAR")
        self._create_indexes()

    # ---- sensitivity ------------------------------------------------------

    SENSITIVITY_FIELDS = (
        "survey",
        "star_class",
        "word",
        "cycles_per_bit",
        "gate",
        "target",
        "frac",
        "is_bound",
        "n_hosts",
        "n_injections",
        "searched_hi",
        "best_rate",
        "run_id",
        "measured_at",
        "model",
    )

    def upsert_sensitivity(self, rows: Iterable[dict]) -> int:
        """Replace the measured threshold for each (survey, class, word, k, gate)."""
        payload = [tuple(r.get(f) for f in self.SENSITIVITY_FIELDS) for r in rows]
        if not payload:
            return 0
        assert self.con is not None
        cols = ", ".join(self.SENSITIVITY_FIELDS)
        marks = ", ".join("?" * len(self.SENSITIVITY_FIELDS))
        with self._transaction():
            self.con.executemany(
                f"INSERT OR REPLACE INTO sensitivity ({cols}) VALUES ({marks})", payload
            )
        return len(payload)

    def sensitivity(self) -> list[dict]:
        """Every measured threshold, newest first per cell."""
        cur = self.execute("SELECT * FROM sensitivity ORDER BY survey, star_class, gate")
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]

    def _create_indexes(self) -> None:
        for sql in (
            "CREATE INDEX IF NOT EXISTS idx_stage0_target ON stage0_targets (target)",
            "CREATE INDEX IF NOT EXISTS idx_stage0_star ON stage0_targets (target, survey)",
            "CREATE INDEX IF NOT EXISTS idx_s1_target ON stage1_summary (target)",
            "CREATE INDEX IF NOT EXISTS idx_oc_star ON stage1_oc (target, survey)",
            "CREATE INDEX IF NOT EXISTS idx_epoch_oc_star ON stage1_epoch_oc (target, survey)",
            "CREATE INDEX IF NOT EXISTS idx_fold_star ON stage1_fold (target, survey)",
            "CREATE INDEX IF NOT EXISTS idx_s2_target ON stage2_vetting (target)",
            "CREATE INDEX IF NOT EXISTS idx_seq_star ON stage1_sequence (target, survey)",
            "CREATE INDEX IF NOT EXISTS idx_score_star ON stage1_score (target, survey)",
            # Ranked paging reads straight off this index instead of sorting
            # the whole table for every page the viewer asks for.
            "CREATE INDEX IF NOT EXISTS idx_score_rank "
            "ON stage1_score (score DESC, p_min ASC, target ASC)",
            "CREATE INDEX IF NOT EXISTS idx_seq_rank "
            "ON stage1_sequence (seq_sort_n DESC, target ASC)",
            "CREATE INDEX IF NOT EXISTS idx_stage0_rank ON stage0_targets (g_mag ASC, target ASC)",
        ):
            try:
                self.execute(sql)
            except sqlite3.Error:
                pass

    def _ensure_stage0_columns(self) -> None:
        """ADD COLUMN for fields introduced after the table already existed."""
        have = {r[1].lower() for r in self.execute("PRAGMA table_info(stage0_targets)").fetchall()}
        wanted = (("bp_rp", "DOUBLE"),)
        for name, decl in wanted:
            if name not in have:
                self.execute(f"ALTER TABLE stage0_targets ADD COLUMN {name} {decl}")

    def _ensure_sequence_columns(self) -> None:
        assert self.con is not None
        ensure_sequence_columns(self.con)

    def execute(self, sql: str, params: Sequence[Any] | None = None) -> sqlite3.Cursor:
        assert self.con is not None
        sql = _sqlite_sql(sql)
        return self.con.execute(sql, params) if params is not None else self.con.execute(sql)

    def query_df(self, sql: str, params: Sequence[Any] | None = None) -> pd.DataFrame:
        assert self.con is not None
        sql = _sqlite_sql(sql)
        if params is None:
            return pd.read_sql_query(sql, self.con)
        return pd.read_sql_query(sql, self.con, params=list(params))

    @contextmanager
    def _transaction(self) -> Iterator[None]:
        assert self.con is not None
        nested = self.con.in_transaction
        if not nested:
            self.con.execute("BEGIN IMMEDIATE")
        try:
            yield
            if not nested:
                self.con.execute("COMMIT")
                self._commits += 1
                if self._commits % self.CHECKPOINT_EVERY == 0:
                    try:
                        self.con.execute("PRAGMA wal_checkpoint(PASSIVE)").fetchone()
                    except sqlite3.Error:
                        pass
        except Exception:
            if not nested:
                self.con.execute("ROLLBACK")
            raise

    def _insert_frame(self, table: str, frame: pd.DataFrame) -> None:
        if frame.empty:
            return
        cols = list(frame.columns)
        placeholders = ", ".join("?" * len(cols))
        sql = f"INSERT INTO {table} ({', '.join(cols)}) VALUES ({placeholders})"
        assert self.con is not None
        self.con.executemany(sql, _frame_rows(frame))

    def _write_score(self, target: str, survey: str, scored: dict) -> None:
        self.execute(
            """INSERT OR REPLACE INTO cep.stage1_score
               (target, survey, score, score_pct, classification, p_min, parts, stage2)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            [
                target,
                survey,
                scored["score"],
                scored["score_pct"],
                scored["classification"],
                scored["p_min"],
                json.dumps(scored["parts"]),
                int(bool(scored["stage2"])),
            ],
        )

    def _refresh_score(self, target: str, survey: str) -> None:
        s1 = self.execute(
            "SELECT error, in_triage, oc_rms_seconds, p_ljung_box, payload, mode "
            "FROM cep.stage1_summary WHERE target = ? AND survey = ?",
            [target, survey],
        ).fetchone()
        if not s1:
            return
        s1_row = {
            "error": s1[0],
            "in_triage": s1[1],
            "oc_rms_seconds": s1[2],
            "p_ljung_box": s1[3],
            "mode": s1[5],
        }
        seq = self.execute(
            "SELECT seq_p, seq_matched FROM cep.stage1_sequence WHERE target = ? AND survey = ?",
            [target, survey],
        ).fetchone()
        if seq:
            s1_row["seq_p"] = seq[0]
            s1_row["seq_matched"] = seq[1]
        s2 = self.execute(
            "SELECT classification, payload FROM cep.stage2_vetting "
            "WHERE target = ? AND survey = ?",
            [target, survey],
        ).fetchone()
        s2_row = {"classification": s2[0]} if s2 else None
        s2_payload = s2[1] if s2 else None
        self._write_score(target, survey, score_from_payloads(s1_row, s1[4], s2_row, s2_payload))

    # ---- stage 0 ------------------------------------------------------------

    def replace_stage0(self, survey: str, star_class: str, df: pd.DataFrame) -> int:
        hits = _stage0_frame(df, survey, star_class)
        with self._transaction():
            self.execute(
                "DELETE FROM cep.stage0_targets WHERE survey = ? AND star_class = ?",
                [survey, star_class],
            )
            self._insert_frame("stage0_targets", hits)
        return len(hits)

    def missing_bp_rp_sources(self) -> pd.DataFrame:
        """Distinct Gaia source_ids that still need a BP−RP color."""
        return self.query_df(
            """SELECT source_id, MIN(ra) AS ra, MIN(dec) AS dec
               FROM cep.stage0_targets
               WHERE source_id IS NOT NULL AND bp_rp IS NULL
               GROUP BY source_id
               ORDER BY source_id"""
        )

    def update_stage0_bp_rp(self, phot: pd.DataFrame) -> int:
        """Set bp_rp on every catalog row sharing a source_id. Returns rows updated."""
        if phot is None or phot.empty:
            return 0
        frame = phot[["source_id", "bp_rp"]].copy()
        frame["source_id"] = pd.to_numeric(frame["source_id"], errors="coerce").astype("Int64")
        frame["bp_rp"] = pd.to_numeric(frame["bp_rp"], errors="coerce")
        frame = frame.dropna(subset=["source_id", "bp_rp"]).drop_duplicates("source_id")
        if frame.empty:
            return 0
        before = int(
            self.execute("SELECT COUNT(*) FROM cep.stage0_targets WHERE bp_rp IS NULL").fetchone()[
                0
            ]
        )
        with self._transaction():
            assert self.con is not None
            self.con.executemany(
                "UPDATE stage0_targets SET bp_rp = ? WHERE source_id = ?",
                [(float(c), int(s)) for s, c in zip(frame["source_id"], frame["bp_rp"])],
            )
        after = int(
            self.execute("SELECT COUNT(*) FROM cep.stage0_targets WHERE bp_rp IS NULL").fetchone()[
                0
            ]
        )
        return max(0, before - after)

    def has_stage0(self, survey: str, star_class: str) -> bool:
        n = self.execute(
            "SELECT COUNT(*) FROM cep.stage0_targets WHERE survey = ? AND star_class = ?",
            [survey, star_class],
        ).fetchone()[0]
        return int(n) > 0

    def list_stage0_targets(
        self, surveys: list[str], star_class: str | None = None
    ) -> pd.DataFrame:
        """Unique (star_class, survey, target) rows from the stage-0 catalog."""
        if not surveys:
            return pd.DataFrame(columns=["star_class", "survey", "target", "period_days"])
        placeholders = ", ".join("?" * len(surveys))
        sql = (
            f"SELECT star_class, survey, target, MAX(period_days) AS period_days "
            f"FROM cep.stage0_targets "
            f"WHERE survey IN ({placeholders}) "
            f"AND target IS NOT NULL AND target != ''"
        )
        params: list = list(surveys)
        if star_class:
            sql += " AND star_class = ?"
            params.append(star_class)
        sql += " GROUP BY star_class, survey, target ORDER BY star_class, survey, target"
        return self.query_df(sql, params)

    def period_hints(self, targets: list[str], surveys: list[str]) -> dict[tuple[str, str], float]:
        """Gaia period_days keyed by (target, survey) for names already in the catalog."""
        if not targets or not surveys:
            return {}
        out: dict[tuple[str, str], float] = {}
        sph = ", ".join("?" * len(surveys))
        for i in range(0, len(targets), 500):
            chunk = [str(t) for t in targets[i : i + 500]]
            tph = ", ".join("?" * len(chunk))
            df = self.query_df(
                f"""SELECT target, survey, period_days FROM cep.stage0_targets
                    WHERE survey IN ({sph}) AND target IN ({tph})
                    AND period_days IS NOT NULL AND period_days > 0""",
                list(surveys) + chunk,
            )
            for row in df.itertuples(index=False):
                out[(str(row.target), str(row.survey))] = float(row.period_days)
        return out

    def ztf_fetch_hints(self, targets: list[str]) -> dict[str, dict]:
        """What stage 1 needs to fetch a ZTF star without asking IRSA where it is.

        `{"pos": (ra, dec), "oids": (oid, ...)}` per target. The oids are
        `ztf_oids_all`, the overlapping-field objects stage 0 already resolved
        -- the same set a cone search returns, and a cone search is the
        expensive call (measured 8-60 s against ~1 s for a query by id).
        """
        if not targets:
            return {}
        out: dict[str, dict] = {}
        for i in range(0, len(targets), 500):
            chunk = [str(t) for t in targets[i : i + 500]]
            tph = ", ".join("?" * len(chunk))
            df = self.query_df(
                f"""SELECT target, ra, dec, ztf_oids_all FROM cep.stage0_targets
                    WHERE survey = 'ztf' AND target IN ({tph})""",
                chunk,
            )
            for row in df.itertuples(index=False):
                oids = tuple(
                    int(o) for o in str(row.ztf_oids_all or "").split(";") if o.strip().isdigit()
                )
                pos = None
                if row.ra is not None and row.dec is not None:
                    try:
                        pos = (float(row.ra), float(row.dec))
                    except (TypeError, ValueError):
                        pos = None
                if pos is None and not oids:
                    continue
                out[str(row.target)] = {"pos": pos, "oids": oids}
        return out

    # ---- stage 1 ------------------------------------------------------------

    def stage1_skip_keys(self, retry_errors: bool = False) -> set[tuple[str, str]]:
        """(target, survey) pairs that --skip-existing should not re-run.

        Any prior stage-1 row counts as done. With retry_errors, a pair is kept
        on the work list when every existing row still has an error (so a later
        successful attempt in another run still wins).
        """
        df = self.query_df("SELECT target, survey, error FROM cep.stage1_summary")
        if df.empty:
            return set()
        have_ok: set[tuple[str, str]] = set()
        have_any: set[tuple[str, str]] = set()
        for r in df.itertuples(index=False):
            key = (str(r.target), str(r.survey))
            have_any.add(key)
            err = r.error
            if err is None or str(err).strip() == "":
                have_ok.add(key)
        return have_ok if retry_errors else have_any

    def upsert_stage1(self, res: dict, oc: pd.DataFrame | None, in_triage: bool) -> None:
        """Replace every prior stage-1 row for this (target, survey)."""
        target = str(res.get("target") or "")
        survey = str(res.get("survey") or "")
        # The fold goes to its own table; keeping it in the payload made every
        # payload read 18 kB instead of under 1 kB.
        fold = fold_blobs(res.get("phase_fold"))
        payload = dumps_payload(payload_fields(res))
        with self._transaction():
            for table in (
                "cep.stage1_summary",
                "cep.stage1_oc",
                "cep.stage1_epoch_oc",
                "cep.stage1_fold",
                "cep.stage1_sequence",
                "cep.stage1_score",
            ):
                self.execute(
                    f"DELETE FROM {table} WHERE target = ? AND survey = ?", [target, survey]
                )
            self.execute(
                """INSERT INTO cep.stage1_summary
                   (target, survey, star_class, mode, period_days, n_series,
                    oc_rms_seconds, p_ljung_box, verdict, error, in_triage, payload)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                [
                    target,
                    survey,
                    res.get("star_class"),
                    res.get("mode"),
                    _float_or_none(res.get("period_days")),
                    _int_or_none(res.get("n_series")),
                    _float_or_none(res.get("oc_rms_seconds")),
                    _float_or_none(res.get("p_ljung_box")),
                    res.get("verdict"),
                    res.get("error"),
                    int(bool(in_triage)),
                    payload,
                ],
            )
            if fold is not None:
                self.execute(
                    """INSERT INTO cep.stage1_fold
                       (target, survey, phase, flux, template_flux)
                       VALUES (?, ?, ?, ?, ?)""",
                    [target, survey, *fold],
                )
            if oc is not None and not oc.empty:
                self._insert_frame("stage1_oc", _oc_frame(oc, target, survey))
            epoch = res.get("epoch_oc")
            if epoch is not None and not getattr(epoch, "empty", True):
                self._insert_frame("stage1_epoch_oc", _epoch_frame(epoch, target, survey))
            # `analyse_oc` already decoded this series; take its row when the
            # record carries one and decode only for a caller that has none
            # (the injector, and a record that never reached the statistics).
            seq = res.get("sequence_row") or score_oc_frame(oc)
            self.execute(seq_insert_sql(), list(seq_persist_tuple(target, survey, seq)))
            self._refresh_score(target, survey)

    def count_stage1(self) -> int:
        return int(self.execute("SELECT COUNT(*) FROM cep.stage1_summary").fetchone()[0])

    def load_triage(self) -> list[dict]:
        """Stage-1 rows currently flagged in_triage.

        Columns only. The payload is ~18 kB per star (96% of it a phase fold)
        and stage 2 reads five fields out of it, so parsing every payload cost
        hundreds of megabytes of JSON before the first star was vetted.
        """
        cols = (
            "target",
            "survey",
            "star_class",
            "mode",
            "period_days",
            "oc_rms_seconds",
            "p_ljung_box",
            "verdict",
        )
        rows = self.execute(
            f"SELECT {', '.join(cols)} FROM cep.stage1_summary "
            "WHERE in_triage ORDER BY target, survey"
        ).fetchall()
        return [dict(zip(cols, r), in_triage=True) for r in rows]

    def stage2_skip_keys(self) -> set[tuple[str, str]]:
        """(target, survey) pairs that already have a stage-2 row."""
        rows = self.execute("SELECT target, survey FROM cep.stage2_vetting").fetchall()
        return {(str(t), str(s)) for t, s in rows}

    def load_stage2(self) -> list[dict]:
        rows = self.execute("SELECT payload FROM cep.stage2_vetting").fetchall()
        recs = [json.loads(r[0]) for r in rows if r[0]]
        recs.sort(key=lambda r: (str(r.get("target") or ""), str(r.get("survey") or "")))
        return recs

    def recompute_triage(self, pred: Callable[[dict], bool]) -> int:
        """Re-apply the triage predicate to every stored row.

        The sequence columns are joined in because they are stored per row
        while the payload may predate them: that is what lets a catalog-match
        route be applied to the whole store without reprocessing. A payload
        that already carries the field wins, since a reprocessed row measured
        it rather than inheriting it. The decode route cannot be recovered this
        way -- its null was never drawn for those rows -- so it needs a real
        stage-1 pass.

        The catalog route's window test lives in the stored code
        (`window_squareness`, stored since SEQ_CODE_VERSION 2). A row `recompute_sequence`
        has not rewritten since carries none, and the route stays closed to
        it until it does.
        """
        rows = self.execute(
            "SELECT s.target, s.survey, s.payload, q.seq_p, q.seq_squareness, q.seq_code_json "
            "FROM cep.stage1_summary s LEFT JOIN cep.stage1_sequence q "
            "ON q.target = s.target AND q.survey = s.survey"
        ).fetchall()
        n = 0
        with self._transaction():
            for target, survey, payload, seq_p, seq_squareness, code_json in rows:
                rec = json.loads(payload) if payload else {}
                if rec.get("seq_p") is None and seq_p is not None:
                    rec["seq_p"] = seq_p
                if rec.get("seq_squareness") is None and seq_squareness is not None:
                    rec["seq_squareness"] = seq_squareness
                if "seq_window_square" not in rec and code_json:
                    code = _json_obj(code_json)
                    rec["seq_window_squareness"] = code.get("window_squareness")
                    rec["seq_window_square"] = is_square_match(code)
                flag = bool(pred(rec))
                n += int(flag)
                self.execute(
                    "UPDATE cep.stage1_summary SET in_triage = ? WHERE target = ? AND survey = ?",
                    [int(flag), target, survey],
                )
                self._refresh_score(str(target), str(survey))
        return n

    # What `recompute_routes` works on: a star it can re-derive the statistics
    # for. A row with an error has no O-C series, and neither has one that
    # stopped at "fewer than 12 O-C points" or the amplitude floor -- together
    # 18,458 of 66,063 rows in the working store, which is why counting the
    # table instead of this made a finished run report 72%.
    _RECOMPUTE_FROM = """
        FROM cep.stage1_summary s
        WHERE (s.error IS NULL OR s.error = '') AND s.payload IS NOT NULL
          AND EXISTS (
              SELECT 1 FROM cep.stage1_oc o
              WHERE o.target = s.target AND o.survey = s.survey
          )
    """

    def count_recompute_routes(self) -> int:
        """How many stars `recompute_routes` will actually touch.

        The denominator a progress bar wants: counting `stage1_summary` counts
        rows the recompute skips without ever reaching the callback, so the
        bar stops short of the end on a run that did everything it could.
        """
        return int(self.execute(f"SELECT COUNT(*) {self._RECOMPUTE_FROM}").fetchone()[0])

    def recompute_routes(
        self,
        worker: Callable[[tuple], tuple],
        pred: Callable[[dict], bool],
        nshuffle: int,
        jobs: int = 1,
        progress: Callable[[int], None] | None = None,
    ) -> int:
        """Re-run the O-C statistics of every stored series and rewrite the row.

        `worker` takes (target, survey, record, oc frame, nshuffle, epoch
        frame or None) and returns
        (target, survey, (new record, [(cycle, detrended days)], sequence row)
        | None | ("error", msg)); it is a top-level function so a process pool
        can run it. The scored sequence row and the cycle indices come back
        from the worker rather than being decoded and re-read here, because
        this loop is the one serial thread of the pass and whatever is left in
        it is paid once per star with every core waiting. The store is touched
        from this process only: the payload, the summary columns the viewer and
        the injector read, the detrended O-C, the sequence row, the triage flag
        and the score are all rewritten together.
        """
        from pools import worker_pool

        stars = self.execute(
            f"SELECT s.target, s.survey, s.payload {self._RECOMPUTE_FROM} "
            "ORDER BY s.target, s.survey"
        ).fetchall()

        def jobs_iter():
            for target, survey, payload in stars:
                oc = self.load_oc(str(target), str(survey))
                if oc is None or oc.empty:
                    continue
                epoch = self.load_epoch_oc(str(target), str(survey))
                yield (str(target), str(survey), json.loads(payload), oc, nshuffle, epoch)

        def write(target: str, survey: str, out) -> bool:
            if not out or (isinstance(out, tuple) and out[0] == "error"):
                return False
            rec, det, seq = out
            flag = bool(pred(rec))
            with self._transaction():
                self.execute(
                    "UPDATE cep.stage1_summary SET payload = ?, verdict = ?, p_ljung_box = ?, "
                    "oc_rms_seconds = ?, in_triage = ? WHERE target = ? AND survey = ?",
                    [
                        dumps_payload(payload_fields(rec)),
                        rec.get("verdict"),
                        _float_or_none(rec.get("p_ljung_box")),
                        _float_or_none(rec.get("oc_rms_seconds")),
                        int(flag),
                        target,
                        survey,
                    ],
                )
                self.con.executemany(
                    "UPDATE stage1_oc SET oc_detrended_days = ? "
                    "WHERE target = ? AND survey = ? AND cycle_index = ?",
                    [(d, target, survey, c) for c, d in det],
                )
                self.execute(
                    "DELETE FROM cep.stage1_sequence WHERE target = ? AND survey = ?",
                    [target, survey],
                )
                self.execute(seq_insert_sql(), list(seq_persist_tuple(target, survey, seq)))
                self._refresh_score(target, survey)
            return True

        n = 0
        i = 0
        if jobs > 1:
            with worker_pool(jobs) as pool:
                for target, survey, out in pool.map(worker, jobs_iter(), chunksize=8):
                    n += int(write(target, survey, out))
                    i += 1
                    if progress is not None:
                        progress(i)
        else:
            for job in jobs_iter():
                target, survey, out = worker(job)
                n += int(write(target, survey, out))
                i += 1
                if progress is not None:
                    progress(i)
        return n

    def recompute_sequence(self, progress: Callable[[int], None] | None = None) -> int:
        """Re-decode every stored O-C and rewrite its stage1_sequence row.

        A change to the sequence scoring cannot be migrated from the stored row
        the way `recompute_triage` migrates the triage flag: `seq_p` depends on
        where in a catalog the word matched, and a deeper catalog can match a
        word that missed before, so the word has to be decoded again. The O-C
        is in the store, so this needs no light curves and no network.

        Run it before `recompute_triage` after a scoring change, since the
        triage gate reads `seq_p`.
        """
        stars = self.execute(
            "SELECT target, survey FROM cep.stage1_summary "
            "WHERE error IS NULL OR error = '' ORDER BY target, survey"
        ).fetchall()
        n = 0
        for i, (target, survey) in enumerate(stars):
            oc = self.load_oc(str(target), str(survey))
            if oc is None or oc.empty:
                continue
            seq = score_oc_frame(oc)
            with self._transaction():
                self.execute(
                    "DELETE FROM cep.stage1_sequence WHERE target = ? AND survey = ?",
                    [target, survey],
                )
                self.execute(
                    seq_insert_sql(), list(seq_persist_tuple(str(target), str(survey), seq))
                )
                self._refresh_score(str(target), str(survey))
            n += 1
            if progress is not None:
                progress(i)
        return n

    def stage1_retry_rows(self) -> pd.DataFrame:
        """Stage-1 rows that died in download or analysis."""
        return self.query_df(
            """SELECT target, survey, star_class, error
               FROM cep.stage1_summary
               WHERE error IS NOT NULL AND (
                   error LIKE 'download failed:%' OR error LIKE 'analysis failed:%'
               )
               ORDER BY star_class, survey, target"""
        )

    def load_fold(self, target: str, survey: str) -> dict | None:
        """The stored phase fold, or None. Falls back to a pre-migration payload."""
        row = self.execute(
            "SELECT phase, flux, template_flux FROM cep.stage1_fold "
            "WHERE target = ? AND survey = ?",
            [target, survey],
        ).fetchone()
        if row:
            fold = fold_from_blobs(*row)
            if fold:
                return fold
        raw = self.execute(
            "SELECT payload FROM cep.stage1_summary WHERE target = ? AND survey = ?",
            [target, survey],
        ).fetchone()
        return _json_obj(raw[0] if raw else None).get("phase_fold")

    def load_oc(self, target: str, survey: str) -> pd.DataFrame:
        return self.query_df(
            """SELECT cycle_index, t_center, shift_days, shift_err_days, npts, chi2_red,
                      amp_ratio, oc_detrended_days
               FROM cep.stage1_oc
               WHERE target = ? AND survey = ?
               ORDER BY cycle_index""",
            [target, survey],
        )

    def load_epoch_oc(self, target: str, survey: str) -> pd.DataFrame | None:
        """The stored per-cycle series of a sparse star, or None for a dense
        one or a row written before the table existed."""
        df = self.query_df(
            """SELECT cycle_index, t_center, shift_days, shift_err_days, npts, chi2_red,
                      oc_detrended_days
               FROM cep.stage1_epoch_oc
               WHERE target = ? AND survey = ?
               ORDER BY cycle_index""",
            [target, survey],
        )
        return None if df.empty else df

    # ---- stage 2 ------------------------------------------------------------

    def upsert_stage2(self, rec: dict) -> None:
        target = str(rec.get("target") or "")
        survey = str(rec.get("survey") or "")
        with self._transaction():
            self.execute(
                "DELETE FROM cep.stage2_vetting WHERE target = ? AND survey = ?", [target, survey]
            )
            self.execute(
                """INSERT INTO cep.stage2_vetting
                   (target, survey, star_class, classification, period_days, payload)
                   VALUES (?, ?, ?, ?, ?, ?)""",
                [
                    target,
                    survey,
                    rec.get("star_class"),
                    rec.get("classification"),
                    _float_or_none(rec.get("period_days")),
                    dumps_payload(rec),
                ],
            )
            self._refresh_score(target, survey)

    # ---- inspect ------------------------------------------------------------

    def inspect(self, name: str) -> dict:
        variants = _name_variants(name)
        digits = _digits(name)
        keys = _sql_in(variants)
        num_keys = _sql_in([digits] if digits else variants)
        stage0 = self.query_df(
            f"""SELECT * FROM cep.stage0_targets
                WHERE upper(replace(target, '_', ' ')) IN ({keys})
                   OR CAST(tic AS VARCHAR) IN ({num_keys})
                   OR CAST(asas_sn_id AS VARCHAR) IN ({num_keys})
                   OR CAST(ztf_oid AS VARCHAR) IN ({num_keys})
                   OR CAST(source_id AS VARCHAR) IN ({num_keys})
                ORDER BY survey, star_class, target"""
        )
        s1 = self.query_df(
            f"""SELECT target, survey, star_class, mode, period_days, n_series,
                       oc_rms_seconds, p_ljung_box, verdict, error, in_triage, payload
                FROM cep.stage1_summary
                WHERE upper(replace(target, '_', ' ')) IN ({keys})
                ORDER BY survey, target"""
        )
        oc_counts = self.query_df(
            f"""SELECT target, survey, COUNT(*) AS n_oc
                FROM cep.stage1_oc
                WHERE upper(replace(target, '_', ' ')) IN ({keys})
                GROUP BY 1, 2"""
        )
        s2 = self.query_df(
            f"""SELECT target, survey, star_class, classification, period_days, payload
                FROM cep.stage2_vetting
                WHERE upper(replace(target, '_', ' ')) IN ({keys})
                ORDER BY survey, target"""
        )
        return {
            "query": name,
            "variants": variants,
            "stage0": stage0,
            "stage1": s1,
            "oc_counts": oc_counts,
            "stage2": s2,
        }


def _stage0_frame(df: pd.DataFrame, survey: str, star_class: str) -> pd.DataFrame:
    out = df.copy()
    out["survey"] = survey
    out["star_class"] = star_class
    if "target" not in out.columns or out["target"].isna().all():
        if survey == "tess":
            out["target"] = "TIC " + pd.to_numeric(out["tic"], errors="coerce").astype(
                "Int64"
            ).astype(str)
        elif survey == "asassn":
            out["target"] = "ASASSN " + pd.to_numeric(out["asas_sn_id"], errors="coerce").astype(
                "Int64"
            ).astype(str)
        elif survey == "ztf":
            out["target"] = "ZTF " + pd.to_numeric(out["ztf_oid"], errors="coerce").astype(
                "Int64"
            ).astype(str)
        else:
            out["target"] = ""
    for col in STAGE0_COLS:
        if col not in out.columns:
            out[col] = pd.NA
    out = out[STAGE0_COLS]
    for col in (
        "source_id",
        "tic",
        "asas_sn_id",
        "asassn_tic_id",
        "ztf_oid",
        "ztf_nobs",
        "ztf_nobs_total",
        "ztf_n_oids",
    ):
        out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")
    for col in (
        "class_score",
        "period_days",
        "g_mag",
        "peak_to_peak_g",
        "bp_rp",
        "ra",
        "dec",
        "l",
        "b",
        "tmag",
        "xmatch_sep_arcsec",
        "ztf_medianmag",
    ):
        out[col] = pd.to_numeric(out[col], errors="coerce")
    for col in (
        "survey",
        "star_class",
        "target",
        "type_best_classification",
        "type2_best_sub_classification",
        "mode_best_classification",
        "multi_mode_best_classification",
        "host_galaxy",
        "asassn_match",
        "ztf_oids_all",
    ):
        out[col] = out[col].astype("string")
    return out.reset_index(drop=True)


def _oc_frame(df: pd.DataFrame, target: str, survey: str) -> pd.DataFrame:
    out = df.copy()
    if "cycle_index" not in out.columns and "index" in out.columns:
        out = out.rename(columns={"index": "cycle_index"})
    if "cycle_index" not in out.columns and "cycle" in out.columns:
        out = out.rename(columns={"cycle": "cycle_index"})
    out["target"] = target
    out["survey"] = survey
    for col in OC_COLS:
        if col not in out.columns:
            out[col] = pd.NA
    return out[OC_COLS].reset_index(drop=True)


def _epoch_frame(df: pd.DataFrame, target: str, survey: str) -> pd.DataFrame:
    out = df.copy()
    if "cycle_index" not in out.columns and "index" in out.columns:
        out = out.rename(columns={"index": "cycle_index"})
    out["target"] = target
    out["survey"] = survey
    for col in EPOCH_OC_COLS:
        if col not in out.columns:
            out[col] = pd.NA
    return out[EPOCH_OC_COLS].reset_index(drop=True)


def _float_or_none(val: object) -> float | None:
    if val is None or val == "":
        return None
    try:
        x = float(val)
    except (TypeError, ValueError):
        return None
    return None if not math.isfinite(x) else x


def _int_or_none(val: object) -> int | None:
    if val is None or val == "":
        return None
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def _name_variants(name: str) -> list[str]:
    raw = name.strip()
    spaced = " ".join(raw.replace("_", " ").split())
    upper = spaced.upper()
    out = {raw, spaced, upper}
    digits = _digits(raw)
    if digits:
        out.update({digits, f"TIC {digits}", f"ASASSN {digits}", f"ZTF {digits}"})
    if upper.startswith("TIC"):
        rest = upper[3:].strip()
        if rest:
            out.add(f"TIC {rest.lstrip('0') or rest}")
    return [v for v in out if v]


def _digits(name: str) -> str:
    import re

    m = re.search(r"(\d{3,})", name)
    return m.group(1).lstrip("0") or m.group(1) if m else ""


def _sql_in(values: list[str]) -> str:
    if not values:
        return "''"
    return ", ".join("'" + v.replace("'", "''").upper().replace("_", " ") + "'" for v in values)
