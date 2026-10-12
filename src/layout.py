"""Where things live on disk (`layout`, not `paths`: several modules already
take a local `paths` argument for a list of FITS files).

Two trees, and the split is the whole point: ``src/`` is what is written and
committed, ``data/`` is what a run produces -- the store, the star lists, the
light-curve cache, the stage-2 reports. ``data/`` is gitignored as one
directory, so nothing a pipeline writes can end up in a commit by accident,
and a checkout is clean without a list of per-file ignore patterns.

Every default path in the CLI comes from here, so moving the outputs is one
environment variable rather than a flag on each stage::

    STARSETI_DATA=D:/starseti uv run python src/stage1_pulse.py ...

Modules read these through the module (``layout.LC_CACHE_DIR``) rather than
importing the names, which is also what lets a test point one at tmp_path.
"""

from __future__ import annotations

import os
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SRC_DIR.parent

_env = os.environ.get("STARSETI_DATA")
DATA_DIR = (Path(_env) if _env else PROJECT_ROOT / "data").resolve()

STARLIST_DIR = DATA_DIR / "starlists"
LC_CACHE_DIR = DATA_DIR / "lc_cache"
RESULTS_DIR = DATA_DIR / "results"
DEFAULT_DB = DATA_DIR / "cepheid.sqlite"


def show(path: str | Path) -> str:
    """`path` written relative to the project when it is inside it.

    Help text and log lines want ``data/cepheid.sqlite``, not sixty
    characters of absolute path that differ on every machine.
    """
    p = Path(path)
    try:
        return p.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(p)


def starlist(name: str | Path) -> Path:
    """Resolve a target-file argument, looking in data/starlists/ as well.

    A path that exists is taken as given, so ``--target-file`` still accepts
    anything. A bare ``tess_cepheids.txt`` is found in the star-list
    directory, which is where stage 0 puts it. A name that is nowhere comes
    back unchanged, so the error names the file the caller actually asked for.
    """
    p = Path(name)
    if p.exists():
        return p
    candidate = STARLIST_DIR / p.name
    return candidate if candidate.exists() else p
