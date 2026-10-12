"""Process pools that can actually be stopped.

`ProcessPoolExecutor` as a context manager answers an exception -- including
the KeyboardInterrupt a Ctrl+C raises -- by *joining* its workers. A worker
inside a numba kernel, or waiting on a socket, does not stop for that, so on
Windows a Ctrl+C left the parent sitting in `shutdown(wait=True)` until every
running task finished on its own, and a task that never finished held the run
forever. Every pool in the pipeline goes through `worker_pool`, which on any
exception terminates the worker processes, cancels what was queued, and
re-raises; the interrupt then reaches the shell within a second.

`terminate_pool` is also what the injector's watchdog calls when a pass has
made no progress for its timeout: the hosts still in flight are written off,
and the processes that were holding them are killed rather than joined.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager

# One BLAS thread per worker. A worker's own numpy would otherwise open a
# thread pool the width of the machine -- sixteen processes each trying to
# fill sixteen cores, which on the least-squares fits in the template loop
# costs more in contention than the threads ever win back. The parent keeps
# whatever it started with (it does almost no linear algebra, and its numpy is
# loaded long before this runs); the children read these at their own import,
# because a spawned process inherits the environment as it stood here.
# setdefault, so anything the shell already set wins.
_THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def single_threaded_workers() -> None:
    """Pin worker BLAS to one thread each. Call before spawning any pool."""
    for var in _THREAD_VARS:
        os.environ.setdefault(var, "1")


def terminate_pool(pool: ProcessPoolExecutor) -> int:
    """Kill the pool's worker processes. Returns how many were signalled.

    Reaches into the executor's process table, which is private API, because
    the public one offers no way to stop a worker that will not return.
    """
    n = 0
    for proc in list(getattr(pool, "_processes", {}).values()):
        try:
            proc.terminate()
            n += 1
        except Exception:
            pass
    return n


@contextmanager
def worker_pool(max_workers: int) -> Iterator[ProcessPoolExecutor]:
    """A ProcessPoolExecutor that is torn down, not joined, on any exception.

    On a clean exit it behaves like the standard context manager and waits
    for the workers. On an exception -- Ctrl+C above all -- it terminates the
    workers, drops whatever was still queued, and lets the exception out.
    """
    single_threaded_workers()
    pool = ProcessPoolExecutor(max_workers=max(1, int(max_workers)))
    try:
        yield pool
    except BaseException:
        terminate_pool(pool)
        pool.shutdown(wait=False, cancel_futures=True)
        raise
    pool.shutdown(wait=True)
