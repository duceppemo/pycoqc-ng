# -*- coding: utf-8 -*-
"""The environment the converters' worker processes start with.

The workers parse text and HDF5, not linear algebra, yet numpy's OpenBLAS
starts a pool of threads per core as it loads, and the jemalloc bundled
with pyarrow, which pandas 3 loads, starts a background thread. In a fork
server those threads make every fork unsafe; in a spawned worker they are
hundreds of idle threads across a pool. Both read their settings from the
environment when they load, so the settings go in before any worker
imports numpy.

Nothing but the standard library may be imported here: the fork server
reads this before it imports numpy.
"""

WORKER_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "JE_ARROW_MALLOC_CONF": "background_thread:false",
}
