# -*- coding: utf-8 -*-
"""Imported once by the converters' fork server, never by the parent.

The workers are forked from the server, so whatever it imports is shared
with them for free. A server holding a thread, though, warns on every fork
and risks the deadlock it was meant to avoid, and two imports start threads
of their own: OpenBLAS starts a pool as soon as numpy loads, and the
jemalloc bundled with pyarrow, which pandas 3 loads, starts a background
thread. Settling both here, before either is imported, keeps the server
single threaded. The settings are in _worker_env.WORKER_ENV, and the
parent also holds them in its environment, for a moment, while it starts
the server and the workers (see _pipeline._worker_environment); its own
numpy loaded long before, so its thread pools are unaffected. The workers
parse text and HDF5, not linear algebra.
"""

import importlib
import os

from pycoQC._worker_env import WORKER_ENV

# The server is started with these already, by ConverterPipeline; set here
# too so the server stays single threaded however it was started.
os.environ.update(WORKER_ENV)

# What a converter's command line imports, so a worker importing the
# parent's main script finds it all loaded.
importlib.import_module("pycoQC.__main__")
