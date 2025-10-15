#!/bin/bash
# Helper script to run Python with conda libraries (avoids GLIBCXX issues)

CONDA_ENV="/home/RDC/simnacma/.conda/envs/dncit-paper"
export LD_LIBRARY_PATH="${CONDA_ENV}/lib:${LD_LIBRARY_PATH}"
export LD_PRELOAD="${CONDA_ENV}/lib/libstdc++.so.6:${CONDA_ENV}/lib/libgcc_s.so.1"

exec "${CONDA_ENV}/bin/python" "$@"