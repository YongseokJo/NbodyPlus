#!/usr/bin/env bash
set -euo pipefail

source "$HOME/intel/oneapi/setvars.sh" >/dev/null 2>&1 || true
module add cuda/cuda-10.1.2 >/dev/null 2>&1 || true

# Avoid mixing OpenMPI (from modules) with Intel MPI (used by mpiicpx at link time).
if [[ -n "${I_MPI_ROOT:-}" ]]; then
	export LD_LIBRARY_PATH="$I_MPI_ROOT/lib:${LD_LIBRARY_PATH:-}"
	MPIRUN="$I_MPI_ROOT/bin/mpirun"
else
	MPIRUN="mpirun"
fi

ulimit -c unlimited || true

NP="${NP:-32}"
"$MPIRUN" -np "$NP" ./abyss.exe -c config.txt >stdout 2>stderr

# If a core file exists, write a quick backtrace.
CORE_FILE="$(ls -1 core core.* 2>/dev/null | head -n 1 || true)"
if [[ -n "$CORE_FILE" ]] && command -v gdb >/dev/null 2>&1; then
	gdb -q -batch -ex "thread apply all bt" ./abyss.exe "$CORE_FILE" >gdb_bt.txt 2>&1 || true
fi
