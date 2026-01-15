#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

load_config
maybe_load_modules

RUN_DIR=${1:-""}
if [[ -z "$RUN_DIR" ]]; then
  echo "Usage: compile.sh <run_dir>" >&2
  exit 2
fi

USE_CUDA_EFFECTIVE=$(detect_use_cuda)

cd "$ABYSS_ROOT/src"

echo "[compile] CXX=$CXX USE_CUDA=$USE_CUDA_EFFECTIVE"

if ! command -v "$CXX" >/dev/null 2>&1; then
  echo "[compile] ERROR: compiler wrapper not found: $CXX" >&2
  echo "[compile] Hint: set CXX=mpicxx (or load OpenMPI/IntelMPI modules)" >&2
  exit 2
fi

make clean || true

if [[ "$USE_CUDA_EFFECTIVE" == "1" ]]; then
  make -j8 CXX="$CXX" USE_CUDA=1 EXTRA_CXXFLAGS="$EXTRA_CXXFLAGS" EXTRA_NVCCFLAGS="$EXTRA_NVCCFLAGS" \
    >"$RUN_DIR/compile_stdout.txt" 2>"$RUN_DIR/compile_stderr.txt"
else
  # Important: Makefile uses `ifdef USE_CUDA`, so even USE_CUDA=0 would enable CUDA.
  # Ensure USE_CUDA is truly undefined for CPU-only builds.
  env -u USE_CUDA make -j8 CXX="$CXX" EXTRA_CXXFLAGS="$EXTRA_CXXFLAGS" EXTRA_NVCCFLAGS="$EXTRA_NVCCFLAGS" \
    >"$RUN_DIR/compile_stdout.txt" 2>"$RUN_DIR/compile_stderr.txt"
fi

if [[ ! -x "$ABYSS_ROOT/src/abyss.exe" ]]; then
  echo "[compile] ERROR: build did not produce src/abyss.exe" >&2
  exit 2
fi

echo "[compile] OK: $ABYSS_ROOT/src/abyss.exe"
