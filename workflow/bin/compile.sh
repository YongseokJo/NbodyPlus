#!/usr/bin/env bash
set -euo pipefail

RUN_DIR=${1:?run_dir}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

workflow_load_config
workflow_setup_env
workflow_sanity

REPO_ROOT="$(workflow_repo_root)"
mkdir -p "$RUN_DIR"

LOG="$RUN_DIR/build.log"
{
  echo "== ABYSS compile =="
  echo "date=$(date)"
  echo "repo_root=$REPO_ROOT"
  echo "mpicxx=$(command -v mpicxx)"
  if [[ "${USE_CUDA:-0}" == "1" ]]; then
    echo "nvcc=$(command -v nvcc)"
    nvcc --version | head -4 || true
    echo "CUDAHOSTCXX=${CUDAHOSTCXX:-}"
  fi
  echo "HDF5_DIR=${HDF5_DIR:-}"
  echo ""
} > "$LOG"

pushd "$REPO_ROOT/src" >/dev/null

# Always clean to keep runs reproducible.
make clean >>"$LOG" 2>&1

MAKE_ARGS=()
if [[ "${USE_CUDA:-0}" == "1" ]]; then
  MAKE_ARGS+=("USE_CUDA=1")
fi
if [[ "${USE_SEVN:-0}" == "1" ]]; then
  MAKE_ARGS+=("USE_SEVN=1")
fi

# Build
set +e
make "${MAKE_ARGS[@]}" >>"$LOG" 2>&1
rc=$?
set -e

popd >/dev/null

if [[ $rc -ne 0 ]]; then
  workflow_die "Build failed (see $LOG)"
fi

if [[ ! -f "$REPO_ROOT/$EXECUTABLE" ]]; then
  workflow_die "Expected executable not found: $REPO_ROOT/$EXECUTABLE"
fi

echo "Build OK: $REPO_ROOT/$EXECUTABLE" >> "$LOG"
