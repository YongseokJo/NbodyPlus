#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

load_config

RUN_DIR=${1:-""}
if [[ -z "$RUN_DIR" ]]; then
  echo "Usage: analyze.sh <run_dir>" >&2
  exit 2
fi

SUMMARY="$RUN_DIR/summary.txt"

{
  echo "ABYSS workflow summary"
  echo "timestamp: $(date)"
  echo "system: ${SYSTEM_NAME}"
  echo "scheduler: ${SCHEDULER}"
  echo "test_dir: ${TEST_DIR}"
  echo "test_config: ${TEST_CONFIG}"
  echo

  if [[ -f "$RUN_DIR/compile_stderr.txt" ]]; then
    echo "---- compile_stderr (first 80 lines) ----"
    sed -n '1,80p' "$RUN_DIR/compile_stderr.txt" || true
    echo
  fi

  if [[ -f "$RUN_DIR/run_stderr.txt" ]]; then
    echo "---- run_stderr (first 120 lines) ----"
    sed -n '1,120p' "$RUN_DIR/run_stderr.txt" || true
    echo
  fi

  if [[ -f "$RUN_DIR/run_stdout.txt" ]]; then
    echo "---- run_stdout (first 120 lines) ----"
    sed -n '1,120p' "$RUN_DIR/run_stdout.txt" || true
    echo
  fi

  echo "---- quick error scan ----"
  for f in "$RUN_DIR/compile_stderr.txt" "$RUN_DIR/run_stderr.txt" "$RUN_DIR/run_stdout.txt"; do
    [[ -f "$f" ]] || continue
    echo "file: $(basename "$f")"
    grep -Ein "error|fatal|segmentation|SIGSEGV|MPI_ERR|cannot open shared object|Abort|CUDA|cuCtx|no GPUs" "$f" | head -n 50 || true
    echo
  done

  if [[ -f "$RUN_DIR/output_files.txt" ]]; then
    echo "---- output files (first 80) ----"
    sed -n '1,80p' "$RUN_DIR/output_files.txt" || true
    echo
  fi
} > "$SUMMARY"

echo "[analyze] wrote $SUMMARY"
