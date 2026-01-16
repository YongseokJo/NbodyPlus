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

TEST_ABS="$REPO_ROOT/$TEST_DIR"
CONFIG_ABS="$TEST_ABS/$RUN_CONFIG"

[[ -d "$TEST_ABS" ]] || workflow_die "TEST_DIR not found: $TEST_ABS"
[[ -f "$CONFIG_ABS" ]] || workflow_die "RUN_CONFIG not found: $CONFIG_ABS"

LOG="$RUN_DIR/run.log"
WORK_DIR="$RUN_DIR/work"
mkdir -p "$WORK_DIR"

# Stage test directory (inputs/config) into a per-run workspace.
# Note: this copies existing output files too; if that becomes large, we can switch to a filtered copy.
rm -rf "$WORK_DIR"/*
cp -a "$TEST_ABS/." "$WORK_DIR/"

# Stage executable
cp "$REPO_ROOT/$EXECUTABLE" "$WORK_DIR/abyss.exe"

{
  echo "== ABYSS run =="
  echo "date=$(date)"
  echo "work_dir=$WORK_DIR"
  echo "config=$RUN_CONFIG"
  echo "scheduler=$SCHEDULER"
  echo "ntasks=$NTASKS"
  echo ""
} > "$LOG"

pushd "$WORK_DIR" >/dev/null

set +e
if [[ "$SCHEDULER" == "slurm" ]]; then
  if [[ -n "${SRUN_MPI:-}" ]]; then
    srun --mpi="$SRUN_MPI" -n "$NTASKS" ./abyss.exe -c "$RUN_CONFIG" >>"$LOG" 2>&1
  else
    srun -n "$NTASKS" ./abyss.exe -c "$RUN_CONFIG" >>"$LOG" 2>&1
  fi
else
  mpirun -np "$NTASKS" ./abyss.exe -c "$RUN_CONFIG" >>"$LOG" 2>&1
fi
rc=$?
set -e

popd >/dev/null

if [[ $rc -ne 0 ]]; then
  workflow_die "Run failed (see $LOG)"
fi

# If the run produced an output directory in work/, keep it there; analyze.sh will summarize.
echo "Run OK" >> "$LOG"
