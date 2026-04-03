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

# Check if IC file was pre-generated (by separate McLuster job)
IC_FILE="$WORK_DIR/mcluster_abyss.dat"
HAS_PREGENERATED_IC=0
if [[ -f "$IC_FILE" ]]; then
  HAS_PREGENERATED_IC=1
  echo "Found pre-generated IC file: $IC_FILE"
fi

# Stage test directory (inputs/config) into a per-run workspace.
# Preserve pre-generated IC file if it exists.
if [[ "$HAS_PREGENERATED_IC" -eq 1 ]]; then
  # Save IC file, clear work dir, restore IC file
  IC_BACKUP=$(mktemp)
  cp "$IC_FILE" "$IC_BACKUP"
  rm -rf "$WORK_DIR"/*
  cp -a "$TEST_ABS/." "$WORK_DIR/"
  mv "$IC_BACKUP" "$IC_FILE"
else
  rm -rf "$WORK_DIR"/*
  cp -a "$TEST_ABS/." "$WORK_DIR/"
fi

# Stage executable
cp "$REPO_ROOT/$EXECUTABLE" "$WORK_DIR/abyss.exe"

# Stage McLuster only if no pre-generated IC and config uses [mcluster] section
if [[ "$HAS_PREGENERATED_IC" -eq 0 ]] && grep -q '^\[mcluster\]' "$WORK_DIR/$RUN_CONFIG" 2>/dev/null; then
  MCLUSTER_BIN="$REPO_ROOT/mcluster/mcluster_sse"
  if [[ -x "$MCLUSTER_BIN" ]]; then
    cp "$MCLUSTER_BIN" "$WORK_DIR/mcluster"
    echo "McLuster staged for IC generation"
  else
    echo "Warning: Config uses [mcluster] but mcluster_sse not found at $MCLUSTER_BIN" >&2
    echo "  Build with USE_MCLUSTER=1 or provide pre-generated IC file" >&2
  fi
fi

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
