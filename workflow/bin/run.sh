#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

load_config
maybe_load_modules

RUN_DIR=${1:-""}
if [[ -z "$RUN_DIR" ]]; then
  echo "Usage: run.sh <run_dir>" >&2
  exit 2
fi

TEST_ABS="$ABYSS_ROOT/$TEST_DIR"
if [[ ! -d "$TEST_ABS" ]]; then
  echo "[run] ERROR: missing test dir: $TEST_ABS" >&2
  exit 2
fi

cd "$TEST_ABS"
cp -f "$ABYSS_ROOT/src/abyss.exe" ./abyss.exe

mkdir -p output

echo "[run] test_dir=$TEST_DIR config=$TEST_CONFIG" | tee "$RUN_DIR/run_meta.txt"

RUN_STDOUT="$RUN_DIR/run_stdout.txt"
RUN_STDERR="$RUN_DIR/run_stderr.txt"
rm -f "$RUN_STDOUT" "$RUN_STDERR"

if [[ -n "${SLURM_JOB_ID:-}" ]] && [[ "$MPI_LAUNCHER" == "srun" ]] && command -v srun >/dev/null 2>&1; then
  tasks="${SLURM_NTASKS:-$NTASKS}"
  tasks_per_node="${SLURM_NTASKS_PER_NODE:-$tasks}"
  echo "[run] launching with srun (ntasks=$tasks)" | tee -a "$RUN_DIR/run_meta.txt"
  srun --mpi="$SRUN_MPI" -n "$tasks" --ntasks-per-node="$tasks_per_node" \
    --output="$RUN_STDOUT" --error="$RUN_STDERR" \
    ./abyss.exe -c "$TEST_CONFIG" || true
elif [[ "$MPI_LAUNCHER" == "mpirun" ]] && command -v mpirun >/dev/null 2>&1; then
  echo "[run] launching with mpirun (ntasks=$NTASKS)" | tee -a "$RUN_DIR/run_meta.txt"
  mpirun -np "$NTASKS" ./abyss.exe -c "$TEST_CONFIG" >"$RUN_STDOUT" 2>"$RUN_STDERR" || true
else
  echo "[run] launching single-process" | tee -a "$RUN_DIR/run_meta.txt"
  ./abyss.exe -c "$TEST_CONFIG" >"$RUN_STDOUT" 2>"$RUN_STDERR" || true
fi

# Snapshot output directory
if [[ -d output ]]; then
  find output -maxdepth 2 -type f | sort >"$RUN_DIR/output_files.txt" || true
fi

echo "[run] done" | tee -a "$RUN_DIR/run_meta.txt"
