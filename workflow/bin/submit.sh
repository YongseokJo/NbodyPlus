#!/usr/bin/env bash
set -euo pipefail

# Submit or run workflow.
# Run `workflow/bin/submit.sh --help` for full usage.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

workflow_load_config

print_help() {
  cat <<'EOF'
ABYSS workflow submit helper

Creates a timestamped run directory under workflow/runs/, then:
  - local: runs compile/run/analyze directly
  - slurm/pbs: writes a job script into the run dir and submits it

Usage:
  workflow/bin/submit.sh [options]

Common examples:
  # Slurm submission (cluster)
  workflow/bin/submit.sh --scheduler slurm --tag run

  # Local run (single rank by default)
  workflow/bin/submit.sh --scheduler local --tag smoke

  # Override test directory + config file
  workflow/bin/submit.sh --test-dir tests/test_10 --config config_smoke.txt --tag smoke

  # Skip steps
  workflow/bin/submit.sh --skip-compile
  workflow/bin/submit.sh --skip-run
  workflow/bin/submit.sh --skip-analyze

  # Enable profiling (adds -DPERFORMANCETRACE)
  workflow/bin/submit.sh --profile

Options:
  --tag <name>               Prefix for run directory name (default: run)
  --scheduler <slurm|pbs|local>
                             Override scheduler from workflow/config.sh
  --test-dir <path>          Test directory relative to repo root (default from config)
  --config <path>            Config file relative to --test-dir (default from config)
  --ntasks <N>               MPI ranks (default from config; local defaults to 1)
  --python <path>            Python interpreter used by analyze/tools
  --summary-file <name>      Summary file name written into the run dir (default: summary.txt)
  --summary-stack-file <path>
                             Stacked summary file path (default: summary_runs.tsv)
  --profile                  Enable profiling (adds -DPERFORMANCETRACE to build)
  --mcluster                 Enable McLuster IC generator build (requires gfortran)
  --no-mcluster              Disable McLuster build
  --skip-compile             Skip the compile step
  --skip-run                 Skip the run step
  --skip-analyze             Skip the analyze step

Stacking / summary_runs.tsv:
  By default the workflow appends ONLY the current run to the stack file.
  To force a full rebuild across all existing runs:
    --stack-rebuild-all

Hardware strings (cpu_arch/gpu_arch):
  These are detected at analyze-time when possible (preferably on compute nodes).
  If no GPU is available/detectable, the stacked GPU column will be 'none'.

Exit codes:
  0 on success; non-zero on argument errors or local run failures.
EOF
}

TAG="run"
SCHED_OVERRIDE=""
TEST_DIR_OVERRIDE=""
RUN_CONFIG_OVERRIDE=""
NTASKS_OVERRIDE=""
PYTHON_OVERRIDE=""
SUMMARY_FILE_OVERRIDE=""
SUMMARY_STACK_FILE_OVERRIDE=""
SKIP_COMPILE=0
SKIP_RUN=0
SKIP_ANALYZE=0
STACK_REBUILD_ALL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tag)
      TAG="$2"; shift 2 ;;
    --scheduler)
      SCHED_OVERRIDE="$2"; shift 2 ;;
    --test-dir)
      TEST_DIR_OVERRIDE="$2"; shift 2 ;;
    --config)
      RUN_CONFIG_OVERRIDE="$2"; shift 2 ;;
    --ntasks)
      NTASKS_OVERRIDE="$2"; shift 2 ;;
    --python)
      PYTHON_OVERRIDE="$2"; shift 2 ;;
    --summary-file)
      SUMMARY_FILE_OVERRIDE="$2"; shift 2 ;;
    --summary-stack-file)
      SUMMARY_STACK_FILE_OVERRIDE="$2"; shift 2 ;;
    --skip-compile)
      SKIP_COMPILE=1; shift 1 ;;
    --skip-run)
      SKIP_RUN=1; shift 1 ;;
    --skip-analyze)
      SKIP_ANALYZE=1; shift 1 ;;
    --profile|--profiling)
      ENABLE_PROFILING=1; shift 1 ;;
    --mcluster)
      USE_MCLUSTER=1; shift 1 ;;
    --no-mcluster)
      USE_MCLUSTER=0; shift 1 ;;
    --stack-rebuild-all)
      STACK_REBUILD_ALL=1; shift 1 ;;
    -h|--help)
      print_help
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

if [[ -n "$SCHED_OVERRIDE" ]]; then
  SCHEDULER="$SCHED_OVERRIDE"
fi
if [[ -n "$TEST_DIR_OVERRIDE" ]]; then
  TEST_DIR="$TEST_DIR_OVERRIDE"
fi
if [[ -n "$RUN_CONFIG_OVERRIDE" ]]; then
  RUN_CONFIG="$RUN_CONFIG_OVERRIDE"
fi
if [[ -n "$NTASKS_OVERRIDE" ]]; then
  NTASKS="$NTASKS_OVERRIDE"
fi
if [[ -n "$PYTHON_OVERRIDE" ]]; then
  PYTHON="$PYTHON_OVERRIDE"
fi
if [[ -n "$SUMMARY_FILE_OVERRIDE" ]]; then
  SUMMARY_FILE="$SUMMARY_FILE_OVERRIDE"
fi
if [[ -n "$SUMMARY_STACK_FILE_OVERRIDE" ]]; then
  SUMMARY_STACK_FILE="$SUMMARY_STACK_FILE_OVERRIDE"
fi
if [[ "$SKIP_COMPILE" -eq 1 ]]; then
  RUN_COMPILE=0
else
  RUN_COMPILE=1
fi
if [[ "$SKIP_RUN" -eq 1 ]]; then
  RUN_RUN=0
else
  RUN_RUN=1
fi
if [[ "$SKIP_ANALYZE" -eq 1 ]]; then
  RUN_ANALYZE=0
else
  RUN_ANALYZE=1
fi

# Keep local runs lightweight by default.
if [[ "$SCHEDULER" == "local" && -z "$NTASKS_OVERRIDE" ]]; then
  NTASKS="1"
fi

# Local runs often execute on login nodes without GPUs.
# If no GPU is detected, default to a CPU-only build/run to avoid failing
# during CUDA initialization.
if [[ "$SCHEDULER" == "local" && "${USE_CUDA:-0}" == "1" ]]; then
  has_gpu=0
  if command -v nvidia-smi &>/dev/null && nvidia-smi -L &>/dev/null; then
    has_gpu=1
  elif [[ -e /dev/nvidia0 ]]; then
    has_gpu=1
  fi

  if [[ $has_gpu -eq 0 ]]; then
    USE_CUDA="0"
  fi
fi

# Propagate resolved config to subprocesses.
# (Without this, compile.sh/run.sh/analyze.sh would reload config defaults.)
export WF_SCHEDULER_OVERRIDE="$SCHEDULER"
export WF_TEST_DIR_OVERRIDE="$TEST_DIR"
export WF_RUN_CONFIG_OVERRIDE="$RUN_CONFIG"
export WF_USE_CUDA_OVERRIDE="$USE_CUDA"
export WF_USE_SEVN_OVERRIDE="$USE_SEVN"
export WF_USE_MCLUSTER_OVERRIDE="${USE_MCLUSTER:-1}"
export WF_ENABLE_PROFILING_OVERRIDE="${ENABLE_PROFILING:-0}"
export WF_NTASKS_OVERRIDE="$NTASKS"
export WF_GPUS_OVERRIDE="$GPUS"
export WF_PYTHON_OVERRIDE="${PYTHON:-}"
export WF_RUN_COMPILE_OVERRIDE="$RUN_COMPILE"
export WF_RUN_RUN_OVERRIDE="$RUN_RUN"
export WF_RUN_ANALYZE_OVERRIDE="$RUN_ANALYZE"
export WF_SUMMARY_FILE_OVERRIDE="${SUMMARY_FILE:-}"
export WF_SUMMARY_STACK_FILE_OVERRIDE="${SUMMARY_STACK_FILE:-}"

# By default, keep stacking fast and incremental (only the current run).
# Users can opt in to a full rebuild across all runs.
export WF_STACK_REBUILD_ALL="$STACK_REBUILD_ALL"

REPO_ROOT="$(workflow_repo_root)"
RUN_DIR="$(workflow_run_dir "$TAG")"
mkdir -p "$RUN_DIR"

# Capture VCS information at submit time.
# (CPU/GPU architecture is detected at runtime during analyze; submit-time values
# can reflect the login node for batch schedulers.)

git_commit=""
git_commit_long=""
git_branch=""
git_tag=""
if command -v git >/dev/null 2>&1; then
  git_commit="$(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || true)"
  git_commit_long="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || true)"
  git_branch="$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || true)"
  git_tag="$(git -C "$REPO_ROOT" describe --tags --exact-match 2>/dev/null || true)"
fi
if [[ -z "$git_tag" && -n "$git_commit" ]]; then
  git_tag="$git_commit"
fi

# Capture resolved settings for reproducibility
{
  echo "timestamp=$(workflow_timestamp)"
  echo "repo_root=$REPO_ROOT"
  echo "scheduler=$SCHEDULER"
  echo "test_dir=$TEST_DIR"
  echo "run_config=$RUN_CONFIG"
  echo "use_cuda=$USE_CUDA"
  echo "use_sevn=$USE_SEVN"
  echo "use_mcluster=${USE_MCLUSTER:-1}"
  echo "enable_profiling=${ENABLE_PROFILING:-0}"
  echo "nodes=$NODES"
  echo "ntasks=$NTASKS"
  echo "gpus=$GPUS"
  echo "python=${PYTHON:-}"
  echo "summary_file=${SUMMARY_FILE:-}"
  echo "summary_stack_file=${SUMMARY_STACK_FILE:-}"
  echo "run_compile=$RUN_COMPILE"
  echo "run_run=$RUN_RUN"
  echo "run_analyze=$RUN_ANALYZE"
  echo "git_commit=$git_commit"
  echo "git_commit_long=$git_commit_long"
  echo "git_branch=$git_branch"
  echo "git_tag=$git_tag"
} > "$RUN_DIR/meta.txt"

cp "$REPO_ROOT/workflow/config.sh" "$RUN_DIR/config.sh"
if [[ -f "$REPO_ROOT/workflow/config.local.sh" ]]; then
  cp "$REPO_ROOT/workflow/config.local.sh" "$RUN_DIR/config.local.sh"
fi

case "$SCHEDULER" in
  local)
    if [[ "$RUN_COMPILE" -eq 1 ]]; then
      "$REPO_ROOT/workflow/bin/compile.sh" "$RUN_DIR"
    else
      echo "Skipping compile (RUN_COMPILE=0)"
    fi

    run_rc=0
    if [[ "$RUN_RUN" -eq 1 ]]; then
      set +e
      "$REPO_ROOT/workflow/bin/run.sh" "$RUN_DIR"
      run_rc=$?
      set -e
    else
      echo "Skipping run (RUN_RUN=0)"
    fi

    if [[ "$RUN_ANALYZE" -eq 1 ]]; then
      # Always analyze (even if run failed) to capture artifacts/log summaries.
      "$REPO_ROOT/workflow/bin/analyze.sh" "$RUN_DIR" || true
    else
      echo "Skipping analyze (RUN_ANALYZE=0)"
    fi

    if [[ $run_rc -ne 0 ]]; then
      workflow_die "Run failed (analyzed artifacts; see $RUN_DIR)"
    fi

    echo "Run complete: $RUN_DIR"
    ;;

  slurm)
    # Check if config has [mcluster] section - requires two-job submission
    CONFIG_PATH="$REPO_ROOT/$TEST_DIR/$RUN_CONFIG"
    NEEDS_MCLUSTER=0
    if [[ -f "$CONFIG_PATH" ]] && grep -q '^\[mcluster\]' "$CONFIG_PATH" 2>/dev/null; then
      NEEDS_MCLUSTER=1
    fi

    MCLUSTER_JOB_ID=""

    if [[ "$NEEDS_MCLUSTER" -eq 1 && "$RUN_RUN" -eq 1 ]]; then
      # Submit McLuster IC generation job first
      MCLUSTER_TEMPLATE="$REPO_ROOT/workflow/templates/mcluster.sbatch.in"
      MCLUSTER_JOB_NAME="mcl_${TAG}"
      MCLUSTER_JOB_SH="$RUN_DIR/mcluster.sbatch"

      "$REPO_ROOT/workflow/bin/render.sh" "$MCLUSTER_TEMPLATE" "$MCLUSTER_JOB_SH" \
        ACCOUNT "$ACCOUNT" \
        MCLUSTER_PARTITION "${MCLUSTER_PARTITION:-ciera-std}" \
        MCLUSTER_WALLTIME "${MCLUSTER_WALLTIME:-02:00:00}" \
        MCLUSTER_CPUS "${MCLUSTER_CPUS:-16}" \
        MCLUSTER_MEM "${MCLUSTER_MEM:-32G}" \
        RUN_COMPILE "$RUN_COMPILE" \
        JOB_NAME "$MCLUSTER_JOB_NAME" \
        REPO_ROOT "$REPO_ROOT" \
        RUN_DIR "$RUN_DIR"

      chmod +x "$MCLUSTER_JOB_SH"

      echo "Submitting McLuster job: $MCLUSTER_JOB_SH"
      MCLUSTER_SUBMIT_OUTPUT=$(sbatch "$MCLUSTER_JOB_SH")
      echo "$MCLUSTER_SUBMIT_OUTPUT" | tee -a "$RUN_DIR/submit.log"

      # Extract job ID from "Submitted batch job 12345"
      MCLUSTER_JOB_ID=$(echo "$MCLUSTER_SUBMIT_OUTPUT" | grep -oP 'Submitted batch job \K\d+')
      echo "McLuster job ID: $MCLUSTER_JOB_ID" >> "$RUN_DIR/submit.log"
    fi

    # Submit main ABYSS simulation job
    TEMPLATE="$REPO_ROOT/workflow/templates/slurm.sbatch.in"
    JOB_NAME="abyss_${TAG}"
    JOB_SH="$RUN_DIR/job.sbatch"

    # If McLuster job was submitted, ABYSS job skips compile (already done)
    ABYSS_RUN_COMPILE="$RUN_COMPILE"
    if [[ -n "$MCLUSTER_JOB_ID" ]]; then
      ABYSS_RUN_COMPILE=0
    fi

    "$REPO_ROOT/workflow/bin/render.sh" "$TEMPLATE" "$JOB_SH" \
      ACCOUNT "$ACCOUNT" \
      PARTITION "$PARTITION" \
      WALLTIME "$WALLTIME" \
      NODES "$NODES" \
      NTASKS "$NTASKS" \
      CPUS_PER_TASK "$CPUS_PER_TASK" \
      GPUS "$GPUS" \
      RUN_COMPILE "$ABYSS_RUN_COMPILE" \
      RUN_RUN "$RUN_RUN" \
      RUN_ANALYZE "$RUN_ANALYZE" \
      JOB_NAME "$JOB_NAME" \
      REPO_ROOT "$REPO_ROOT" \
      RUN_DIR "$RUN_DIR"

    chmod +x "$JOB_SH"

    # Submit with dependency if McLuster job was submitted
    echo "Submitting ABYSS job: $JOB_SH"
    if [[ -n "$MCLUSTER_JOB_ID" ]]; then
      echo "  (depends on McLuster job $MCLUSTER_JOB_ID)"
      sbatch --dependency=afterok:"$MCLUSTER_JOB_ID" "$JOB_SH" | tee -a "$RUN_DIR/submit.log"
    else
      sbatch "$JOB_SH" | tee -a "$RUN_DIR/submit.log"
    fi
    ;;

  pbs)
    TEMPLATE="$REPO_ROOT/workflow/templates/pbs.pbs.in"
    JOB_NAME="abyss_${TAG}"

    JOB_SH="$RUN_DIR/job.pbs"
    "$REPO_ROOT/workflow/bin/render.sh" "$TEMPLATE" "$JOB_SH" \
      ACCOUNT "$ACCOUNT" \
      PARTITION "$PARTITION" \
      WALLTIME "$WALLTIME" \
      NODES "$NODES" \
      CPUS_PER_TASK "$CPUS_PER_TASK" \
      GPUS "$GPUS" \
      RUN_COMPILE "$RUN_COMPILE" \
      RUN_RUN "$RUN_RUN" \
      RUN_ANALYZE "$RUN_ANALYZE" \
      JOB_NAME "$JOB_NAME" \
      REPO_ROOT "$REPO_ROOT" \
      RUN_DIR "$RUN_DIR"

    chmod +x "$JOB_SH"

    echo "Submitting: $JOB_SH"
    qsub "$JOB_SH" | tee "$RUN_DIR/submit.log"
    ;;

  *)
    workflow_die "Unknown scheduler: $SCHEDULER"
    ;;
esac
