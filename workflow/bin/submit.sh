#!/usr/bin/env bash
set -euo pipefail

# Submit or run workflow.
# Usage examples:
#   workflow/bin/submit.sh --scheduler slurm --tag test1
#   workflow/bin/submit.sh --scheduler local --tag test1
# Optional overrides:
#   --config test/test1/config.toml
#   --test-dir test/test1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

workflow_load_config

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
    -h|--help)
      sed -n '1,120p' "$0"
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
export WF_NTASKS_OVERRIDE="$NTASKS"
export WF_GPUS_OVERRIDE="$GPUS"
export WF_PYTHON_OVERRIDE="${PYTHON:-}"
export WF_RUN_COMPILE_OVERRIDE="$RUN_COMPILE"
export WF_RUN_RUN_OVERRIDE="$RUN_RUN"
export WF_RUN_ANALYZE_OVERRIDE="$RUN_ANALYZE"
export WF_SUMMARY_FILE_OVERRIDE="${SUMMARY_FILE:-}"
export WF_SUMMARY_STACK_FILE_OVERRIDE="${SUMMARY_STACK_FILE:-}"

REPO_ROOT="$(workflow_repo_root)"
RUN_DIR="$(workflow_run_dir "$TAG")"
mkdir -p "$RUN_DIR"

# Capture resolved settings for reproducibility
{
  echo "timestamp=$(workflow_timestamp)"
  echo "repo_root=$REPO_ROOT"
  echo "scheduler=$SCHEDULER"
  echo "test_dir=$TEST_DIR"
  echo "run_config=$RUN_CONFIG"
  echo "use_cuda=$USE_CUDA"
  echo "use_sevn=$USE_SEVN"
  echo "nodes=$NODES"
  echo "ntasks=$NTASKS"
  echo "gpus=$GPUS"
  echo "python=${PYTHON:-}"
  echo "summary_file=${SUMMARY_FILE:-}"
  echo "summary_stack_file=${SUMMARY_STACK_FILE:-}"
  echo "run_compile=$RUN_COMPILE"
  echo "run_run=$RUN_RUN"
  echo "run_analyze=$RUN_ANALYZE"
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
    TEMPLATE="$REPO_ROOT/workflow/templates/slurm.sbatch.in"
    JOB_NAME="abyss_${TAG}"

    JOB_SH="$RUN_DIR/job.sbatch"
    "$REPO_ROOT/workflow/bin/render.sh" "$TEMPLATE" "$JOB_SH" \
      ACCOUNT "$ACCOUNT" \
      PARTITION "$PARTITION" \
      WALLTIME "$WALLTIME" \
      NODES "$NODES" \
      NTASKS "$NTASKS" \
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
    sbatch "$JOB_SH" | tee "$RUN_DIR/submit.log"
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
