#!/usr/bin/env bash
set -euo pipefail

ABYSS_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

load_config() {
  # Avoid surprises from inherited environment variables (e.g. NTASKS) by clearing
  # the workflow variables before loading config defaults and local overrides.
  unset SYSTEM_NAME SCHEDULER ACCOUNT PARTITION QUEUE
  unset NODES NTASKS CPUS_PER_TASK GPUS WALLTIME
  unset TEST_DIR TEST_CONFIG
  unset USE_CUDA CXX MPI_LAUNCHER SRUN_MPI OMPI_BACKING_DIR WORK_DIR
  unset MODULES

  # shellcheck disable=SC1091
  source "$ABYSS_ROOT/workflow/config.sh"
  if [[ -f "$ABYSS_ROOT/workflow/config.local.sh" ]]; then
    # shellcheck disable=SC1091
    source "$ABYSS_ROOT/workflow/config.local.sh"
  fi
}

maybe_load_modules() {
  if ! command -v module >/dev/null 2>&1; then
    return 0
  fi

  # MODULES may be empty
  if [[ ${#MODULES[@]} -eq 0 ]]; then
    return 0
  fi

  module purge || true
  for m in "${MODULES[@]}"; do
    module load "$m" || {
      echo "[workflow] WARNING: failed to load module: $m" >&2
    }
  done
}

detect_use_cuda() {
  # Echo 0/1
  case "${USE_CUDA}" in
    0|1) echo "${USE_CUDA}"; return 0 ;;
    auto)
      if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
        echo 1
      else
        echo 0
      fi
      return 0
      ;;
    *)
      echo "[workflow] Invalid USE_CUDA=${USE_CUDA} (expected 0|1|auto)" >&2
      exit 2
      ;;
  esac
}

timestamp_tag() {
  date +"%Y%m%d_%H%M%S"
}
