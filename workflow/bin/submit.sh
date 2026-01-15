#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

load_config

JOB_KIND=${1:-"run"}

JOB_SCRIPT=$(bash "$ABYSS_ROOT/workflow/bin/render_job.sh" "$JOB_KIND")
echo "[workflow] Rendered: $JOB_SCRIPT"

case "$SCHEDULER" in
  slurm)
    if ! command -v sbatch >/dev/null 2>&1; then
      echo "[workflow] sbatch not found" >&2
      exit 2
    fi
    sbatch "$JOB_SCRIPT"
    ;;
  pbs)
    if ! command -v qsub >/dev/null 2>&1; then
      echo "[workflow] qsub not found" >&2
      exit 2
    fi
    qsub "$JOB_SCRIPT"
    ;;
  *)
    echo "[workflow] Unsupported SCHEDULER=$SCHEDULER" >&2
    exit 2
    ;;
esac
