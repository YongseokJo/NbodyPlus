#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

load_config

JOB_KIND=${1:-"run"} # run|compile|custom (currently informational)

TAG=$(timestamp_tag)
JOB_NAME="abyss-${JOB_KIND}-${TAG}"

RUN_DIR="$ABYSS_ROOT/${WORK_DIR}/${TAG}"
mkdir -p "$RUN_DIR" "$ABYSS_ROOT/logs" "$ABYSS_ROOT/${WORK_DIR}"

LOG_OUT="$ABYSS_ROOT/logs/${JOB_NAME}.out"
LOG_ERR="$ABYSS_ROOT/logs/${JOB_NAME}.err"

TEMPLATE=""
OUT_SCRIPT="$RUN_DIR/job.sh"

case "$SCHEDULER" in
  slurm) TEMPLATE="$ABYSS_ROOT/workflow/templates/slurm.sbatch.in" ;;
  pbs)   TEMPLATE="$ABYSS_ROOT/workflow/templates/pbs.pbs.in" ;;
  *)
    echo "[workflow] render_job.sh supports SCHEDULER=slurm|pbs only (got $SCHEDULER)" >&2
    exit 2
    ;;
esac

if [[ ! -f "$TEMPLATE" ]]; then
  echo "[workflow] Missing template: $TEMPLATE" >&2
  exit 2
fi

PBS_QUEUE_LINE="#PBS -q ${QUEUE}"
if [[ -z "$QUEUE" ]]; then
  PBS_QUEUE_LINE=""
fi

CPUS_TOTAL=$(( CPUS_PER_TASK * NTASKS ))

# Simple placeholder substitution (no external deps)
content=$(cat "$TEMPLATE")
content=${content//\{\{JOB_NAME\}\}/$JOB_NAME}
content=${content//\{\{ACCOUNT\}\}/$ACCOUNT}
content=${content//\{\{PARTITION\}\}/$PARTITION}
content=${content//\{\{NODES\}\}/$NODES}
content=${content//\{\{NTASKS\}\}/$NTASKS}
content=${content//\{\{CPUS_PER_TASK\}\}/$CPUS_PER_TASK}
content=${content//\{\{CPUS_TOTAL\}\}/$CPUS_TOTAL}
content=${content//\{\{GPUS\}\}/$GPUS}
content=${content//\{\{WALLTIME\}\}/$WALLTIME}
content=${content//\{\{LOG_OUT\}\}/$LOG_OUT}
content=${content//\{\{LOG_ERR\}\}/$LOG_ERR}
content=${content//\{\{ABYSS_ROOT\}\}/$ABYSS_ROOT}
content=${content//\{\{RUN_DIR\}\}/$RUN_DIR}
content=${content//\{\{OMPI_BACKING_DIR\}\}/$OMPI_BACKING_DIR}
content=${content//\{\{PBS_QUEUE_LINE\}\}/$PBS_QUEUE_LINE}

printf "%s\n" "$content" > "$OUT_SCRIPT"
chmod +x "$OUT_SCRIPT" || true

echo "$OUT_SCRIPT"
