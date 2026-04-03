#!/usr/bin/env bash
set -euo pipefail

# Run McLuster IC generation as a standalone step.
# Called by SLURM job before ABYSS simulation.

RUN_DIR=${1:?run_dir}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

workflow_load_config
workflow_setup_env

REPO_ROOT="$(workflow_repo_root)"

TEST_ABS="$REPO_ROOT/$TEST_DIR"
CONFIG_ABS="$TEST_ABS/$RUN_CONFIG"

[[ -d "$TEST_ABS" ]] || workflow_die "TEST_DIR not found: $TEST_ABS"
[[ -f "$CONFIG_ABS" ]] || workflow_die "RUN_CONFIG not found: $CONFIG_ABS"

WORK_DIR="$RUN_DIR/work"
LOG="$RUN_DIR/mcluster.log"

{
  echo "== McLuster IC Generation =="
  echo "date=$(date)"
  echo "work_dir=$WORK_DIR"
  echo "OMP_NUM_THREADS=${OMP_NUM_THREADS:-not set}"
  echo ""
} > "$LOG"

# Create and stage work directory
mkdir -p "$WORK_DIR"
rm -rf "$WORK_DIR"/*
cp -a "$TEST_ABS/." "$WORK_DIR/"

# Stage McLuster binary
MCLUSTER_SRC="$REPO_ROOT/mcluster/mcluster_sse"
if [[ ! -x "$MCLUSTER_SRC" ]]; then
  workflow_die "McLuster binary not found: $MCLUSTER_SRC (build with USE_MCLUSTER=1)"
fi
cp "$MCLUSTER_SRC" "$WORK_DIR/mcluster"

MCLUSTER_BIN="$WORK_DIR/mcluster"

# Parse mcluster parameters from config file
CONFIG_FILE="$WORK_DIR/$RUN_CONFIG"

# Extract mcluster parameters using a simple parser
parse_mcluster_param() {
  local param=$1
  local default=$2
  grep -A50 '^\[mcluster\]' "$CONFIG_FILE" | grep -E "^${param}\s*=" | head -1 | sed 's/.*=\s*//' | tr -d ' ' || echo "$default"
}

N=$(parse_mcluster_param "N" "0")
M=$(parse_mcluster_param "M" "0")
P=$(parse_mcluster_param "P" "0")
R=$(parse_mcluster_param "R" "0.8")
f=$(parse_mcluster_param "f" "1")
Z=$(parse_mcluster_param "Z" "0.02")
b=$(parse_mcluster_param "b" "0")
e=$(parse_mcluster_param "e" "0")

# Build McLuster command
MCLUSTER_ARGS=()

# M takes precedence over N
if [[ "$M" != "0" && "$M" != "0.0" ]]; then
  MCLUSTER_ARGS+=("-M" "$M")
else
  MCLUSTER_ARGS+=("-N" "$N")
fi

MCLUSTER_ARGS+=("-P" "$P")
MCLUSTER_ARGS+=("-R" "$R")
MCLUSTER_ARGS+=("-f" "$f")
MCLUSTER_ARGS+=("-Z" "$Z")
MCLUSTER_ARGS+=("-b" "$b")
MCLUSTER_ARGS+=("-e" "$e")
MCLUSTER_ARGS+=("-C" "3")   # ASCII table format
MCLUSTER_ARGS+=("-u" "1")   # Astrophysical units
MCLUSTER_ARGS+=("-o" "mcluster_ic")

echo "Running: $MCLUSTER_BIN ${MCLUSTER_ARGS[*]}" >> "$LOG"
echo "" >> "$LOG"

pushd "$WORK_DIR" >/dev/null

set +e
"$MCLUSTER_BIN" "${MCLUSTER_ARGS[@]}" >> "$LOG" 2>&1
rc=$?
set -e

popd >/dev/null

if [[ $rc -ne 0 ]]; then
  echo "McLuster failed with exit code: $rc" >> "$LOG"
  workflow_die "McLuster failed (see $LOG)"
fi

# Verify output file exists
MCLUSTER_OUTPUT="$WORK_DIR/mcluster_ic.txt"
if [[ ! -f "$MCLUSTER_OUTPUT" ]]; then
  workflow_die "McLuster output not found: $MCLUSTER_OUTPUT"
fi

# Transform McLuster output to ABYSS format
# McLuster format: mass x y z vx vy vz [extras...]
# ABYSS format:    x y z vx vy vz mass
# Unit conversions: pc -> kpc (divide by 1000), Msun -> 1e-9 Msun (divide by 1e9)

ABYSS_IC="$WORK_DIR/mcluster_abyss.dat"

echo "" >> "$LOG"
echo "Transforming McLuster output to ABYSS format..." >> "$LOG"

awk '
BEGIN { count = 0 }
/^#/ { next }  # Skip header lines
NF >= 7 {
  mass = $1
  x = $2 / 1000.0    # pc -> kpc
  y = $3 / 1000.0
  z = $4 / 1000.0
  vx = $5            # km/s unchanged
  vy = $6
  vz = $7
  mass_out = mass / 1e9  # Msun -> 1e-9 Msun units
  printf "%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", x, y, z, vx, vy, vz, mass_out
  count++
}
END { print "Transformed " count " particles" > "/dev/stderr" }
' "$MCLUSTER_OUTPUT" > "$ABYSS_IC" 2>> "$LOG"

if [[ ! -s "$ABYSS_IC" ]]; then
  workflow_die "Failed to create ABYSS IC file: $ABYSS_IC"
fi

NPARTICLES=$(wc -l < "$ABYSS_IC")
echo "Created ABYSS IC file: $ABYSS_IC ($NPARTICLES particles)" >> "$LOG"
echo "" >> "$LOG"
echo "McLuster IC generation complete" >> "$LOG"
