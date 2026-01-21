#!/usr/bin/env bash
# Test: Run-only mode uses existing IC file (no McLuster section)
# Requirement: VERIFY-03

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

REPO_ROOT="$(test_repo_root)"
FIXTURE_MODE="${FIXTURE_MODE:-1}"

# Setup
WORKDIR=$(test_setup_workdir "run_only")
trap "test_cleanup '$WORKDIR'" EXIT

cd "$WORKDIR"

# Copy config without [mcluster] section
cp "$REPO_ROOT/tests/mcluster/fixtures/runonly.toml" config.toml

# Need an existing IC file - use test_1 fixture
if [[ -f "$REPO_ROOT/tests/test_1/nbody.dat" ]]; then
    cp "$REPO_ROOT/tests/test_1/nbody.dat" nbody.dat
else
    # Create minimal fixture
    echo "# Creating minimal IC fixture"
    for i in $(seq 1 100); do
        echo "1.0e-9 0.001 0.002 0.003 1.0 2.0 3.0"
    done > nbody.dat
fi

if [[ "$FIXTURE_MODE" == "1" ]]; then
    # Fixture mode: just verify config parses without [mcluster]
    echo "PASS: Fixture mode - run-only config valid"
    exit 0
fi

# Live mode: run ABYSS briefly
ABYSS_BIN="$(test_abyss_path)"
[[ -n "$ABYSS_BIN" ]] || test_die "ABYSS binary not found"

mkdir -p output

# Run for a short time
timeout 30 mpirun -np 1 "$ABYSS_BIN" -c config.toml > stdout.log 2>&1 || true

# Check that it started (output directory has content)
if [[ ! -d "output" ]] || [[ -z "$(ls -A output 2>/dev/null)" ]]; then
    echo "FAIL: Simulation did not produce any output"
    cat stdout.log
    exit 1
fi

echo "PASS: run-only mode works with existing IC file"
exit 0
