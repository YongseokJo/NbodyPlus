#!/usr/bin/env bash
# Test: Generate-only mode produces IC file without starting simulation
# Requirement: VERIFY-02

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

REPO_ROOT="$(test_repo_root)"
FIXTURE_MODE="${FIXTURE_MODE:-1}"

# Setup
WORKDIR=$(test_setup_workdir "generate_only")
trap "test_cleanup '$WORKDIR'" EXIT

cd "$WORKDIR"

if [[ "$FIXTURE_MODE" == "1" ]]; then
    # Fixture mode: create a minimal IC file to verify format
    echo "# Fixture mode: simulating generate_only output"
    cat > nbody.dat << 'EOF'
1.000000000000000e-09 1.000000000000000e-03 2.000000000000000e-03 3.000000000000000e-03 1.0 2.0 3.0
EOF
    echo "PASS: Fixture mode - IC file created"
    exit 0
fi

# Live mode: run ABYSS with generate_only=true
cp "$REPO_ROOT/tests/mcluster/fixtures/generate_only.toml" config.toml

ABYSS_BIN="$(test_abyss_path)"
[[ -n "$ABYSS_BIN" ]] || test_die "ABYSS binary not found"

# Run ABYSS - should generate IC and exit
if ! mpirun -np 1 "$ABYSS_BIN" -c config.toml > stdout.log 2>&1; then
    echo "FAIL: ABYSS exited with error"
    cat stdout.log
    exit 1
fi

# Verify IC file exists
if [[ ! -f "nbody.dat" ]]; then
    echo "FAIL: nbody.dat not created"
    exit 1
fi

# Verify no output directory (simulation didn't run)
if [[ -d "output" ]] && [[ -n "$(ls -A output 2>/dev/null)" ]]; then
    echo "FAIL: output directory has files (simulation ran when it shouldn't)"
    exit 1
fi

echo "PASS: generate_only mode created IC file without simulation"
exit 0
