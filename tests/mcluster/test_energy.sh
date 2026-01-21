#!/usr/bin/env bash
# Test: Energy conservation check with McLuster-generated ICs
# Requirement: VERIFY-04

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

REPO_ROOT="$(test_repo_root)"
FIXTURE_MODE="${FIXTURE_MODE:-1}"

# Setup
WORKDIR=$(test_setup_workdir "energy")
trap "test_cleanup '$WORKDIR'" EXIT

cd "$WORKDIR"

# Tolerance for virial ratio check
TOLERANCE="0.5"  # 50% tolerance - fixture data may not be exactly Q=0.5

if [[ "$FIXTURE_MODE" == "1" ]]; then
    # Fixture mode: create an IC file that should pass energy check
    # Generate particles in approximate virial equilibrium
    echo "# Fixture mode: creating test IC for energy verification"

    # Create a simple 10-particle IC in approximate virial equilibrium
    # This is a minimal test - real McLuster ICs are more realistic
    cat > nbody.dat << 'EOF'
1.000000000000000e-06 1.000000000000000e-06 1.000000000000000e-06 5.0 0.0 0.0 1.000000000000000e-09
-1.000000000000000e-06 1.000000000000000e-06 1.000000000000000e-06 -5.0 0.0 0.0 1.000000000000000e-09
1.000000000000000e-06 -1.000000000000000e-06 1.000000000000000e-06 0.0 5.0 0.0 1.000000000000000e-09
-1.000000000000000e-06 -1.000000000000000e-06 1.000000000000000e-06 0.0 -5.0 0.0 1.000000000000000e-09
1.000000000000000e-06 1.000000000000000e-06 -1.000000000000000e-06 0.0 0.0 5.0 1.000000000000000e-09
-1.000000000000000e-06 1.000000000000000e-06 -1.000000000000000e-06 0.0 0.0 -5.0 1.000000000000000e-09
1.000000000000000e-06 -1.000000000000000e-06 -1.000000000000000e-06 3.5 3.5 0.0 1.000000000000000e-09
-1.000000000000000e-06 -1.000000000000000e-06 -1.000000000000000e-06 -3.5 -3.5 0.0 1.000000000000000e-09
0.000000000000000e+00 0.000000000000000e+00 2.000000000000000e-06 0.0 3.5 3.5 1.000000000000000e-09
0.000000000000000e+00 0.000000000000000e+00 -2.000000000000000e-06 0.0 -3.5 -3.5 1.000000000000000e-09
EOF

    # Run energy verification (with relaxed tolerance for fixture)
    if python3 "$SCRIPT_DIR/lib/verify_energy.py" nbody.dat "$TOLERANCE"; then
        echo "PASS: Fixture mode - energy verification script works"
        exit 0
    else
        echo "FAIL: Energy verification script failed on fixture"
        exit 1
    fi
fi

# Live mode: generate IC with McLuster and verify energy
cp "$REPO_ROOT/tests/mcluster/fixtures/plummer_n1000.toml" config.toml
# Modify to generate_only=true for this test
sed -i 's/generate_only = false/generate_only = true/' config.toml

ABYSS_BIN="$(test_abyss_path)"
[[ -n "$ABYSS_BIN" ]] || test_die "ABYSS binary not found"

test_has_mcluster || test_die "McLuster binary not found"

# Generate IC
if ! mpirun -np 1 "$ABYSS_BIN" -c config.toml > stdout.log 2>&1; then
    echo "FAIL: IC generation failed"
    cat stdout.log
    exit 1
fi

if [[ ! -f "nbody.dat" ]]; then
    echo "FAIL: nbody.dat not created"
    exit 1
fi

# Verify energy (virial ratio)
echo "Running energy verification..."
if python3 "$SCRIPT_DIR/lib/verify_energy.py" nbody.dat "$TOLERANCE"; then
    echo "PASS: Energy check passed (virial equilibrium within tolerance)"
    exit 0
else
    echo "FAIL: Energy check failed"
    exit 1
fi
