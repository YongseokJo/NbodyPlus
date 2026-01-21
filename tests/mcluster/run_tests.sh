#!/usr/bin/env bash
# McLuster Integration Test Runner
#
# Runs all McLuster integration tests and reports results.
#
# Environment variables:
#   FIXTURE_MODE=1  Use pre-generated fixtures (default, no gfortran needed)
#   FIXTURE_MODE=0  Run live McLuster execution (requires gfortran)
#   FAIL_FAST=1     Stop on first failure
#
# Exit codes:
#   0 = All tests passed
#   1 = One or more tests failed

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

# Test list in execution order
TESTS=(
    "test_generate_only"
    "test_run_only"
    "test_generate_run"
    "test_energy"
)

# Configuration
FIXTURE_MODE="${FIXTURE_MODE:-1}"
FAIL_FAST="${FAIL_FAST:-0}"
OUTPUT_DIR="$SCRIPT_DIR/output"

# Counters
PASSED=0
FAILED=0
SKIPPED=0

# Setup
mkdir -p "$OUTPUT_DIR"

echo ""
echo "McLuster Integration Tests"
echo "=========================="
echo "Mode: $([ "$FIXTURE_MODE" == "1" ] && echo "fixture" || echo "live")"
echo "Fail fast: $([ "$FAIL_FAST" == "1" ] && echo "yes" || echo "no")"
echo ""

# Run tests
for test in "${TESTS[@]}"; do
    test_script="$SCRIPT_DIR/${test}.sh"
    log_file="$OUTPUT_DIR/${test}.log"

    if [[ ! -x "$test_script" ]]; then
        printf "%-35s SKIP (not found)\n" "$test"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Run test with environment
    printf "%-35s " "$test..."
    if FIXTURE_MODE="$FIXTURE_MODE" "$test_script" > "$log_file" 2>&1; then
        echo "PASS"
        PASSED=$((PASSED + 1))
    else
        echo "FAIL"
        FAILED=$((FAILED + 1))

        # Show error output
        echo "  --- Error output ---"
        tail -15 "$log_file" | sed 's/^/  /'
        echo "  --- End output ---"

        if [[ "$FAIL_FAST" == "1" ]]; then
            echo ""
            echo "Stopping (FAIL_FAST=1)"
            break
        fi
    fi
done

# Summary
echo ""
echo "=========================="
echo "Results: $PASSED passed, $FAILED failed, $SKIPPED skipped"

if [[ $FAILED -gt 0 ]]; then
    echo ""
    echo "Failed test logs in: $OUTPUT_DIR/"
    exit 1
fi

if [[ $SKIPPED -gt 0 ]] && [[ $PASSED -eq 0 ]]; then
    echo "WARNING: All tests skipped"
    exit 1
fi

echo "All tests passed!"
exit 0
