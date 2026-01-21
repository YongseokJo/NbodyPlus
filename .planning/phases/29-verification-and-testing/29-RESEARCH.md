# Phase 29: Verification and Testing - Research

**Researched:** 2026-01-20
**Domain:** Shell script testing, energy verification, McLuster integration testing
**Confidence:** HIGH

## Summary

Phase 29 implements comprehensive testing for the McLuster integration completed in Phases 25-28. The testing approach uses shell scripts with Make targets, following the existing workflow patterns in the ABYSS codebase. Tests verify three operational modes: full pipeline (generate+run), generate-only, and run-only (regression).

Energy verification focuses on initial conditions (IC) checking, not full simulation runs. The virial theorem provides the theoretical basis: for a system in virial equilibrium, 2K + U = 0 (where K is kinetic energy and U is potential energy). McLuster generates ICs in virial equilibrium by default (Q=0.5), so the ratio 2K/|U| should equal 1.0.

**Primary recommendation:** Implement fixture-based quick tests (pre-generated IC files) alongside live McLuster tests (full end-to-end). Use Python for energy computation (leverage existing analyze_energy.py patterns) called from shell scripts. Target 1e-4 relative tolerance on virial ratio check.

## Standard Stack

The established tools for this testing domain:

### Core
| Tool | Version | Purpose | Why Standard |
|------|---------|---------|--------------|
| Bash | 4.x+ | Test script logic | Universal, set -euo pipefail patterns |
| GNU Make | 3.8+ | Test invocation | Already used in ABYSS build system |
| Python | 3.x | Energy computation | Existing analyze_energy.py provides patterns |
| h5py | 3.x | HDF5 reading | Already used in ABYSS tools |
| NumPy | 1.x | Numerical computation | Already used in ABYSS tools |

### Supporting
| Tool | Version | Purpose | When to Use |
|------|---------|---------|-------------|
| bc | any | Floating-point comparison in bash | Virial ratio tolerance check |
| awk | any | Text parsing | Extract values from McLuster/ABYSS output |
| diff | any | File comparison | Compare IC file formats |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Shell scripts | Bats-core | More features but adds dependency |
| Python energy | Pure bash | Bash cannot do N^2 energy calculation efficiently |
| Make targets | npm/pytest | Already using Make, no need to change |

**No additional installation required** - all tools already present in the ABYSS environment.

## Architecture Patterns

### Recommended Test Structure
```
tests/
  mcluster/                      # McLuster integration tests
    fixtures/                    # Pre-generated test data
      plummer_n1000.toml        # Config for Plummer profile N=1000
      plummer_n1000.dat         # Pre-generated IC (fixture mode)
      king_n1000.toml           # Config for King profile N=1000
      king_n1000.dat            # Pre-generated IC
      kroupa_imf.toml           # Config with Kroupa IMF
      kroupa_imf.dat            # Pre-generated IC
    run_tests.sh                # Main test runner script
    test_generate_run.sh        # Test: config -> McLuster -> ABYSS (VERIFY-01)
    test_generate_only.sh       # Test: config -> McLuster -> IC file (VERIFY-02)
    test_run_only.sh            # Test: existing IC -> ABYSS (VERIFY-03)
    test_energy.sh              # Test: energy conservation check (VERIFY-04)
    lib/                        # Shared test utilities
      common.sh                 # Test helper functions
      verify_energy.py          # Energy computation script
```

### Pattern 1: Test Runner with Status Reporting
**What:** Central script that runs all tests with consistent output format
**When to use:** For the main test entry point
**Example:**
```bash
#!/usr/bin/env bash
# Source: Based on workflow/bin/common.sh patterns
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/common.sh"

TESTS=(
    "test_generate_only"
    "test_generate_run"
    "test_run_only"
    "test_energy"
)

FAILED=0
PASSED=0
FAIL_FAST="${FAIL_FAST:-0}"

for test in "${TESTS[@]}"; do
    printf "%-40s" "$test..."
    if "$SCRIPT_DIR/${test}.sh" > "$SCRIPT_DIR/output/${test}.log" 2>&1; then
        echo "PASS"
        ((PASSED++))
    else
        echo "FAIL"
        echo "  Error output:"
        tail -20 "$SCRIPT_DIR/output/${test}.log" | sed 's/^/    /'
        ((FAILED++))
        [[ "$FAIL_FAST" == "1" ]] && break
    fi
done

echo ""
echo "Results: $PASSED passed, $FAILED failed"
exit $((FAILED > 0 ? 1 : 0))
```

### Pattern 2: Fixture-Based vs Live McLuster Test Modes
**What:** Support both pre-generated fixtures (quick) and live McLuster execution (full)
**When to use:** For all tests that involve IC generation
**Example:**
```bash
#!/usr/bin/env bash
set -euo pipefail

# Mode selection: FIXTURE_MODE=1 uses pre-generated ICs, FIXTURE_MODE=0 runs McLuster
FIXTURE_MODE="${FIXTURE_MODE:-1}"

if [[ "$FIXTURE_MODE" == "1" ]]; then
    # Quick mode: use pre-generated fixture
    cp "$FIXTURES_DIR/plummer_n1000.dat" "$WORK_DIR/nbody.dat"
else
    # Full mode: run McLuster
    "$MCLUSTER_BINARY" -N 1000 -P 0 -R 0.8 -C 3 -u 1 -o "$WORK_DIR/mcluster_ic"
    # Transform to ABYSS format (using existing pipeline)
fi
```

### Pattern 3: Energy Verification
**What:** Compute kinetic and potential energy, check virial ratio
**When to use:** For VERIFY-04 (energy conservation check)
**Example (Python, based on analyze_energy.py):**
```python
#!/usr/bin/env python3
"""Verify energy of McLuster-generated initial conditions."""

import sys
import numpy as np

# Unit conversions (from src/def.h)
POSITION_UNIT_PC = 4.0
MASS_UNIT_MSUN = 0.0001424198

def compute_virial_ratio(filename):
    """Compute 2K/|U| for particles in IC file."""
    data = np.loadtxt(filename)
    # ABYSS format: x y z vx vy vz mass
    pos = data[:, 0:3]  # kpc in file, need to convert
    vel = data[:, 3:6]  # km/s
    mass = data[:, 6]   # 1e-9 Msun units

    # Convert to code units for G=1 calculation
    m_code = mass * 1e9 / MASS_UNIT_MSUN  # to Msun then code
    pos_code = pos * 1000 / POSITION_UNIT_PC  # kpc to pc to code

    # Kinetic: 0.5 * sum(m * v^2)
    # Note: velocity needs conversion km/s -> code units
    # ... (full implementation)

    kinetic = 0.5 * np.sum(m_code * np.sum(vel_code**2, axis=1))

    # Potential: -sum(mi*mj/rij) for i<j
    potential = 0.0
    n = len(mass)
    for i in range(n-1):
        for j in range(i+1, n):
            r = np.linalg.norm(pos_code[i] - pos_code[j])
            if r > 0:
                potential -= m_code[i] * m_code[j] / r

    virial_ratio = 2 * kinetic / abs(potential)
    return kinetic, potential, virial_ratio

if __name__ == "__main__":
    ic_file = sys.argv[1]
    tolerance = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4

    K, U, Q = compute_virial_ratio(ic_file)
    print(f"Kinetic:  {K:.10e}")
    print(f"Potential: {U:.10e}")
    print(f"Virial ratio (2K/|U|): {Q:.10f}")
    print(f"Expected: 1.0 (tolerance: {tolerance})")

    if abs(Q - 1.0) > tolerance:
        print(f"FAIL: |Q - 1.0| = {abs(Q-1.0):.6e} > {tolerance}")
        sys.exit(1)
    print("PASS: Energy check within tolerance")
    sys.exit(0)
```

### Anti-Patterns to Avoid
- **Hardcoded paths:** Use `$SCRIPT_DIR` and relative paths, not absolute
- **Missing error handling:** Always use `set -euo pipefail`
- **Silent failures:** Capture and report stderr on failure
- **Floating-point comparison in bash:** Use Python or `bc` for tolerance checks
- **Full simulation runs:** Tests verify IC generation only, not simulation

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Energy calculation | Custom bash math | Python numpy | O(N^2) computation, floating point |
| TOML generation | String concatenation | Here-doc templates | Proper escaping, readability |
| Test output format | Custom format | TAP-style (test name + PASS/FAIL) | Standard, easy to parse |
| Floating-point compare | `[ $a == $b ]` | bc or Python | Bash can't do float comparison |
| Path finding | Hardcoded `/path/to` | `dirname "${BASH_SOURCE[0]}"` | Portable across invocation methods |

**Key insight:** Bash is for orchestration and flow control; numerical computation and complex parsing should use Python.

## Common Pitfalls

### Pitfall 1: Floating-Point Tolerance in Bash
**What goes wrong:** Tests fail because bash `[ 0.9999 == 1.0 ]` is false
**Why it happens:** Bash doesn't do floating-point comparison
**How to avoid:** Use bc or Python for numeric comparisons
**Warning signs:** Tests pass/fail inconsistently around boundary values
```bash
# Wrong
if [[ "$ratio" == "1.0" ]]; then ...

# Right
if python3 -c "import sys; sys.exit(0 if abs($ratio - 1.0) < 1e-4 else 1)"; then ...
# Or
if (( $(echo "$ratio - 1.0" | bc -l | tr -d '-') < 0.0001 )); then ...
```

### Pitfall 2: Unquoted Variables with Spaces
**What goes wrong:** `$filename` with spaces breaks commands
**Why it happens:** Word splitting on unquoted variables
**How to avoid:** Always quote: `"$filename"`
**Warning signs:** Tests fail on paths with spaces, or inconsistently

### Pitfall 3: Missing set -e in Subshells
**What goes wrong:** Errors in `$()` subshells don't propagate
**Why it happens:** Subshells don't inherit errexit
**How to avoid:** Check return codes explicitly or use `|| exit 1`
```bash
# Risky
result=$(failing_command)

# Safe
result=$(failing_command) || { echo "Failed"; exit 1; }
```

### Pitfall 4: Race Conditions in Parallel Tests
**What goes wrong:** Tests interfere with each other's files
**Why it happens:** Tests write to same directory without isolation
**How to avoid:** Each test uses its own temp directory
```bash
WORK_DIR=$(mktemp -d)
trap "rm -rf $WORK_DIR" EXIT
```

### Pitfall 5: Wrong Energy Units
**What goes wrong:** Energy check fails even with correct physics
**Why it happens:** Unit conversion mismatch between McLuster and ABYSS
**How to avoid:** Verify unit chain: McLuster (Msun, pc, km/s) -> ABYSS input (1e-9 Msun, kpc, km/s) -> code units
**Warning signs:** Energy values off by factors of 1e3, 1e6, 1e9

## Code Examples

Verified patterns from existing ABYSS codebase:

### Test Helper Functions (based on workflow/bin/common.sh)
```bash
#!/usr/bin/env bash
# tests/mcluster/lib/common.sh
set -euo pipefail

# Get repository root
test_repo_root() {
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    cd "$script_dir/../../.." && pwd
}

# Print test status
test_status() {
    local name="$1"
    local status="$2"  # PASS or FAIL
    local detail="${3:-}"

    printf "%-40s %s\n" "$name" "$status"
    if [[ -n "$detail" && "$status" == "FAIL" ]]; then
        echo "  $detail"
    fi
}

# Create isolated work directory
test_setup_workdir() {
    local test_name="$1"
    local workdir
    workdir=$(mktemp -d -t "abyss_test_${test_name}_XXXXXX")
    echo "$workdir"
}

# Cleanup work directory
test_cleanup() {
    local workdir="$1"
    [[ -d "$workdir" ]] && rm -rf "$workdir"
}

# Check if McLuster binary is available
test_has_mcluster() {
    local repo_root
    repo_root="$(test_repo_root)"
    [[ -x "$repo_root/src/mcluster" ]] || [[ -x "$repo_root/mcluster/mcluster_sse" ]]
}

# Get McLuster binary path
test_mcluster_path() {
    local repo_root
    repo_root="$(test_repo_root)"
    if [[ -x "$repo_root/src/mcluster" ]]; then
        echo "$repo_root/src/mcluster"
    elif [[ -x "$repo_root/mcluster/mcluster_sse" ]]; then
        echo "$repo_root/mcluster/mcluster_sse"
    else
        echo ""
    fi
}
```

### TOML Config Fixture Template
```bash
# Generate test config with McLuster section
generate_test_config() {
    local output_file="$1"
    local n_particles="${2:-1000}"
    local profile="${3:-0}"  # 0=Plummer, 1=King

    cat > "$output_file" << EOF
# Test configuration for McLuster integration
Filename = "nbody.dat"
StopTime = 1.0e5
OutputDirectory = "output"

[numerics]
eta = 0.01
FixNumNeighbor = 100
InitialRadius = 0.2

[output]
dtOutput = 1.0e4

[mcluster]
N = $n_particles
P = $profile
R = 0.8
f = 1
Z = 0.02
generate_only = false
EOF
}
```

### Makefile Test Targets
```makefile
# tests/mcluster/Makefile
.PHONY: test test-quick test-full test-generate test-run test-energy clean

# Default: run quick tests (fixture-based)
test: test-quick

# Quick tests with pre-generated fixtures (no gfortran needed)
test-quick:
	FIXTURE_MODE=1 ./run_tests.sh

# Full tests with live McLuster execution (requires gfortran)
test-full:
	@if ! which gfortran >/dev/null 2>&1; then \
		echo "ERROR: gfortran required for full tests"; exit 1; \
	fi
	FIXTURE_MODE=0 ./run_tests.sh

# Individual test targets
test-generate:
	./test_generate_only.sh

test-run:
	./test_run_only.sh

test-energy:
	./test_energy.sh

clean:
	rm -rf output/
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Inline test scripts | Separate test files + runner | Best practice | Maintainability, isolation |
| Manual test execution | Make targets | Standard | CI/CD integration |
| Full simulation tests | IC-only verification | Phase 29 decision | Fast tests (seconds vs hours) |

**Deprecated/outdated:**
- Bats 0.x: Use bats-core if adopting Bats (not required for this phase)
- `test` command: Prefer `[[ ]]` for bash conditionals

## Open Questions

Things that couldn't be fully resolved:

1. **Exact energy tolerance for virial check**
   - What we know: McLuster generates virial equilibrium (Q=0.5), theory says 2K/|U|=1
   - What's unclear: How much numerical precision loss in unit conversions?
   - Recommendation: Start with 1e-4, adjust based on empirical results

2. **MPI test execution**
   - What we know: ABYSS requires MPI (mpirun -np N)
   - What's unclear: Can tests run with mpirun -np 1 for simplicity? Or need multi-rank?
   - Recommendation: Test with np=1 first; add multi-rank if issues found

3. **HDF5 output validation**
   - What we know: Full pipeline produces HDF5 output
   - What's unclear: Exact format to verify beyond file existence
   - Recommendation: Check file exists, non-empty, has expected groups

## Sources

### Primary (HIGH confidence)
- `/gpfs/home/vjl4366/pkg/ABYSS/workflow/bin/common.sh` - Existing bash patterns
- `/gpfs/home/vjl4366/pkg/ABYSS/tools/analyze_energy.py` - Energy computation patterns
- `/gpfs/home/vjl4366/pkg/ABYSS/src/def.h` - Unit definitions (MASS_UNIT, POSITION_UNIT, etc.)
- `/gpfs/home/vjl4366/pkg/ABYSS/src/mcluster_runner.cpp` - Unit conversion implementation
- `/gpfs/home/vjl4366/pkg/ABYSS/mcluster/README` - McLuster documentation

### Secondary (MEDIUM confidence)
- [Bash Shell Script Best Practices](https://sharats.me/posts/shell-script-best-practices/) - Shell script patterns
- [Bats-core GitHub](https://github.com/bats-core/bats-core) - Test framework reference
- [N-body simulations (Scholarpedia)](http://www.scholarpedia.org/article/N-body_simulations_(gravitational)) - Virial equilibrium theory

### Tertiary (LOW confidence)
- WebSearch results on bash testing patterns - General guidance, verified against codebase

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Using existing ABYSS tools and patterns
- Architecture: HIGH - Based on existing workflow/ structure
- Pitfalls: HIGH - Documented from common bash/testing issues
- Energy verification: MEDIUM - Physics clear, exact tolerance needs empirical validation

**Research date:** 2026-01-20
**Valid until:** 2026-02-20 (30 days - stable domain)
