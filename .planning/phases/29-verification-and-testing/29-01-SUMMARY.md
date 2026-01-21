---
phase: 29-verification-and-testing
plan: 01
subsystem: testing
tags: [bash, python, makefile, virial-ratio, energy-verification]

# Dependency graph
requires:
  - phase: 28-output-format-handling
    provides: "Unit conversion in transformMclusterOutput() (pc->kpc, Msun->1e-9 Msun)"
provides:
  - "Test helper functions (repo_root, status, workdir, cleanup)"
  - "Energy verification computing virial ratio from IC files"
  - "Make targets for test invocation (test, test-quick, test-full)"
affects:
  - "29-02 (energy verification test script)"
  - "29-03 (integration tests)"

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Test helper library (tests/mcluster/lib/common.sh)"
    - "Virial ratio verification (2K/|U|)"
    - "Fixture vs live test modes (FIXTURE_MODE env var)"

key-files:
  created:
    - "tests/mcluster/lib/common.sh"
    - "tests/mcluster/lib/verify_energy.py"
    - "tests/mcluster/Makefile"
  modified: []

key-decisions:
  - "test_* prefix for all helper functions (namespace collision avoidance)"
  - "FIXTURE_MODE env var distinguishes quick (fixture) vs full (live) tests"
  - "Default tolerance 1e-4 for virial ratio verification (adjustable via CLI)"

patterns-established:
  - "Test helpers in tests/mcluster/lib/ directory"
  - "Shell test functions return paths via echo, not setting variables"
  - "Python verification scripts use exit codes (0=PASS, 1=FAIL)"

# Metrics
duration: 2min
completed: 2026-01-21
---

# Phase 29 Plan 01: Test Infrastructure Summary

**Test helper library (common.sh), virial ratio verification (verify_energy.py), and Make targets for test invocation**

## Performance

- **Duration:** 2 min 31 sec
- **Started:** 2026-01-21T02:10:25Z
- **Completed:** 2026-01-21T02:12:56Z
- **Tasks:** 3
- **Files created:** 3

## Accomplishments

- Test helper library with 8 functions for consistent test execution
- Energy verification script computing virial ratio Q = 2K/|U| from IC files
- Makefile with test, test-quick, test-full targets for different modes

## Task Commits

Each task was committed atomically:

1. **Task 1: Create test helper library (common.sh)** - `0f9f113` (feat)
2. **Task 2: Create energy verification script (verify_energy.py)** - `d0b8bb1` (feat)
3. **Task 3: Create test Makefile** - `f82e88b` (feat)

## Files Created

- `tests/mcluster/lib/common.sh` - Test helper functions (repo_root, status, workdir, cleanup, mcluster_path, abyss_path, die)
- `tests/mcluster/lib/verify_energy.py` - Virial ratio computation for IC file verification (uses src/def.h unit constants)
- `tests/mcluster/Makefile` - Test invocation targets (test, test-quick, test-full, clean)

## Decisions Made

- **test_* function prefix:** All helper functions use test_ prefix to avoid namespace collisions with workflow scripts
- **FIXTURE_MODE environment variable:** Separates quick tests (pre-generated fixtures) from full tests (live McLuster execution)
- **Default tolerance 1e-4:** Virial ratio verification uses 1e-4 tolerance by default, adjustable via CLI argument
- **Exit codes for verification:** Python verification scripts use exit 0 for PASS, exit 1 for FAIL (shell script friendly)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Test infrastructure ready for Plan 29-02 (energy verification test script)
- common.sh can be sourced by test scripts: `source lib/common.sh`
- verify_energy.py can be invoked: `python lib/verify_energy.py <ic_file> [tolerance]`
- Make targets ready: `make test-quick` (fixture mode) or `make test-full` (live mode)

---
*Phase: 29-verification-and-testing*
*Completed: 2026-01-21*
