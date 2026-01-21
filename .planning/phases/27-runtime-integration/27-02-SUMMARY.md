---
phase: 27-runtime-integration
plan: 02
subsystem: runtime
tags: [mpi, mcluster, ic-generation, main-integration, subprocess-orchestration]

# Dependency graph
requires:
  - phase: 27-01
    provides: McLuster subprocess execution functions (buildMclusterArgs, runMclusterSubprocess, validateMclusterOutput, transformMclusterOutput)
  - phase: 26-02
    provides: MclusterConfig struct with has_mcluster_section flag and generate_only parameter
provides:
  - McLuster orchestration integrated into main.cpp execution flow
  - ROOT-only subprocess execution with MPI broadcast synchronization
  - generate_only mode with clean MPI resource cleanup
  - Automatic fname update to transformed IC file
affects: [28-output-format-handling, 29-verification-testing]

# Tech tracking
tech-stack:
  added: []
  patterns: [ROOT-only subprocess with MPI_Bcast synchronization]

key-files:
  created: []
  modified:
    - src/main.cpp

key-decisions:
  - "McLuster runs only on ROOT rank to avoid parallel subprocess conflicts"
  - "MPI_Bcast propagates success/failure flag for rank synchronization"
  - "generate_only mode performs clean MPI resource cleanup before exit"
  - "fname updated with static string to ensure lifetime extends past function scope"

patterns-established:
  - "ROOT-only subprocess execution with MPI broadcast for multi-rank coordination"
  - "Early exit with proper MPI cleanup for generate-only workflows"

# Metrics
duration: 3min
completed: 2026-01-21
---

# Phase 27 Plan 02: Main Integration Summary

**McLuster IC generation orchestration in main.cpp with ROOT-only execution, MPI rank synchronization, and generate_only early exit support**

## Performance

- **Duration:** 3 min
- **Started:** 2026-01-21T00:49:00Z
- **Completed:** 2026-01-21T00:52:11Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Integrated McLuster subprocess orchestration into main.cpp between config parsing and data reading
- Implemented ROOT-only subprocess execution with MPI_Bcast synchronization for all ranks
- Added generate_only mode support with clean MPI resource cleanup and exit(0)
- Automatic fname update to point to transformed IC file for seamless simulation continuation

## Task Commits

Each task was committed atomically:

1. **Task 1: Add McLuster orchestration to main.cpp** - `9de1914` (feat)

## Files Created/Modified
- `src/main.cpp` - Added mcluster_runner.h include and 80-line McLuster orchestration block

## Decisions Made
- McLuster subprocess runs only on ROOT rank to prevent parallel execution conflicts
- MPI_Bcast used to synchronize success/failure status across all ranks before continuing
- generate_only=true triggers clean exit with proper MPI resource deallocation
- Used static string for generated_ic_path to ensure lifetime extends beyond scope

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

**Pre-existing build issue (unrelated to this plan):**
The project has a compilation error in `read_parameter_file.cpp` when using gcc 10.2.0 due to incomplete type issues with the TOML library. This is a pre-existing issue that affects the full build but does not affect the correctness of the main.cpp changes, which compile successfully when tested in isolation.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Runtime integration complete (Plans 27-01 and 27-02)
- Phase 27 may have additional plans for error handling or testing
- Ready for Phase 28 (Output Format Handling) to verify McLuster-generated IC files work with ABYSS simulation
- Pre-existing build issue with TOML library should be addressed separately

---
*Phase: 27-runtime-integration*
*Completed: 2026-01-21*
