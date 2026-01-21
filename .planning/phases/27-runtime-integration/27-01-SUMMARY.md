---
phase: 27-runtime-integration
plan: 01
subsystem: runtime
tags: [subprocess, posix, fork-exec, mcluster, ic-generation]

# Dependency graph
requires:
  - phase: 26-config-parser-extension
    provides: MclusterConfig struct with parsed TOML config
provides:
  - McLuster subprocess execution module (mcluster_runner.h/cpp)
  - buildMclusterArgs() for config-to-CLI conversion
  - runMclusterSubprocess() with fork/exec and output capture
  - validateMclusterOutput() for file/format validation
  - transformMclusterOutput() for column reordering
affects: [27-02 main.cpp integration, 28-output-format-handling]

# Tech tracking
tech-stack:
  added: []
  patterns: [fork-exec subprocess with pipe capture]

key-files:
  created:
    - src/mcluster_runner.h
    - src/mcluster_runner.cpp
  modified:
    - src/Makefile

key-decisions:
  - "Fork/exec pattern over system()/popen() for proper exit code and stderr capture"
  - "Output transformation in dedicated function rather than modifying readData()"
  - "Scientific notation with 8 decimal places for transformed output"

patterns-established:
  - "RunResult struct pattern for subprocess execution results"
  - "POSIX fork/exec with pipe-based output capture"

# Metrics
duration: 2min
completed: 2026-01-21
---

# Phase 27 Plan 01: McLuster Runner Module Summary

**POSIX fork/exec subprocess execution module with output validation and McLuster-to-ABYSS format transformation**

## Performance

- **Duration:** 2 min
- **Started:** 2026-01-21T00:37:33Z
- **Completed:** 2026-01-21T00:39:38Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments
- Created mcluster_runner.h with RunResult struct and 4 function declarations
- Implemented runMclusterSubprocess() with fork/exec and pipe-based stdout/stderr capture
- Implemented validateMclusterOutput() with file existence, size, header, and line count checks
- Implemented transformMclusterOutput() for column reordering (McLuster: mass x y z vx vy vz -> ABYSS: x y z vx vy vz mass)
- Added mcluster_runner.cpp to build system (311 lines of implementation)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create mcluster_runner.h header** - `4756cac` (feat)
2. **Task 2: Implement mcluster_runner.cpp** - `1810f33` (feat)
3. **Task 3: Add mcluster_runner to build system** - `4e4731f` (chore)

## Files Created/Modified
- `src/mcluster_runner.h` - Header with RunResult struct and function declarations
- `src/mcluster_runner.cpp` - Implementation (311 lines) with fork/exec subprocess execution
- `src/Makefile` - Added mcluster_runner.cpp to CXX_SRCS

## Decisions Made
- Used fork/exec pattern instead of system() or popen() for proper exit code handling and stderr capture
- Implemented output transformation as dedicated function rather than modifying existing readData()
- Used scientific notation with 8 decimal places for transformed output to maintain precision
- M takes precedence over N in argument building (consistent with Phase 26 decision)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- McLuster subprocess execution module complete and ready for integration
- Plan 27-02 will integrate into main.cpp (call from ROOT process after config parsing)
- All 4 functions exported and ready for use
- Build system updated to compile new module

---
*Phase: 27-runtime-integration*
*Completed: 2026-01-21*
