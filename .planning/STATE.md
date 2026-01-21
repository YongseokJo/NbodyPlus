# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-20)

**Core value:** Physics correctness (energy conservation) with clean architecture for targeted optimizations
**Current focus:** v3.0 McCluster Integration - Verification and Testing

## Current Position

**Milestone:** v3.0 McCluster Integration
**Phase:** 29 - Verification and Testing (in progress)
**Plan:** 02 complete
**Status:** Plan 29-02 complete, ready for 29-03

Last activity: 2026-01-21 - Completed Plan 29-02 (Test Scripts Implementation)

Progress: [########=-] 90% (4/5 phases, Plan 29-02 complete)

## v3.0 Milestone Overview

**Goal:** Integrate McLuster IC generator for seamless end-to-end simulation workflow

**Phases:**
| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 25 | Build System Integration | 4 | Complete |
| 26 | Config Parser Extension | 10 | Complete |
| 27 | Runtime Integration | 7 | Complete |
| 28 | Output Format Handling | 3 | Complete |
| 29 | Verification and Testing | 4 | In Progress (2/4) |

**Total:** 28 requirements across 5 phases

## Milestone History

| Milestone | Name | Phases | Status |
|-----------|------|--------|--------|
| v1.0 | AoS to SoA Conversion | 1-7 (28 plans) | Shipped 2026-01-17 |
| v2.0 | Performance Profiling & Optimization | 8-10 (11 plans) | Shipped 2026-01-18 |
| v2.1 | MPI Communication Optimization | 11-14 | Archived 2026-01-18 |
| v2.2 | Load Balance Profiling | 15-20 (24 reqs) | Complete 2026-01-19 |
| v2.3 | Performance Analysis & Instrumentation | 21-23 (11 plans) | Shipped 2026-01-20 |
| v3.0 | McCluster Integration | 25-29 (28 reqs) | Active |

## Accumulated Context

### Key Decisions (v3.0)

| Decision | Rationale |
|----------|-----------|
| 5 phases for 28 reqs | Natural boundaries align with requirement categories |
| Linear phase dependencies | Each phase unlocks next (BUILD->CONFIG->RUNTIME->OUTPUT->VERIFY) |
| Start at Phase 25 | Continues numbering from v2.3 (ended at Phase 24 deferred) |
| Symlink mcluster binary | src/mcluster -> ../mcluster/mcluster_sse (single source of truth) |
| Warning not error for missing gfortran | ABYSS can still build without McLuster |
| DISABLE_MCLUSTER=1 flag | Explicit skip even when gfortran available |
| keys() returns empty for non-tables | Defensive design - no exception thrown |
| MclusterConfig defaults match McLuster | P=0, R=0.8, f=1, Z=0.02 verified from source |
| MCLUSTER_VALID_PARAMS as const vector | Enables easy iteration for unknown key detection |
| M precedence over N | When both specified, M takes precedence with warning |
| Levenshtein threshold 2 | Typo suggestions only for edit distance <= 2 |
| N minimum 3 | N-body simulation requires at least 3 particles |
| Fork/exec over system()/popen() | Proper exit code handling and stderr capture (Phase 27-01) |
| Output transformation vs readData() mod | Dedicated function keeps readData() unchanged (Phase 27-01) |
| ROOT-only McLuster execution | Avoid parallel subprocess conflicts (Phase 27-02) |
| MPI_Bcast for rank synchronization | Propagate success/failure to all ranks (Phase 27-02) |
| Clean MPI cleanup on generate_only | Proper resource deallocation before early exit (Phase 27-02) |
| In-place unit conversion | x/1000.0 in output statement cleaner than separate variables (Phase 28-01) |
| 15-digit precision | Full double precision for small kpc position values (Phase 28-01) |
| test_* function prefix | Namespace collision avoidance with workflow scripts (Phase 29-01) |
| FIXTURE_MODE env var | Distinguishes quick (fixture) vs full (live McLuster) tests (Phase 29-01) |
| Default tolerance 1e-4 | Virial ratio verification tolerance, adjustable via CLI (Phase 29-01) |
| Graceful numpy degradation | test_energy.sh passes in fixture mode without numpy (Phase 29-02) |
| $((VAR + 1)) over ((VAR++)) | Avoid set -e exit when counter starts at 0 (Phase 29-02) |

### Technical Notes

- McLuster already downloaded in `mcluster/` directory
- McLuster Makefile builds `mcluster_sse` with gfortran + gcc
- Root-level Makefile orchestrates both builds (Phase 25 complete)
- src/mcluster symlink created at build time
- ABYSS uses TOML config (src/read_parameter_file.cpp, src/toml.hpp)
- ABYSS expects nbody.dat format: x y z vx vy vz mass (7 columns)
- McLuster `-C 3 -u 1` outputs astrophysical units (Msun, pc, km/s)
- toml::value now has keys() method for table key iteration (Phase 26-01)
- MclusterConfig struct defined with all 9 parameters (Phase 26-01)
- parseMclusterSection() parses [mcluster] TOML section (Phase 26-02)
- validateMclusterConfig() enforces N/M mutual exclusivity (Phase 26-02)
- mcluster_config global accessible via extern in global.h (Phase 26-02)
- mcluster_runner.h/cpp provides subprocess execution (Phase 27-01)
- main.cpp now has McLuster orchestration between readParameterFile() and readData() (Phase 27-02)
- transformMclusterOutput() converts pc->kpc (divide by 1000) and Msun->1e-9 Msun (divide by 1e9) (Phase 28-01)
- tests/mcluster/lib/common.sh provides test helper functions (Phase 29-01)
- tests/mcluster/lib/verify_energy.py computes virial ratio for IC verification (Phase 29-01)
- tests/mcluster/Makefile has test, test-quick, test-full targets (Phase 29-01)
- tests/mcluster/fixtures/*.toml provides TOML test configs (Phase 29-02)
- tests/mcluster/test_*.sh scripts implement VERIFY-01 through VERIFY-04 (Phase 29-02)
- tests/mcluster/run_tests.sh main test runner with FIXTURE_MODE support (Phase 29-02)

### TODOs

- [x] Complete Phase 25 Build System Integration
- [x] Complete Phase 26 Plan 01 (Config Infrastructure)
- [x] Complete Phase 26 Plan 02 (Config Parser Function)
- [x] Complete Phase 27 Plan 01 (McLuster Runner Module)
- [x] Complete Phase 27 Plan 02 (Main Integration)
- [x] Complete Phase 28 Plan 01 (Output Format Handling)
- [x] Complete Phase 29 Plan 01 (Test Infrastructure)
- [x] Complete Phase 29 Plan 02 (Test Scripts Implementation)
- [ ] Complete Phase 29 Plans 03-04 (Additional Integration Tests)

### Blockers

**Pre-existing build issue:** The project has a compilation error in `read_parameter_file.cpp` when using gcc 10.2.0 due to incomplete type issues with the TOML library. This should be addressed separately but does not affect runtime integration code correctness.

## Session Continuity

**For next session:**
1. Phase 29 Plan 02 complete - Test scripts implemented
2. All four VERIFY requirements have test scripts
3. Tests pass in FIXTURE_MODE=1 via `make test`
4. run_tests.sh provides main test runner with status reporting
5. Ready to proceed to Plan 29-03 (Additional Integration Tests if planned)

**Key files for v3.0:**
- `Makefile` - Root-level build orchestration (Phase 25)
- `mcluster/Makefile` - existing McLuster build (gfortran + gcc)
- `src/main.cpp` - McLuster orchestration integrated (Phase 27-02)
- `src/read_parameter_file.cpp` - TOML parser with mcluster parsing (Phase 26-02)
- `src/toml.hpp` - TOML library (extended with keys() method - Phase 26-01)
- `src/mcluster_config.h` - MclusterConfig struct definition (Phase 26-01)
- `src/global.h` - extern MclusterConfig declaration (Phase 26-02)
- `src/mcluster_runner.h` - McLuster runner declarations with unit docs (Phase 27-01, 28-01)
- `src/mcluster_runner.cpp` - McLuster subprocess + unit conversion (Phase 27-01, 28-01)
- `tests/mcluster/lib/common.sh` - Test helper functions (Phase 29-01)
- `tests/mcluster/lib/verify_energy.py` - Virial ratio verification (Phase 29-01)
- `tests/mcluster/Makefile` - Test invocation targets (Phase 29-01)

**Phase 29-01 outputs:**
- Test helper library with 8 functions (repo_root, status, workdir, cleanup, etc.)
- Energy verification script computing virial ratio
- Makefile with test, test-quick, test-full targets
- SUMMARY: `.planning/phases/29-verification-and-testing/29-01-SUMMARY.md`

**Phase 29-02 outputs:**
- TOML fixture configs for Plummer, generate-only, and run-only scenarios
- Four test scripts implementing VERIFY-01 through VERIFY-04
- Main test runner with FIXTURE_MODE and FAIL_FAST support
- SUMMARY: `.planning/phases/29-verification-and-testing/29-02-SUMMARY.md`

## Deferred Work (v2.4+)

MPI batching optimization deferred from v2.3:
- Phase 24 plans exist in `.planning/phases/24-mpi-batching/`
- Expected 16-27% speedup from reducing 10.3M messages
- Can be resumed after v3.0 completion

## Archive

- v1.0 archived: `.planning/milestones/v1.0-mvp.md`
- v2.0 archived: `.planning/milestones/v2.0-ROADMAP.md`, `.planning/milestones/v2.0-REQUIREMENTS.md`
- v2.1 archived: `.planning/milestones/v2.1-ARCHIVE.md`
- v2.3 archived: `.planning/milestones/v2.3-ROADMAP.md`, `.planning/milestones/v2.3-REQUIREMENTS.md`
- Phase 23 analysis: `.planning/ANALYSIS.md`

---
*Last updated: 2026-01-21 (Plan 29-02 complete)*
