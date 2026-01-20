# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-20)

**Core value:** Physics correctness (energy conservation) with clean architecture for targeted optimizations
**Current focus:** v3.0 McCluster Integration - Config Parser Extension Complete

## Current Position

**Milestone:** v3.0 McCluster Integration
**Phase:** 26 - Config Parser Extension (complete)
**Plan:** 02 complete, phase complete
**Status:** Phase 26 complete, ready for Phase 27

Last activity: 2026-01-20 - Completed 26-02-PLAN.md (Config Parser Function)

Progress: [###-------] 32% (1.6/5 phases complete)

## v3.0 Milestone Overview

**Goal:** Integrate McLuster IC generator for seamless end-to-end simulation workflow

**Phases:**
| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 25 | Build System Integration | 4 | Complete |
| 26 | Config Parser Extension | 10 | Complete |
| 27 | Runtime Integration | 7 | Pending |
| 28 | Output Format Handling | 3 | Pending |
| 29 | Verification and Testing | 4 | Pending |

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

### Technical Notes

- McLuster already downloaded in `mcluster/` directory
- McLuster Makefile builds `mcluster_sse` with gfortran + gcc
- Root-level Makefile orchestrates both builds (Phase 25 complete)
- src/mcluster symlink created at build time
- ABYSS uses TOML config (src/read_parameter_file.cpp, src/toml.hpp)
- ABYSS expects nbody.dat format: x y z vx vy vz mass (7 columns)
- McLuster `-C 3 -u 0` outputs compatible ASCII format
- toml::value now has keys() method for table key iteration (Phase 26-01)
- MclusterConfig struct defined with all 9 parameters (Phase 26-01)
- parseMclusterSection() parses [mcluster] TOML section (Phase 26-02)
- validateMclusterConfig() enforces N/M mutual exclusivity (Phase 26-02)
- mcluster_config global accessible via extern in global.h (Phase 26-02)

### TODOs

- [x] Complete Phase 25 Build System Integration
- [x] Complete Phase 26 Plan 01 (Config Infrastructure)
- [x] Complete Phase 26 Plan 02 (Config Parser Function)
- [ ] Complete Phase 27 Runtime Integration

### Blockers

None

## Session Continuity

**For next session:**
1. Phase 26 complete - all mcluster config parsing in place
2. parseMclusterSection() ready to be called
3. mcluster_config global accessible from main.cpp
4. Ready for Phase 27: Runtime Integration (mcluster binary invocation)

**Key files for v3.0:**
- `Makefile` - Root-level build orchestration (Phase 25)
- `mcluster/Makefile` - existing McLuster build (gfortran + gcc)
- `src/main.cpp` - entry point (add mcluster invocation in Phase 27)
- `src/read_parameter_file.cpp` - TOML parser with mcluster parsing (Phase 26-02)
- `src/toml.hpp` - TOML library (extended with keys() method - Phase 26-01)
- `src/mcluster_config.h` - MclusterConfig struct definition (Phase 26-01)
- `src/global.h` - extern MclusterConfig declaration (Phase 26-02)

**Phase 26 outputs:**
- toml::value keys() method for table key iteration
- MclusterConfig struct with 9 parameters and defaults
- MCLUSTER_VALID_PARAMS array for unknown key detection
- parseMclusterSection() for TOML parsing
- validateMclusterConfig() for N/M mutual exclusivity
- levenshteinDistance() and suggestSimilarParam() for typo detection
- mcluster_config extern declaration in global.h
- SUMMARY: `.planning/phases/26-config-parser-extension/26-02-SUMMARY.md`

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
*Last updated: 2026-01-20 (Phase 26 complete)*
