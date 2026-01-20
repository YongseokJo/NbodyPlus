# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-20)

**Core value:** Physics correctness (energy conservation) with clean architecture for targeted optimizations
**Current focus:** v3.0 McCluster Integration - Config Parser Extension

## Current Position

**Milestone:** v3.0 McCluster Integration
**Phase:** 25 - Build System Integration (complete)
**Plan:** 01 complete
**Status:** Phase 25 complete, ready for Phase 26

Last activity: 2026-01-20 - Completed 25-01-PLAN.md (Build System Integration)

Progress: [##--------] 20% (1/5 phases complete)

## v3.0 Milestone Overview

**Goal:** Integrate McLuster IC generator for seamless end-to-end simulation workflow

**Phases:**
| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 25 | Build System Integration | 4 | Complete |
| 26 | Config Parser Extension | 10 | Pending |
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

### Technical Notes

- McLuster already downloaded in `mcluster/` directory
- McLuster Makefile builds `mcluster_sse` with gfortran + gcc
- Root-level Makefile orchestrates both builds (Phase 25 complete)
- src/mcluster symlink created at build time
- ABYSS uses TOML config (src/read_parameter_file.cpp, src/toml.hpp)
- ABYSS expects nbody.dat format: x y z vx vy vz mass (7 columns)
- McLuster `-C 3 -u 0` outputs compatible ASCII format

### TODOs

- [x] Complete Phase 25 Build System Integration
- [ ] Plan Phase 26 Config Parser Extension

### Blockers

None

## Session Continuity

**For next session:**
1. Phase 25 complete - root Makefile builds ABYSS + McLuster
2. Phase 26 ready for planning (Config Parser Extension)
3. Key file: `src/read_parameter_file.cpp` for config parser extension
4. McLuster binary accessible at src/mcluster after build

**Key files for v3.0:**
- `Makefile` - Root-level build orchestration (NEW - Phase 25)
- `mcluster/Makefile` - existing McLuster build (gfortran + gcc)
- `src/main.cpp` - entry point (add mcluster detection)
- `src/read_parameter_file.cpp` - TOML parser (add [mcluster] section)
- `src/toml.hpp` - TOML library

**Phase 25 outputs:**
- Root Makefile with gfortran detection and conditional mcluster build
- src/mcluster symlink to mcluster_sse binary
- SUMMARY: `.planning/phases/25-build-system-integration/25-01-SUMMARY.md`

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
*Last updated: 2026-01-20 (Phase 25-01 complete)*
