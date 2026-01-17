# Project State

## Current Focus

**Milestone:** 1 — AoS to SoA Conversion (v1.0)
**Phase:** 2 — MPI Integration
**Status:** Complete

## Quick Context

Converting ABYSS from Array of Structures (AoS) to Structure of Arrays (SoA) for performance. Full SoA approach with accessor functions. MPI shared memory now uses multiple MPI_Win objects (66 windows for all SoA arrays).

## Baseline Metrics

| Metric | Value | Source |
|--------|-------|--------|
| dE/E0 mean | 3.20665e-05 | baseline_20260116_234947 |
| Wall time | 28.33s | baseline_20260116_234947 |
| Git commit | 7e90c74 | stable_candidate branch |

## Phase Progress

| Phase | Status | Notes |
|-------|--------|-------|
| 1. Core SoA Container | Complete ✓ | 3 plans, 3 commits |
| 2. MPI Integration | Complete ✓ | 4 plans, 8 commits |
| 3. GPU Integration | Complete ✓ | 5 plans, 5 commits |
| 4. CPU Routines | Not Started | — |
| 5. SDAR Compatibility | Not Started | — |
| 6. I/O Updates | Not Started | — |
| 7. Validation | Not Started | Baseline captured |

## Phase 2 Summary

Created MPI shared memory infrastructure for SoA:
- `particle_data_mpi.h/cpp` — ParticleDataMPI class with 66 MPI_Win handles
- Updated `global.h`, `default_global.cpp` with extern particle_data
- Updated `mpi_routines.cpp` with allocate_shared(), sync_all(), timing
- Updated `main.cpp` with clean shutdown (deallocate_shared)
- DEBUG_MPI verification for cross-rank access

Commits: `b57b92e`, `3d06888`, `1b666dd`, `cc34c57`, `dae2e28`, `1c82ebc`, `9000806`, `a48af2c`, `665a8ad`

## Key Decisions

- Full SoA (not hybrid)
- Accessor functions for encapsulation
- SoA for MPI shared memory (multiple MPI_Win)
- Keep SDAR Group as-is (use proxy pattern)
- GPU: AoS→SoA bridge for backward compat (deferred to Phase 4)

## Next Action

Run `/gsd:plan-phase 4` to create detailed plan for CPU Routines phase.

---
*Last updated: 2026-01-17 (Phase 2 complete)*
