# Project State

## Current Focus

**Milestone:** 1 — AoS to SoA Conversion (v1.0)
**Phase:** 4 — CPU Routines
**Status:** Complete

## Quick Context

Converting ABYSS from Array of Structures (AoS) to Structure of Arrays (SoA) for performance. Full SoA approach with accessor functions. MPI shared memory now uses multiple MPI_Win objects (66 windows for all SoA arrays). CPU routines updated with SoA helpers using bridge pattern.

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
| 4. CPU Routines | Complete ✓ | 5 plans, 5 commits |
| 5. SDAR Compatibility | Not Started | — |
| 6. I/O Updates | Not Started | — |
| 7. Validation | Not Started | Baseline captured |

## Phase 4 Summary

Added SoA-compatible CPU routines using bridge pattern:
- `particle_data.h/cpp` — sync_from_particle(), sync_to_particle(), bulk accessors
- `timestep_routines.cpp` — SoA overloads for getNewTimeStepReg/Irr
- `update_particle.cpp` — predict_second_order(), correct_fourth_order() free functions
- `compute_acceleration.cpp` — AccumulatorSoA struct, predict_neighbor_soa()
- `regular_routines.cpp`, `irregular_routines.cpp` — SoA integration

Commits: `938a7bc`, `91066de`, `5d8bcbe`, `b926bda`, `03ba526`

## Key Decisions

- Full SoA (not hybrid)
- Accessor functions for encapsulation
- SoA for MPI shared memory (multiple MPI_Win)
- Keep SDAR Group as-is (use proxy pattern)
- GPU: AoS→SoA bridge for backward compat (deferred to Phase 4)

## Next Action

Run `/gsd:plan-phase 5` to create detailed plan for SDAR Compatibility phase.

---
*Last updated: 2026-01-17 (Phase 4 complete)*
