# Project State

## Current Focus

**Milestone:** 1 — AoS to SoA Conversion (v1.0)
**Phase:** 6 — I/O Updates
**Status:** Complete

## Quick Context

Converting ABYSS from Array of Structures (AoS) to Structure of Arrays (SoA) for performance. Full SoA approach with accessor functions. MPI shared memory now uses multiple MPI_Win objects (66 windows for all SoA arrays). CPU routines updated with SoA helpers using bridge pattern. SDAR compatibility maintained via exit-only sync pattern (AoS -> SoA after FewBody operations). I/O updated with SoA sync after initialization and checkpoint restore.

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
| 5. SDAR Compatibility | Complete ✓ | 5 plans, exit-only sync |
| 6. I/O Updates | Complete ✓ | 3 plans, 3 commits |
| 7. Validation | Not Started | Baseline captured |

## Phase 6 Summary

Added SoA sync to initialization and checkpoint restore:
- `root_routines.cpp` — sync AoS to SoA after InitializationRoutines()
- `read_write.cpp` — sync AoS to SoA at end of readCheckpoint()
- HDF5 output verified: writeParticle() reads from AoS (no changes needed)

**Key insight:** HDF5 output works without modification because it reads from AoS, which is kept up-to-date via the exit-only sync pattern from Phase 5.

## Key Decisions

- Full SoA (not hybrid)
- Accessor functions for encapsulation
- SoA for MPI shared memory (multiple MPI_Win)
- Keep SDAR Group as-is (use exit-only sync pattern to bridge)
- GPU: AoS→SoA bridge for backward compat
- I/O: Sync at initialization entry points, read from AoS for output

## Next Action

Run `/gsd:plan-phase 7` to create detailed plan for Validation phase.

---
*Last updated: 2026-01-17 (Phase 6 complete)*
