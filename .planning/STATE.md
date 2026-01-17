# Project State

## Current Focus

**Milestone:** 1 — AoS to SoA Conversion (v1.0)
**Phase:** 2 — MPI Integration
**Status:** Not Started

## Quick Context

Converting ABYSS from Array of Structures (AoS) to Structure of Arrays (SoA) for performance. Full SoA approach with accessor functions. MPI shared memory will use multiple `MPI_Win` objects.

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
| 2. MPI Integration | Not Started | — |
| 3. GPU Integration | Not Started | — |
| 4. CPU Routines | Not Started | — |
| 5. SDAR Compatibility | Not Started | — |
| 6. I/O Updates | Not Started | — |
| 7. Validation | Not Started | Baseline captured |

## Key Decisions

- Full SoA (not hybrid)
- Accessor functions for encapsulation
- SoA for MPI shared memory (multiple MPI_Win)
- Keep SDAR Group as-is (use proxy pattern)

## Next Action

Run `/gsd:plan-phase 2` to create detailed plan for MPI Integration phase.

---
*Last updated: 2026-01-17 (Phase 1 complete)*
