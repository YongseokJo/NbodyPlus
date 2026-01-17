# Project State

## Current Focus

**Milestone:** 1 — AoS to SoA Conversion (v1.0)
**Phase:** 5 — SDAR Compatibility
**Status:** Complete

## Quick Context

Converting ABYSS from Array of Structures (AoS) to Structure of Arrays (SoA) for performance. Full SoA approach with accessor functions. MPI shared memory now uses multiple MPI_Win objects (66 windows for all SoA arrays). CPU routines updated with SoA helpers using bridge pattern. SDAR compatibility maintained via exit-only sync pattern (AoS -> SoA after FewBody operations).

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
| 6. I/O Updates | Not Started | — |
| 7. Validation | Not Started | Baseline captured |

## Phase 5 Summary

Added SoA sync to FewBody operations for SDAR compatibility:
- `particle_data.h/cpp` — batch sync helpers (sync_from_particles, sync_to_particles, etc.)
- `fb_initialization.cpp` — exit sync: CM to SoA, members marked inactive
- `fb_termination.cpp` — exit sync: members to SoA, CM marked inactive
- `fb_integration.cpp` — exit sync: members after AR integration, merger products

**Key design decision:** Exit-only sync pattern. particles[] is authoritative during FewBody operations. Entry syncs (SoA -> AoS) were removed because they corrupted particle data.

## Key Decisions

- Full SoA (not hybrid)
- Accessor functions for encapsulation
- SoA for MPI shared memory (multiple MPI_Win)
- Keep SDAR Group as-is (use exit-only sync pattern to bridge)
- GPU: AoS→SoA bridge for backward compat

## Next Action

Run `/gsd:plan-phase 6` to create detailed plan for I/O Updates phase.

---
*Last updated: 2026-01-17 (Phase 5 complete)*
