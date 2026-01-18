# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-17)

**Core value:** Identify and eliminate irregular force bottlenecks through systematic profiling
**Current focus:** Phase 8 — Enhanced Profiler Infrastructure

## Current Position

**Milestone:** v2.0 — Performance Profiling & Optimization
**Phase:** 8 (first of v2.0)
**Status:** Planned (6 plans in 3 waves)

Last activity: 2026-01-17 — v2.0 milestone initialized

Progress: ░░░░░░░░░░ 0%

## Phase Overview

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 8 | Enhanced Profiler Infrastructure | 16 | ◆ Planned |
| 9 | Profile & Analyze | 3 | ○ Pending |
| 10 | Targeted Optimization | 3 | ○ Pending |

## v2.0 Focus

**Goal:** Build detailed profiler to identify irregular force bottlenecks, then optimize.

**Phase 1 targets:**
- Per-rank MPI statistics (min/avg/max)
- Histogram support for call-time distributions
- Work counters (particles, neighbor pairs)
- Irregular force sub-timers (NeighborLoop, CMLoop, Correction)
- Worker-side timing (recv, dispatch, send)
- Queue scheduler timing (assign, run)

**Key files to modify:**
- `src/profiler.h`
- `src/Particle/compute_acceleration.cpp`
- `src/worker_routines.cpp`
- `src/queue_scheduler.h`

## Previous Milestone

**v1.0 AoS to SoA Conversion** (Shipped 2026-01-17)
- 7 phases, 28 plans, 57 commits
- Energy conservation: ✓ Pass (2.52e-5 vs 3.21e-5 baseline)
- Performance: No change (bottleneck elsewhere)

## Next Action

Run `/gsd:plan-phase 8` to create execution plan for Phase 8.

---
*Last updated: 2026-01-17 (v2.0 initialized)*
