# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-17)

**Core value:** Identify and eliminate irregular force bottlenecks through systematic profiling
**Current focus:** v2.0 Complete — Awaiting next milestone

## Current Position

**Milestone:** v2.0 — Performance Profiling & Optimization
**Phase:** 10 (complete)
**Status:** Milestone v2.0 complete

Last activity: 2026-01-18 — Phase 10 complete, v2.0 milestone complete

Progress: ██████████ 100%

## Phase Overview

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 8 | Enhanced Profiler Infrastructure | 16 | ✓ Complete |
| 9 | Profile & Analyze | 3 | ✓ Complete |
| 10 | Targeted Optimization | 3 | ✓ Complete |

## Phase 10 Summary

**Targeted Optimization** (Completed 2026-01-18)

3 plans in 2 waves:
- Wave 1: Create AVX-512 force kernel, integrate vectorized kernel
- Wave 2: Validation and profiling comparison

Key results:
- **AVX-512 vectorization implemented** for irregular force neighbor loop
- IrregularForce: 69.5s → 67.3s (**-3.2%** improvement)
- Wall time: 130s → 126.4s (-2.8% improvement)
- **Energy conservation verified**: dE/E0 = 1.34e-6 (PASS)
- Target was 20%, achieved 3.2% due to gather overhead and memory-bound workload

## v2.0 Milestone Complete

**Goal:** Build detailed profiler to identify irregular force bottlenecks, then optimize.

**Phase 8 delivered:**
- Per-rank MPI statistics (min/avg/max) ✓
- Histogram support for call-time distributions ✓
- Work counters (particles, neighbor pairs) ✓
- Irregular force sub-timers ✓
- Worker-side timing ✓
- Queue scheduler timing ✓

**Phase 9 delivered:**
- Profiling data collected from test1 simulation ✓
- Primary bottleneck identified: IrregularForce (53.5%) ✓
- Optimization targets documented ✓

**Phase 10 delivered:**
- AVX-512 vectorization of irregular force loop ✓
- Performance improvement: 3.2% (below 20% target) ⚠
- Energy conservation maintained ✓
- Recommendations for v2.1 documented ✓

## Milestone Statistics

- **Phases:** 3
- **Plans executed:** 11 (6 + 2 + 3)
- **Requirements completed:** 22/22
- **Performance gain:** 3.2% (limited by memory-bound workload)
- **Energy conservation:** ✓ Maintained

## Recommendations for v2.1

1. **MPI batching** — Reduce ~105M messages/interval overhead
2. **SIMD gather intrinsics** — Use `_mm512_i64gather_pd` to avoid pre-gather copies
3. **GPU irregular forces** — Port irregular force kernel to GPU

## Previous Milestone

**v1.0 AoS to SoA Conversion** (Shipped 2026-01-17)
- 7 phases, 28 plans, 57 commits
- Energy conservation: ✓ Pass (2.52e-5 vs 3.21e-5 baseline)
- Performance: No change (bottleneck elsewhere)

## Next Action

Run `/gsd:new-milestone` to start v2.1 or conclude the project.

---
*Last updated: 2026-01-18 (v2.0 milestone complete)*
