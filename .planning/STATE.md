# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-18)

**Core value:** Physics correctness (energy conservation) with clean architecture for targeted optimizations
**Current focus:** Load Balance Profiling

## Current Position

**Milestone:** v2.2 — Load Balance Profiling
**Phase:** 17 (Worker Distribution)
**Status:** Ready to plan

Last activity: 2026-01-18 — Phase 16 complete

Progress: ████░░░░░░ 33% (2/6 phases)

## v2.2 Overview

**Load Balance Profiling** (In Progress)

6 phases, 24 requirements:
- Phase 15: Neighbor Profiling (4 requirements) ✓
- Phase 16: Queue Dispatch Profiling (4 requirements) ✓
- Phase 17: Worker Distribution (4 requirements)
- Phase 18: Particle Type Breakdown (5 requirements)
- Phase 19: Memory Access Profiling (3 requirements)
- Phase 20: Analysis & Reporting (4 requirements)

**Goal:** Understand where load imbalance comes from before committing to optimization approach

**Target measurements:**
- Neighbor count variance ✓
- Queue dispatch overhead ✓
- Worker distribution
- CM particle breakdown
- Few-body overhead
- Memory access patterns

## Phase Status

| Phase | Name | Status | Requirements |
|-------|------|--------|--------------|
| 15 | Neighbor Profiling | ✓ Complete | 4/4 |
| 16 | Queue Dispatch Profiling | ✓ Complete | 4/4 |
| 17 | Worker Distribution | ○ Pending | 0/4 |
| 18 | Particle Type Breakdown | ○ Pending | 0/5 |
| 19 | Memory Access Profiling | ○ Pending | 0/3 |
| 20 | Analysis & Reporting | ○ Pending | 0/4 |

## Phase 16 Summary

Completed 2026-01-18 with 4 plans:
- 16-01: Queue depth tracking infrastructure (QueueDepthTracker class)
- 16-02: Instrument queue scheduler (depth sampling, starvation detection)
- 16-03: Dispatch latency tracking (OnlineStats for latency distribution)
- 16-04: Queue dispatch statistics output (CSV/JSON/console)

Key deliverables:
- `QueueDepthTracker` class for pending task monitoring
- `PROFILE_QUEUE_DEPTH`, `PROFILE_STARVATION_EVENT`, `PROFILE_DISPATCH_LATENCY` macros
- Root-side timing breakdown (assign vs wait ratio)
- Bottleneck detection (`isDispatchBottleneck()`)

## Phase 15 Summary

Completed 2026-01-18 with 4 plans:
- 15-01: Online statistics infrastructure (OnlineStats, CorrelationTracker)
- 15-02: Per-particle neighbor counting instrumentation
- 15-03: Neighbor statistics in CSV/JSON/console output
- 15-04: Neighbor count histogram (15 buckets, 0-50K+)

Key deliverables:
- `OnlineStats` class with Welford's algorithm
- `CorrelationTracker` for neighbor-time correlation
- `NeighborHistogram` for distribution analysis
- `PROFILE_NEIGHBOR_TIME` macro for instrumentation

## Milestone History

| Milestone | Name | Phases | Status |
|-----------|------|--------|--------|
| v1.0 | AoS to SoA Conversion | 1-7 (28 plans) | ✓ Shipped 2026-01-17 |
| v2.0 | Performance Profiling & Optimization | 8-10 (11 plans) | ✓ Shipped 2026-01-18 |
| v2.1 | MPI Communication Optimization | 11-14 | ✗ Archived 2026-01-18 |
| v2.2 | Load Balance Profiling | 15-20 (24 reqs) | ◆ In Progress |

## Archive

- v1.0 archived: `.planning/milestones/v1.0-mvp.md`
- v2.0 archived: `.planning/milestones/v2.0-ROADMAP.md`, `.planning/milestones/v2.0-REQUIREMENTS.md`
- v2.1 archived: `.planning/milestones/v2.1-ARCHIVE.md`

## Next Action

Run `/gsd:plan-phase 17` to plan Worker Distribution.

Note: Phases 17-18 can be planned/executed in parallel (17 has no dependencies, 18 depends on 15 which is complete).

---
*Last updated: 2026-01-18 (Phase 16 complete)*
