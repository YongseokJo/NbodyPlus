# Roadmap: ABYSS v2.0 Performance Profiling & Optimization

**Created:** 2026-01-17
**Core Value:** Identify and eliminate irregular force bottlenecks through systematic profiling

## Phase Overview

| Phase | Name | Requirements | Status |
|-------|------|--------------|--------|
| 8 | Enhanced Profiler Infrastructure | PROF-01 to PROF-05, IRR-01 to IRR-05, MPI-01 to MPI-06 | Complete |
| 9 | Profile & Analyze | ANLZ-01 to ANLZ-03 | Pending |
| 10 | Targeted Optimization | OPT-01 to OPT-03 | Pending |

---

## Phase 8: Enhanced Profiler Infrastructure

**Goal:** Build comprehensive profiling infrastructure with per-rank statistics, sub-timers, work counters, and load balance metrics.

**Requirements:** PROF-01, PROF-02, PROF-03, PROF-04, PROF-05, IRR-01, IRR-02, IRR-03, IRR-04, IRR-05, MPI-01, MPI-02, MPI-03, MPI-04, MPI-05, MPI-06

**Key Files:**
- `src/profiler.h` — Core profiler enhancements (per-rank, histograms, counters)
- `src/Particle/compute_acceleration.cpp` — Irregular force sub-timers
- `src/worker_routines.cpp` — Worker-side timing
- `src/queue_scheduler.h` — Queue timing

**Success Criteria:**
1. Running with PERFORMANCETRACE produces per-rank min/avg/max statistics in output
2. Irregular force breakdown shows time spent in NeighborLoop vs CMLoop vs Correction
3. Work counter shows neighbor pairs evaluated per timestep
4. Load balance report shows work distribution across MPI ranks
5. All new timers compile and run without affecting simulation correctness

**Estimated Plans:** 4-6

**Completed:** 2026-01-17 (6 plans in 3 waves)

---

## Phase 9: Profile & Analyze

**Goal:** Collect profiling data from representative simulations and identify the primary bottleneck with quantitative evidence.

**Requirements:** ANLZ-01, ANLZ-02, ANLZ-03

**Depends On:** Phase 8 complete

**Key Activities:**
- Run standard test simulation with enhanced profiler
- Collect timing breakdown across all components
- Analyze load imbalance across ranks
- Identify which component dominates wall time

**Success Criteria:**
1. Profiling data collected from at least one representative simulation run
2. Primary bottleneck identified (component taking >X% of wall time)
3. Load imbalance quantified (max/avg ratio documented)
4. Analysis documented with specific optimization targets

**Estimated Plans:** 2-3

---

## Phase 10: Targeted Optimization

**Goal:** Implement optimization for the identified bottleneck and validate with before/after comparison.

**Requirements:** OPT-01, OPT-02, OPT-03

**Depends On:** Phase 9 complete (bottleneck identified)

**Potential Optimizations (depending on Phase 2 findings):**
- **If force loop is bottleneck:** SIMD vectorization of neighbor loop
- **If MPI communication is bottleneck:** Reduce message count, batch sends
- **If load imbalance is bottleneck:** Improve work distribution algorithm
- **If queue scheduling is bottleneck:** Optimize assignment algorithm

**Success Criteria:**
1. Optimization implemented for primary bottleneck
2. Before/after profiling shows measurable improvement
3. Energy conservation verified (within baseline tolerance: dE/E0 < 1e-4)
4. No regression in simulation correctness

**Estimated Plans:** 3-5

---

## Milestone Success Criteria

v2.0 is complete when:
1. All 22 requirements marked complete
2. Bottleneck identified and addressed
3. Performance improvement quantified (or documented if no improvement possible)
4. Energy conservation maintained

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Profiler overhead affects measurements | Use lightweight instrumentation, verify overhead < 5% |
| Bottleneck is in GPU code | Document finding, defer GPU optimization to v2.1 |
| No single clear bottleneck | Document balanced profile, prioritize highest-impact area |
| Optimization breaks physics | Always verify energy conservation after changes |

---
*Roadmap created: 2026-01-17*
*Last updated: 2026-01-17*
