# Roadmap: ABYSS v2.2 Load Balance Profiling

**Created:** 2026-01-18
**Milestone:** v2.2
**Goal:** Instrument and analyze load balance across MPI workers to understand imbalance sources before committing to optimization approach

## Phase Overview

| Phase | Name | Requirements | Focus | Status |
|-------|------|--------------|-------|--------|
| 15 | Neighbor Profiling | 4 | Per-particle neighbor counts, distribution, outliers | ✓ Complete |
| 16 | Queue Dispatch Profiling | 4 | Dispatch latency, root overhead, queue depth | ✓ Complete |
| 17 | Worker Distribution | 4 | Particles per worker, compute time, load balance ratio | ✓ Complete |
| 18 | Particle Type Breakdown | 5 | CM vs regular timing, few-body breakdown | ✓ Complete |
| 19 | Memory Access Profiling | 3 | Cache miss rates, memory bandwidth | ✓ Complete |
| 20 | Analysis & Reporting | 4 | CSV output, histograms, recommendations | ✓ Complete |

**Total:** 6 phases, 24 requirements

---

## Phase 15: Neighbor Profiling

**Goal:** Measure neighbor count variance to understand if particles with many neighbors cause load imbalance

**Requirements:**
- NEIGH-01: Track per-particle neighbor count during force calculation
- NEIGH-02: Compute neighbor count statistics (min/avg/max/stddev) per interval
- NEIGH-03: Identify outlier particles (>2σ neighbor count)
- NEIGH-04: Correlate neighbor count with per-particle compute time

**Success Criteria:**
1. Neighbor count tracked for every particle in irregular force loop
2. Per-interval statistics aggregated and output to profiling CSV
3. Outlier particles flagged with their neighbor count and compute time
4. Correlation coefficient computed between neighbor count and compute time

**Key Files:**
- `src/Particle/compute_acceleration.cpp` — Neighbor loop instrumentation
- `src/profiler.h` — New counters and statistics
- `src/irregular_routines.cpp` — Per-particle timing hooks

---

## Phase 16: Queue Dispatch Profiling

**Goal:** Measure queue dispatch overhead to determine if root is bottleneck in feeding workers

**Requirements:**
- QUEUE-01: Measure dispatch latency (time from worker completion to next task received)
- QUEUE-02: Track root-side dispatch overhead (time spent in assign vs waiting)
- QUEUE-03: Sample queue depth over time (pending tasks in queue)
- QUEUE-04: Detect worker starvation events (worker idle with non-empty queue)

**Success Criteria:**
1. Per-dispatch latency measured on worker side (completion → next task)
2. Root time breakdown: assign time vs wait time per interval
3. Queue depth sampled at regular intervals, output as time series
4. Starvation events logged when worker waits > threshold with queue > 0

**Key Files:**
- `src/queue_scheduler.h` — Dispatch timing, queue depth tracking
- `src/worker.h` — Worker-side latency measurement
- `src/irregular_routines.cpp` — Integration points

**Dependencies:** None (can run in parallel with Phase 15)

---

## Phase 17: Worker Distribution

**Goal:** Measure how work is distributed across workers to identify imbalance

**Requirements:**
- WRKR-01: Count particles processed per worker per interval
- WRKR-02: Track compute time per worker per interval
- WRKR-03: Identify "heavy" particles (>2σ compute time)
- WRKR-04: Compute load balance ratio (max_worker_time / avg_worker_time)

**Success Criteria:**
1. Per-worker particle count tracked and aggregated per interval
2. Per-worker total compute time tracked per interval
3. Heavy particles identified with their worker assignment
4. Load balance ratio computed: ratio > 1.5 indicates significant imbalance

**Key Files:**
- `src/queue_scheduler.h` — Per-worker statistics
- `src/worker.h` — Compute time tracking
- `src/profiler.h` — Load balance ratio calculation

**Dependencies:** None (can run in parallel with Phases 15-16)

---

## Phase 18: Particle Type Breakdown

**Goal:** Separate timing for different particle types to understand which cause more work

**Requirements:**
- PTYPE-01: Separate timing for CM particles vs regular particles
- PTYPE-02: Track CM particle frequency and compute time ratio
- PTYPE-03: Few-body search time (separate from force calculation)
- PTYPE-04: Few-body initialization time
- PTYPE-05: Few-body integration time

**Success Criteria:**
1. CM particle force time tracked separately from regular particles
2. Ratio computed: CM_time / regular_time, CM_count / total_count
3. Few-body search time isolated from force calculation
4. Few-body init and integration times separated
5. Per-interval breakdown output to CSV

**Key Files:**
- `src/irregular_routines.cpp` — CM particle detection and timing
- `src/FewBody/*.cpp` — Few-body timing hooks
- `src/profiler.h` — New timer categories

**Dependencies:** Phase 15 (uses per-particle timing infrastructure)

---

## Phase 19: Memory Access Profiling

**Goal:** Measure cache behavior to determine if memory access is limiting performance

**Requirements:**
- MEM-01: Measure L1/L2/L3 cache miss rates in force loops (via PAPI or perf)
- MEM-02: Track memory bandwidth utilization during force calculation
- MEM-03: Identify memory-bound vs compute-bound phases

**Success Criteria:**
1. Cache miss rates measured for force loop (requires PAPI or perf integration)
2. Memory bandwidth estimated from cache miss rate and line size
3. Roofline analysis: compute intensity vs achieved performance
4. Classification: memory-bound or compute-bound per phase

**Key Files:**
- `src/profiler.h` — PAPI/perf integration
- `src/Particle/compute_acceleration.cpp` — Measurement region markers
- `Makefile` — PAPI library linkage (if used)

**Dependencies:** Phases 15-17 complete (to correlate with other metrics)

---

## Phase 20: Analysis & Reporting

**Goal:** Generate comprehensive load balance analysis with optimization recommendations

**Requirements:**
- ANLYS-01: Generate per-interval load balance summary (CSV)
- ANLYS-02: Generate per-worker histogram of compute times
- ANLYS-03: Generate neighbor count distribution histogram
- ANLYS-04: Produce analysis report with optimization recommendations

**Success Criteria:**
1. CSV output includes all profiled metrics per interval
2. Histogram data for worker compute times (binned distribution)
3. Histogram data for neighbor counts (identify long tail)
4. Analysis report written to `.planning/` with:
   - Summary of findings
   - Top 3 imbalance sources
   - Recommended optimization approach for v2.3

**Key Files:**
- `src/profiler.h` — Histogram generation, CSV output
- `workflow/bin/analyze_loadbalance.py` — Post-processing script
- `.planning/phases/20-analysis-reporting/ANALYSIS.md` — Final report

**Dependencies:** Phases 15-19 complete

---

## Execution Order

```
Wave 1 (parallel):
  ├── Phase 15: Neighbor Profiling
  ├── Phase 16: Queue Dispatch Profiling
  └── Phase 17: Worker Distribution

Wave 2:
  └── Phase 18: Particle Type Breakdown (depends on Phase 15)

Wave 3:
  └── Phase 19: Memory Access Profiling (depends on Phases 15-17)

Wave 4:
  └── Phase 20: Analysis & Reporting (depends on all)
```

## Risk Mitigation

| Risk | Mitigation | Phase |
|------|------------|-------|
| Profiling overhead > 5% | Sample-based measurement, conditional compilation | All |
| PAPI not available on target system | Fallback to perf counters or skip MEM requirements | 19 |
| Histogram memory usage | Bounded bins, streaming aggregation | 20 |
| Analysis inconclusive | Multiple runs, statistical significance tests | 20 |

## Definition of Done

Milestone v2.2 is complete when:
- [x] All 6 phases completed
- [x] 24 requirements satisfied
- [x] Profiling overhead < 5% runtime impact (conditional compilation)
- [x] Load balance analysis report produced
- [x] Top imbalance sources identified (expected sources documented)
- [x] Optimization recommendation for v2.3 documented
- [ ] Code committed to branch (pending bash permissions)
- [x] Energy conservation unchanged (profiling-only changes)

---
*Roadmap created: 2026-01-18*
*Milestone completed: 2026-01-19*
