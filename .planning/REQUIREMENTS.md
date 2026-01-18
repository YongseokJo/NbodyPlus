# Requirements: ABYSS v2.0 Performance Profiling & Optimization

**Defined:** 2026-01-17
**Core Value:** Identify and eliminate irregular force bottlenecks through systematic profiling

## v2.0 Requirements

Requirements for performance profiling infrastructure and optimization.

### Profiler Core

- [ ] **PROF-01**: Profiler aggregates per-rank statistics via MPI_Reduce (min/avg/max across all ranks)
- [ ] **PROF-02**: Profiler supports histogram buckets for call-time distributions
- [ ] **PROF-03**: Profiler tracks work counters (particles processed, neighbor pairs evaluated)
- [ ] **PROF-04**: Profiler outputs load balance report showing work distribution across ranks
- [ ] **PROF-05**: Profiler can compute throughput metrics (pairs/second, particles/second)

### Irregular Force Profiling

- [ ] **IRR-01**: Sub-timer captures time in main neighbor loop (lines 99-161 in compute_acceleration.cpp)
- [ ] **IRR-02**: Sub-timer captures time in CM particle loop (lines 163-205)
- [ ] **IRR-03**: Sub-timer captures time in 4th-order correction (lines 221-267)
- [ ] **IRR-04**: Counter tracks total neighbor pairs evaluated per timestep
- [ ] **IRR-05**: Per-rank irregular force timing visible in aggregated output

### MPI/Worker Profiling

- [ ] **MPI-01**: Worker-side timer captures time waiting in MPI_Recv for work
- [ ] **MPI-02**: Worker-side timer captures task dispatch/execution time
- [ ] **MPI-03**: Worker-side timer captures MPI_Isend/Wait completion time
- [ ] **MPI-04**: Queue scheduler timer captures assignQueueAuto() time
- [ ] **MPI-05**: Queue scheduler timer captures runQueueAuto() time
- [ ] **MPI-06**: Idle time per rank is trackable (time between task completion and next assignment)

### Analysis

- [ ] **ANLZ-01**: Profile data collected from representative simulation runs
- [ ] **ANLZ-02**: Bottleneck identified with quantitative evidence (percentage of wall time)
- [ ] **ANLZ-03**: Load imbalance quantified (max/avg ratio across ranks)

### Optimization

- [ ] **OPT-01**: Primary bottleneck addressed based on profiling findings
- [ ] **OPT-02**: Optimization validated with before/after profiling comparison
- [ ] **OPT-03**: Energy conservation verified after optimization (within baseline tolerance)

## v2.1+ Requirements

Deferred to future release. Tracked but not in current roadmap.

### Advanced Profiling

- **PROF-06**: Hardware counter integration (cache misses via PAPI/likwid)
- **PROF-07**: Flamegraph-compatible output format
- **PROF-08**: Real-time profiling dashboard

### Additional Optimizations

- **OPT-04**: GPU kernel optimization if GPU path is bottleneck
- **OPT-05**: Memory allocation optimization (pool allocators)
- **OPT-06**: NUMA-aware data placement

## Out of Scope

| Feature | Reason |
|---------|--------|
| Algorithm changes | v2.0 is infrastructure + targeted optimization, not algorithmic redesign |
| SDAR profiling internals | External library, profiler wraps calls only |
| GPU kernel rewrite | Defer unless profiling shows GPU as primary bottleneck |
| Production profiler mode | v2.0 is detailed analysis mode only |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| PROF-01 | Phase 1 | Pending |
| PROF-02 | Phase 1 | Pending |
| PROF-03 | Phase 1 | Pending |
| PROF-04 | Phase 1 | Pending |
| PROF-05 | Phase 1 | Pending |
| IRR-01 | Phase 1 | Pending |
| IRR-02 | Phase 1 | Pending |
| IRR-03 | Phase 1 | Pending |
| IRR-04 | Phase 1 | Pending |
| IRR-05 | Phase 1 | Pending |
| MPI-01 | Phase 1 | Pending |
| MPI-02 | Phase 1 | Pending |
| MPI-03 | Phase 1 | Pending |
| MPI-04 | Phase 1 | Pending |
| MPI-05 | Phase 1 | Pending |
| MPI-06 | Phase 1 | Pending |
| ANLZ-01 | Phase 2 | Pending |
| ANLZ-02 | Phase 2 | Pending |
| ANLZ-03 | Phase 2 | Pending |
| OPT-01 | Phase 3 | Pending |
| OPT-02 | Phase 3 | Pending |
| OPT-03 | Phase 3 | Pending |

**Coverage:**
- v2.0 requirements: 22 total
- Mapped to phases: 22
- Unmapped: 0 ✓

---
*Requirements defined: 2026-01-17*
*Last updated: 2026-01-17 after initial definition*
