# Requirements: ABYSS v2.0 Performance Profiling & Optimization

**Defined:** 2026-01-17
**Core Value:** Identify and eliminate irregular force bottlenecks through systematic profiling

## v2.0 Requirements

Requirements for performance profiling infrastructure and optimization.

### Profiler Core

- [x] **PROF-01**: Profiler aggregates per-rank statistics via MPI_Reduce (min/avg/max across all ranks)
- [x] **PROF-02**: Profiler supports histogram buckets for call-time distributions
- [x] **PROF-03**: Profiler tracks work counters (particles processed, neighbor pairs evaluated)
- [x] **PROF-04**: Profiler outputs load balance report showing work distribution across ranks
- [x] **PROF-05**: Profiler can compute throughput metrics (pairs/second, particles/second)

### Irregular Force Profiling

- [x] **IRR-01**: Sub-timer captures time in main neighbor loop (lines 99-161 in compute_acceleration.cpp)
- [x] **IRR-02**: Sub-timer captures time in CM particle loop (lines 163-205)
- [x] **IRR-03**: Sub-timer captures time in 4th-order correction (lines 221-267)
- [x] **IRR-04**: Counter tracks total neighbor pairs evaluated per timestep
- [x] **IRR-05**: Per-rank irregular force timing visible in aggregated output

### MPI/Worker Profiling

- [x] **MPI-01**: Worker-side timer captures time waiting in MPI_Recv for work
- [x] **MPI-02**: Worker-side timer captures task dispatch/execution time
- [x] **MPI-03**: Worker-side timer captures MPI_Isend/Wait completion time
- [x] **MPI-04**: Queue scheduler timer captures assignQueueAuto() time
- [x] **MPI-05**: Queue scheduler timer captures runQueueAuto() time
- [x] **MPI-06**: Idle time per rank is trackable (time between task completion and next assignment)

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
| PROF-01 | Phase 8 | Complete |
| PROF-02 | Phase 8 | Complete |
| PROF-03 | Phase 8 | Complete |
| PROF-04 | Phase 8 | Complete |
| PROF-05 | Phase 8 | Complete |
| IRR-01 | Phase 8 | Complete |
| IRR-02 | Phase 8 | Complete |
| IRR-03 | Phase 8 | Complete |
| IRR-04 | Phase 8 | Complete |
| IRR-05 | Phase 8 | Complete |
| MPI-01 | Phase 8 | Complete |
| MPI-02 | Phase 8 | Complete |
| MPI-03 | Phase 8 | Complete |
| MPI-04 | Phase 8 | Complete |
| MPI-05 | Phase 8 | Complete |
| MPI-06 | Phase 8 | Complete |
| ANLZ-01 | Phase 9 | Pending |
| ANLZ-02 | Phase 9 | Pending |
| ANLZ-03 | Phase 9 | Pending |
| OPT-01 | Phase 10 | Pending |
| OPT-02 | Phase 10 | Pending |
| OPT-03 | Phase 10 | Pending |

**Coverage:**
- v2.0 requirements: 22 total
- Mapped to phases: 22
- Unmapped: 0 ✓

---
*Requirements defined: 2026-01-17*
*Last updated: 2026-01-17 Phase 8 complete*
