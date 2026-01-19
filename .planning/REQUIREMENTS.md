# Requirements: ABYSS v2.2 Load Balance Profiling

**Defined:** 2026-01-18
**Core Value:** Physics correctness (energy conservation) with clean architecture for targeted optimizations

## v2.2 Requirements

Requirements for load balance profiling milestone. Each maps to roadmap phases.

### Neighbor Profiling

- [x] **NEIGH-01**: Track per-particle neighbor count during force calculation
- [x] **NEIGH-02**: Compute neighbor count statistics (min/avg/max/stddev) per interval
- [x] **NEIGH-03**: Identify outlier particles (>2σ neighbor count)
- [x] **NEIGH-04**: Correlate neighbor count with per-particle compute time

### Queue Dispatch

- [ ] **QUEUE-01**: Measure dispatch latency (time from worker completion to next task received)
- [ ] **QUEUE-02**: Track root-side dispatch overhead (time spent in assign vs waiting)
- [ ] **QUEUE-03**: Sample queue depth over time (pending tasks in queue)
- [ ] **QUEUE-04**: Detect worker starvation events (worker idle with non-empty queue)

### Worker Distribution

- [ ] **WRKR-01**: Count particles processed per worker per interval
- [ ] **WRKR-02**: Track compute time per worker per interval
- [ ] **WRKR-03**: Identify "heavy" particles (>2σ compute time)
- [ ] **WRKR-04**: Compute load balance ratio (max_worker_time / avg_worker_time)

### Particle Type Breakdown

- [ ] **PTYPE-01**: Separate timing for CM particles vs regular particles
- [ ] **PTYPE-02**: Track CM particle frequency and compute time ratio
- [ ] **PTYPE-03**: Few-body search time (separate from force calculation)
- [ ] **PTYPE-04**: Few-body initialization time
- [ ] **PTYPE-05**: Few-body integration time

### Memory Access

- [ ] **MEM-01**: Measure L1/L2/L3 cache miss rates in force loops (via PAPI or perf)
- [ ] **MEM-02**: Track memory bandwidth utilization during force calculation
- [ ] **MEM-03**: Identify memory-bound vs compute-bound phases

### Analysis Output

- [ ] **ANLYS-01**: Generate per-interval load balance summary (CSV)
- [ ] **ANLYS-02**: Generate per-worker histogram of compute times
- [ ] **ANLYS-03**: Generate neighbor count distribution histogram
- [ ] **ANLYS-04**: Produce analysis report with optimization recommendations

## v2.3+ Requirements

Deferred to future release based on v2.2 analysis findings.

### Optimization (TBD based on analysis)

- **OPT-01**: MPI batching (if queue dispatch is bottleneck)
- **OPT-02**: OpenMP hybrid parallelism (if workers underutilized)
- **OPT-03**: Work stealing (if particle distribution uneven)
- **OPT-04**: GPU irregular forces (if compute-bound)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Optimization implementation | v2.2 is analysis only; optimization deferred to v2.3+ |
| Async MPI revisit | Archived in v2.1; single-message async doesn't help |
| Algorithm changes | Profiling work is measurement, not modification |
| SDAR Group conversion | Kept as-is with exit-only sync |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| NEIGH-01 | Phase 15 | Complete |
| NEIGH-02 | Phase 15 | Complete |
| NEIGH-03 | Phase 15 | Complete |
| NEIGH-04 | Phase 15 | Complete |
| QUEUE-01 | Phase 16 | Pending |
| QUEUE-02 | Phase 16 | Pending |
| QUEUE-03 | Phase 16 | Pending |
| QUEUE-04 | Phase 16 | Pending |
| WRKR-01 | Phase 17 | Pending |
| WRKR-02 | Phase 17 | Pending |
| WRKR-03 | Phase 17 | Pending |
| WRKR-04 | Phase 17 | Pending |
| PTYPE-01 | Phase 18 | Pending |
| PTYPE-02 | Phase 18 | Pending |
| PTYPE-03 | Phase 18 | Pending |
| PTYPE-04 | Phase 18 | Pending |
| PTYPE-05 | Phase 18 | Pending |
| MEM-01 | Phase 19 | Pending |
| MEM-02 | Phase 19 | Pending |
| MEM-03 | Phase 19 | Pending |
| ANLYS-01 | Phase 20 | Pending |
| ANLYS-02 | Phase 20 | Pending |
| ANLYS-03 | Phase 20 | Pending |
| ANLYS-04 | Phase 20 | Pending |

**Coverage:**
- v2.2 requirements: 24 total
- Mapped to phases: 24
- Complete: 4 (Phase 15)
- Pending: 20

---
*Requirements defined: 2026-01-18*
*Last updated: 2026-01-18 (Phase 15 complete)*
