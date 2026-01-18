# Requirements: ABYSS v2.1 MPI Communication Optimization

**Defined:** 2026-01-18
**Core Value:** Physics correctness (energy conservation) with clean architecture for targeted optimizations

## v2.1 Requirements

Requirements for MPI async optimization. Each maps to roadmap phases.

### Async Infrastructure

- [ ] **ASYNC-01**: Worker struct has MPI_Request arrays for send/recv operations
- [ ] **ASYNC-02**: Worker has dedicated result buffer (no buffer sharing)
- [ ] **ASYNC-03**: Worker has `send_task_async()` using MPI_Isend
- [ ] **ASYNC-04**: Worker has `post_receive()` using MPI_Irecv
- [ ] **ASYNC-05**: Request handles initialized to MPI_REQUEST_NULL
- [ ] **ASYNC-06**: Request reuse only after wait/test completion

### QueueScheduler Async

- [ ] **QSCH-01**: QueueScheduler manages recv_requests array for all workers
- [ ] **QSCH-02**: `waitQueueAsync()` uses MPI_Waitany for first completion
- [ ] **QSCH-03**: `postAllReceives()` pre-posts receives for all workers
- [ ] **QSCH-04**: Completion processing re-posts receive for worker
- [ ] **QSCH-05**: Backward compatibility maintained (blocking methods retained)

### IrregularRoutines Integration

- [ ] **IRRG-01**: Main loop restructured for async pattern
- [ ] **IRRG-02**: Multiple tasks sent before waiting for completions
- [ ] **IRRG-03**: CM particle dependencies preserved (correct ordering)
- [ ] **IRRG-04**: Skip list updates work with async completions
- [ ] **IRRG-05**: Few-body termination/initialization unaffected

### Communication-Computation Overlap

- [ ] **OVLP-01**: Local work performed between sends and waits
- [ ] **OVLP-02**: MPI_Testany used for opportunistic completion checks
- [ ] **OVLP-03**: Profiler measures overlap effectiveness

### Profiling

- [ ] **PROF-01**: Timer for MPI_Isend operations
- [ ] **PROF-02**: Timer for MPI_Irecv operations
- [ ] **PROF-03**: Timer for MPI_Waitany operations
- [ ] **PROF-04**: Timer for MPI_Testany operations
- [ ] **PROF-05**: Overlap time measured (work done during async window)
- [ ] **PROF-06**: Before/after comparison documented

### Validation

- [ ] **VALD-01**: Energy conservation maintained (dE/E0 within tolerance)
- [ ] **VALD-02**: Physics results match blocking version
- [ ] **VALD-03**: CM particle handling verified with few-body test cases
- [ ] **VALD-04**: No message ordering issues detected

### Batching Research

- [ ] **RSCH-01**: CM particle dependency handling strategy documented
- [ ] **RSCH-02**: Load balancing approach for variable work documented
- [ ] **RSCH-03**: Callback mechanism changes for batch completion documented
- [ ] **RSCH-04**: Concrete implementation recommendation for v2.2
- [ ] **RSCH-05**: Expected impact quantified (message reduction, speedup)

## v2.2+ Requirements

Deferred to future release. Tracked but not in current roadmap.

### Batching Implementation

- **BTCH-01**: Dependency-aware batching for simple particles
- **BTCH-02**: Batch message format (multiple particles per message)
- **BTCH-03**: Batch completion callbacks
- **BTCH-04**: Adaptive batch sizing based on work estimates

### Further Optimization

- **SIMD-01**: SIMD gather intrinsics for neighbor access
- **GPU-01**: GPU irregular force kernel

## Out of Scope

| Feature | Reason |
|---------|--------|
| Batching implementation | Research only for v2.1, deferred to v2.2 |
| GPU irregular forces | Major architectural change, deferred |
| Persistent requests | Variable patterns, overhead not justified |
| One-sided communication (RMA) | Different programming model, too invasive |
| Worker-side async | Explore only if root-side shows benefit first |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| ASYNC-01 | Phase 11 | Complete |
| ASYNC-02 | Phase 11 | Complete |
| ASYNC-03 | Phase 11 | Complete |
| ASYNC-04 | Phase 11 | Complete |
| ASYNC-05 | Phase 11 | Complete |
| ASYNC-06 | Phase 11 | Complete |
| QSCH-01 | Phase 11 | Complete |
| QSCH-02 | Phase 11 | Complete |
| QSCH-03 | Phase 11 | Complete |
| QSCH-04 | Phase 11 | Complete |
| QSCH-05 | Phase 11 | Complete |
| IRRG-01 | Phase 12 | Pending |
| IRRG-02 | Phase 12 | Pending |
| IRRG-03 | Phase 12 | Pending |
| IRRG-04 | Phase 12 | Pending |
| IRRG-05 | Phase 12 | Pending |
| OVLP-01 | Phase 12 | Pending |
| OVLP-02 | Phase 12 | Pending |
| OVLP-03 | Phase 12 | Pending |
| PROF-01 | Phase 11 | Complete |
| PROF-02 | Phase 11 | Complete |
| PROF-03 | Phase 11 | Complete |
| PROF-04 | Phase 11 | Complete |
| PROF-05 | Phase 13 | Pending |
| PROF-06 | Phase 13 | Pending |
| VALD-01 | Phase 13 | Pending |
| VALD-02 | Phase 13 | Pending |
| VALD-03 | Phase 13 | Pending |
| VALD-04 | Phase 13 | Pending |
| RSCH-01 | Phase 14 | Pending |
| RSCH-02 | Phase 14 | Pending |
| RSCH-03 | Phase 14 | Pending |
| RSCH-04 | Phase 14 | Pending |
| RSCH-05 | Phase 14 | Pending |

**Coverage:**
- v2.1 requirements: 32 total
- Mapped to phases: 32
- Unmapped: 0 ✓

---
*Requirements defined: 2026-01-18*
*Last updated: 2026-01-18 — Phase 11 requirements marked Complete*
