# Roadmap: ABYSS v2.1 MPI Communication Optimization

**Created:** 2026-01-18
**Milestone:** v2.1
**Goal:** Reduce MPI overhead in irregular force communication through async operations

## Phase Overview

| Phase | Name | Requirements | Focus | Status |
|-------|------|--------------|-------|--------|
| 11 | Async Infrastructure | 15 | Worker, QueueScheduler, profiler timers | ✓ Complete |
| 12 | Integration & Overlap | 8 | IrregularRoutines, overlap implementation | ✓ Complete |
| 13 | Validation & Measurement | 6 | Correctness, performance comparison | Pending |
| 14 | Batching Research | 5 | Document strategy for v2.2 | Pending |

**Total:** 4 phases, 34 requirements

---

## Phase 11: Async Infrastructure

**Goal:** Add non-blocking MPI infrastructure to Worker and QueueScheduler

**Requirements:**
- ASYNC-01: Worker struct has MPI_Request arrays for send/recv operations
- ASYNC-02: Worker has dedicated result buffer (no buffer sharing)
- ASYNC-03: Worker has `send_task_async()` using MPI_Isend
- ASYNC-04: Worker has `post_receive()` using MPI_Irecv
- ASYNC-05: Request handles initialized to MPI_REQUEST_NULL
- ASYNC-06: Request reuse only after wait/test completion
- QSCH-01: QueueScheduler manages recv_requests array for all workers
- QSCH-02: `waitQueueAsync()` uses MPI_Waitany for first completion
- QSCH-03: `postAllReceives()` pre-posts receives for all workers
- QSCH-04: Completion processing re-posts receive for worker
- QSCH-05: Backward compatibility maintained (blocking methods retained)
- PROF-01: Timer for MPI_Isend operations
- PROF-02: Timer for MPI_Irecv operations
- PROF-03: Timer for MPI_Waitany operations
- PROF-04: Timer for MPI_Testany operations

**Success Criteria:**
1. Worker can send tasks asynchronously without blocking
2. QueueScheduler can process completions via MPI_Waitany
3. All async methods have profiler instrumentation
4. Blocking methods still work (no regression)
5. Unit tests pass for new async methods

**Key Files:**
- `src/worker.h`
- `src/queue_scheduler.h`
- `src/profiler.h`

---

## Phase 12: Integration & Overlap

**Goal:** Integrate async MPI into IrregularRoutines with communication-computation overlap

**Requirements:**
- IRRG-01: Main loop restructured for async pattern
- IRRG-02: Multiple tasks sent before waiting for completions
- IRRG-03: CM particle dependencies preserved (correct ordering)
- IRRG-04: Skip list updates work with async completions
- IRRG-05: Few-body termination/initialization unaffected
- OVLP-01: Local work performed between sends and waits
- OVLP-02: MPI_Testany used for opportunistic completion checks
- OVLP-03: Profiler measures overlap effectiveness

**Success Criteria:**
1. Irregular force loop uses async send/recv
2. Multiple workers have tasks in flight simultaneously
3. CM particles processed correctly (no physics errors)
4. Overlap achieved between sends and waits
5. Simulation completes without hangs or deadlocks

**Key Files:**
- `src/irregular_routines.cpp`
- `src/queue_scheduler.h`

**Dependencies:** Phase 11 (async infrastructure)

---

## Phase 13: Validation & Measurement

**Goal:** Verify correctness and measure performance improvement

**Requirements:**
- VALD-01: Energy conservation maintained (dE/E0 within tolerance)
- VALD-02: Physics results match blocking version
- VALD-03: CM particle handling verified with few-body test cases
- VALD-04: No message ordering issues detected
- PROF-05: Overlap time measured (work done during async window)
- PROF-06: Before/after comparison documented

**Success Criteria:**
1. Energy conservation: dE/E0 < 1e-5 (matches v2.0 baseline)
2. MPI wait time reduced compared to blocking version
3. Overlap effectiveness quantified (% of time doing useful work)
4. No new physics bugs introduced
5. Performance comparison documented in analysis report

**Validation Workflow:**
- Run: `workflow/bin/submit.sh --tag async-v21 --scheduler slurm`
- Compare: `summary_runs.tsv` vs v2.0 baseline

**Dependencies:** Phase 12 (integration)

---

## Phase 14: Batching Research

**Goal:** Research and document batching strategies for v2.2

**Requirements:**
- RSCH-01: CM particle dependency handling strategy documented
- RSCH-02: Load balancing approach for variable work documented
- RSCH-03: Callback mechanism changes for batch completion documented
- RSCH-04: Concrete implementation recommendation for v2.2
- RSCH-05: Expected impact quantified (message reduction, speedup)

**Success Criteria:**
1. Clear strategy for handling CM dependencies with batching
2. Load balancing approach defined (work-weighted or adaptive)
3. Callback mechanism redesign documented
4. Implementation phases outlined for v2.2
5. Expected impact: 10-30% speedup, 100-200x message reduction

**Deliverable:** `.planning/research/v2.1/BATCHING-FINAL.md`

**Dependencies:** Phase 13 (async baseline established)

---

## Risk Mitigation

| Risk | Mitigation | Phase |
|------|------------|-------|
| Buffer corruption from early reuse | Dedicated buffers per worker, strict lifecycle | 11 |
| Deadlock from misordered operations | Post receives before sends | 11 |
| CM dependency breakage | Preserve existing CM handling logic | 12 |
| No actual overlap benefit | Profile overlap, identify local work | 12-13 |
| Performance regression | A/B comparison with blocking version | 13 |

## Definition of Done

Milestone v2.1 is complete when:
- [ ] All 4 phases completed
- [ ] 32 requirements satisfied
- [ ] Energy conservation validated
- [ ] MPI wait time reduced (measurable improvement)
- [ ] Batching research documented with v2.2 recommendation
- [ ] Code committed to `MPI_async` branch
- [ ] Ready for merge to `abyss-stable`

---
*Roadmap created: 2026-01-18*
