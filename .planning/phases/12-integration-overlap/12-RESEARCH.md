# Phase 12 Research: Integration & Overlap

**Phase Goal:** Integrate async MPI into IrregularRoutines with communication-computation overlap

## Executive Summary

The current `IrregularRoutines` function uses a synchronous dispatch pattern where the root sends one task and blocks waiting for completion before processing results. Phase 11 added async MPI infrastructure (Worker::send_task_async(), QueueScheduler::waitQueueAsync(), testQueueAsync()) but IrregularRoutines doesn't use them yet.

The main integration challenge is preserving CM particle dependencies while achieving overlap. The current code already handles CM particles specially - separating them from single particles in initializeIrr() and processing them via a secondary iterator that waits for neighbor updates.

## Current IrregularRoutines Flow Analysis

### Main Loop Structure (lines 94-792)

```
while (skiplist->getFirstNode() != nullptr):
    1. Get particles for this timestep from skip list
    2. Filter inactive particles
    3. IRREGULAR FORCE:
       - initializeIrr() separates single vs CM particles
       - Loop: assign → run → wait(non-blocking) → callback
       - CM particles processed when neighbors are up-to-date
    4. IRREGULAR UPDATE: serial particle updates
    5. FEW-BODY TERMINATION: handle binary mergers/terminations
    6. FEW-BODY INITIALIZATION: form new binaries
    7. SKIP LIST UPDATE: insert particles for next timestep
    8. Delete processed node
```

### Current Dispatch Pattern (lines 156-198)

```cpp
do {
    queue_scheduler.assignQueueAuto();   // Assign queued particles to free workers
    queue_scheduler.runQueueAuto();      // Send tasks (blocking MPI_Send)
    do {
        worker = queue_scheduler.waitQueue(1);  // Non-blocking probe
        // Handle CM particle iteration (lines 162-191)
        if (worker != nullptr) { }
    } while (worker == nullptr);
    queue_scheduler.callback(worker);    // Blocking MPI_Recv for result
} while (queue_scheduler.isComplete());
```

**Key observations:**
1. `waitQueue(1)` uses MPI_Iprobe - already non-blocking
2. `runQueueAuto()` uses blocking `sendTask()` internally
3. `callback()` does blocking MPI_Recv for result
4. CM particles are processed opportunistically during wait loops

### CM Particle Handling (lines 151-191)

CM particles require special ordering:
1. They're separated into `CMPtcls` set during initializeIrr()
2. During the wait loop, code iterates through CMPtcls
3. Each CM particle needs neighbors to be up-to-date (commented out check at line 172-178)
4. CM particles get TASK_AR_INTEGRATION sent to their assigned worker
5. Result: CM particles processed interleaved with single particles

**Critical insight:** The commented neighbor check (lines 172-178) was removed, suggesting CM particles can now be sent without waiting for neighbors. This simplifies async integration.

## Phase 11 Async Infrastructure

### Worker Async Methods
- `send_task_async()`: Uses MPI_Isend with send_request handle
- `post_receive()`: Pre-posts MPI_Irecv with recv_request handle
- `callback_async()`: Handles completion without MPI_Recv (data already in buffer)

### QueueScheduler Async Methods
- `postAllReceives()`: Pre-posts receives for all workers
- `waitQueueAsync()`: Uses MPI_Waitany on active worker requests
- `testQueueAsync()`: Uses MPI_Testany for non-blocking completion check
- `callbackAsync()`: Processes completion and re-posts receive
- `runQueueAsync()`: Uses Worker::send_task_async() instead of blocking send

## Proposed Integration Strategy

### Strategy: Async Send + Waitany Pattern

Replace the current pattern with:

```cpp
queue_scheduler.postAllReceives();  // Pre-post all receives
do {
    queue_scheduler.assignQueueAuto();
    queue_scheduler.runQueueAsync();    // Async sends

    // Overlap window: do local work here
    processReadyCMParticles();

    do {
        worker = queue_scheduler.testQueueAsync();  // Non-blocking check
        if (worker != nullptr) {
            queue_scheduler.callbackAsync(worker);
        }
    } while (worker == nullptr);
    // Or use waitQueueAsync() if no local work available

} while (queue_scheduler.isComplete());
```

### Local Work for Overlap Window (OVLP-01)

Potential local work during async window:
1. **CM particle iteration** - already happens during wait loop
2. **Particle list filtering** - erase inactive particles
3. **Skip list preparation** - pre-compute next timestep values
4. **FewBody checks** - check binary interrupt states

The current CM particle processing loop is ideal local work - it iterates through CMPtcls and assigns AR integration tasks.

### MPI_Testany for Opportunistic Checks (OVLP-02)

`testQueueAsync()` already implements MPI_Testany. Use it:
1. During CM particle iteration loop
2. After each async send batch
3. Before falling back to waitQueueAsync()

## Skip List Mechanism

The skip list (`skip_list.h`) is a probabilistic data structure for efficient time-step management:
- Keys are `next_block_irr` (next irregular time block)
- Values are particle IDs grouped at that time
- `getFirstNode()` returns particles for earliest timestep
- `updateSkipList()` re-inserts particles after time advancement
- `deleteFirstNode()` removes processed timestep

**Integration impact:** Skip list operations are serial and local - no MPI interaction. They can happen during the async window as local work.

## Few-Body Subsystem Interactions

### FBTermination (lines 315-535)
- Handles binary mergers, terminations
- Updates particle lists, maps, and indices
- Contains some synchronous worker calls (lines 377-381, 654-656)
- These are single-particle operations, not bulk async candidates

### FBInitialization (lines 599-742)
- Handles new binary formation
- Contains synchronous worker calls for TASK_MAKE_GROUP, TASK_DELETE_GROUP
- These are event-driven, not in the hot path

**Recommendation:** Keep FewBody synchronous operations as-is initially. They're relatively rare events (new binary formation, termination) and not the primary MPI overhead source.

## Deadlock Risk Analysis

### Potential Deadlock Scenarios

1. **Missing receive posts**: If receives aren't posted before sends, workers may block sending results back.
   - **Mitigation:** Call `postAllReceives()` at start of async section

2. **Request reuse before completion**: Reusing MPI_Request before previous operation completes.
   - **Mitigation:** Phase 11 infrastructure handles this (wait/test before reuse)

3. **Circular wait**: Root waiting for Worker A, Worker A waiting for data from Worker B.
   - **Not applicable:** Workers don't communicate with each other, only with root

4. **Progress deadlock**: MPI library not making progress on background operations.
   - **Mitigation:** Regular calls to MPI_Testany/Waitany ensure progress

### Safe Ordering Rules

1. Post receives before sends
2. Always check/wait for send completion before reusing send buffer
3. Use non-blocking tests when doing local work
4. Fall back to blocking wait when no local work available

## Implementation Approach

### Plan 12-01: Main Loop Restructure (IRRG-01, IRRG-02)

Modify the main force calculation loop:
1. Call `postAllReceives()` before entering dispatch loop
2. Replace `runQueueAuto()` with `runQueueAsync()`
3. Replace `waitQueue(1)` + `callback()` with `testQueueAsync()` + `callbackAsync()`
4. Fall back to `waitQueueAsync()` when testQueueAsync returns nullptr

### Plan 12-02: CM Particle Dependency Handling (IRRG-03)

The current CM handling already works with async:
- CM particles are in `CMPtcls` set
- Iterator processes them during wait loop
- Each CM gets TASK_AR_INTEGRATION sent to its worker
- No blocking dependencies remain (neighbor check was commented out)

### Plan 12-03: Skip List + Local Work (IRRG-04, OVLP-01)

Identify local work opportunities:
1. CM particle iteration (existing)
2. Skip list update preparation
3. Particle state updates that don't need results

### Plan 12-04: Testany Integration (OVLP-02)

Add opportunistic completion checks:
1. After each CM particle task assignment
2. In inner wait loop
3. Before falling back to blocking wait

### Plan 12-05: Profiler Integration (OVLP-03)

Add overlap measurement:
1. Timer for async window (time between last send and first completion)
2. Counter for local work items processed during overlap
3. Timer for time spent in testQueueAsync vs waitQueueAsync

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Physics errors from ordering | Low | High | Preserve existing CM iteration pattern |
| Deadlock | Low | High | Pre-post receives, regular progress calls |
| No performance gain | Medium | Medium | Measure overlap effectiveness |
| Regression in blocking path | Low | Medium | Keep blocking methods as fallback |

## Success Metrics

1. **Correctness**: Energy conservation matches blocking version (dE/E0 < 1e-5)
2. **Overlap**: >10% of async window spent on local work
3. **No deadlocks**: Simulation completes without hangs
4. **MPI wait time reduction**: Measurable reduction in QueueWait timer

## Files to Modify

| File | Changes |
|------|---------|
| `src/irregular_routines.cpp` | Main loop restructure, async calls |
| `src/queue_scheduler.h` | Minor (already has async methods) |

## References

- Phase 11 implementation: `.planning/phases/11-async-infrastructure/`
- Profiler timers: `src/profiler.h` (MPIIsend, MPIIrecv, MPIWaitany, MPITestany)
- Worker async: `src/worker.h` (send_task_async, post_receive, callback_async)
- QueueScheduler async: `src/queue_scheduler.h` (waitQueueAsync, testQueueAsync, callbackAsync)

---
*Research completed: 2026-01-18*
