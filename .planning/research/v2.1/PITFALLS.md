# Pitfalls Research: MPI Async Optimization

## Critical Pitfalls

### 1. Buffer Modification Before Completion

**What it is:** Modifying a send buffer between `MPI_Isend` and `MPI_Wait`

**Warning signs:**
- Intermittent data corruption
- Non-deterministic results
- Works with small message counts, fails with large

**Prevention:**
- Each worker has dedicated buffers
- Clear ownership semantics
- Wait before any buffer reuse

**Phase to address:** Phase 11 (Async Infrastructure)

**Recovery:** If detected, add buffer copies or wait points

### 2. Request Handle Leaks

**What it is:** Failing to wait/test on all requests, or using `MPI_Request_free` incorrectly

**Warning signs:**
- Memory growth over time
- MPI internal errors after many iterations
- "Too many pending requests" errors

**Prevention:**
- Track all outstanding requests
- Initialize to `MPI_REQUEST_NULL`
- Always wait/test before request reuse
- Never use `MPI_Request_free` on active receives

**Phase to address:** Phase 11-12 (Infrastructure + QueueScheduler)

**Recovery:** Add request tracking, audit all paths

### 3. Deadlock from Misordered Operations

**What it is:** All processes sending before receiving, filling buffers

**Warning signs:**
- Hang with all ranks in MPI_Send/Isend
- Works with few messages, hangs with many
- Different behavior with different MPI implementations

**Prevention:**
- Post receives before sends
- Use non-blocking on both sides
- Never assume MPI buffer capacity

**Phase to address:** Phase 12-13 (QueueScheduler + Integration)

**Recovery:** Reorder to post receives first

### 4. Lost Completions

**What it is:** Missing a completion notification, leaving work unprocessed

**Warning signs:**
- Tasks never complete
- Workers stuck waiting
- Particle counts mismatch

**Prevention:**
- Track expected vs completed count
- Always process Waitany/Testany results
- Re-post receives after each completion

**Phase to address:** Phase 13 (IrregularRoutines Integration)

**Recovery:** Add completion tracking assertions

### 5. CM Particle Dependency Breakage

**What it is:** Processing CM particles before their neighbors are updated

**Warning signs:**
- Energy conservation failure
- SDAR integration errors
- Physics results differ from blocking version

**Prevention:**
- Maintain CM particle special handling
- Don't process CM result until neighbors done
- Test with few-body cases specifically

**Phase to address:** Phase 13 (IrregularRoutines Integration)

**Recovery:** Defer CM processing, add dependency checks

## Performance Pitfalls

### 6. No Actual Overlap

**What it is:** Async calls but no work done between send and wait

**Warning signs:**
- No speedup vs blocking
- Profiler shows 100% wait time
- All work in async "window" is negligible

**Prevention:**
- Profile overlap effectiveness
- Add local work between sends and waits
- Use Testany to check during computation

**Phase to address:** Phase 14 (Measurement & Tuning)

**Recovery:** Restructure to create overlap opportunity

### 7. Excessive Progress Checking

**What it is:** Calling MPI_Testany too frequently, overhead exceeds benefit

**Warning signs:**
- High CPU usage with little throughput
- Testany appears hot in profiler
- More time testing than computing

**Prevention:**
- Test only when local work is exhausted
- Batch test calls (every N iterations)
- Fall back to Waitany when no local work

**Phase to address:** Phase 14 (Measurement & Tuning)

**Recovery:** Reduce test frequency, add batching

### 8. Pipeline Depth Mismatch

**What it is:** Too few or too many in-flight messages

**Too few:** Workers idle waiting for tasks
**Too many:** Memory pressure, request management overhead

**Warning signs:**
- Workers with low utilization (too few)
- Memory growth, slow request handling (too many)

**Prevention:**
- Match pipeline depth to worker count
- Profile worker idle time
- Monitor request array sizes

**Phase to address:** Phase 14 (Measurement & Tuning)

**Recovery:** Tune pipeline depth based on measurements

## Correctness Pitfalls

### 9. Non-Determinism from Processing Order

**What it is:** Results differ between runs due to completion order

**Warning signs:**
- Different results each run
- Energy conservation varies
- Hard to reproduce bugs

**Prevention:**
- Physics should be order-independent within timestep
- Use stable processing for CM dependencies
- Document any order-sensitive code

**Phase to address:** Phase 13 (Integration)

**Recovery:** Enforce ordering where physics requires it

### 10. Mixed Blocking/Non-Blocking Issues

**What it is:** Mixing async and sync calls incorrectly during transition

**Warning signs:**
- Deadlocks during incremental migration
- Some paths work, others hang
- Regression in previously working code

**Prevention:**
- Clear separation of blocking/non-blocking paths
- Test each path independently
- Maintain backward compatibility during migration

**Phase to address:** Phase 12-13 (Transition period)

**Recovery:** Isolate async code, test in controlled settings

## Pitfall Checklist by Phase

| Phase | Critical Pitfalls | Tests |
|-------|-------------------|-------|
| 11 | #1 Buffer, #2 Request handles | Unit tests for worker async |
| 12 | #3 Deadlock, #2 Request handles | QueueScheduler tests |
| 13 | #4 Lost completions, #5 CM deps, #9 Determinism | Full simulation, energy check |
| 14 | #6 No overlap, #7 Excessive testing, #8 Pipeline | Profiling, benchmarks |

## Sources

- [ENCCS Intermediate MPI - Non-blocking](https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/)
- [Cornell Virtual Workshop - Wait/Test](https://cvw.cac.cornell.edu/mpip2p/nonblocking-communication/wait-test-free)
- [MPI Request_free Issues](https://github.com/mpi-forum/mpi-forum-historic/issues/47)

---
*Research completed: 2026-01-18*
