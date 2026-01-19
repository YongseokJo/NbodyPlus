# Research: Phase 17 — Worker Distribution

## Overview

Phase 17 instruments worker-level metrics to understand how work is distributed across MPI workers and identify sources of load imbalance.

## Current Architecture Analysis

### Worker Task Flow

**Root side (queue_scheduler.h):**
1. `initializeIrr()` → Populates queue with particle IDs
2. `assignQueueAuto()` → Assigns tasks from `_queue_list` to `_FreeWorkers`
3. `runQueueAuto()` → Calls `worker->runQueue()` → MPI_Send
4. `waitQueue()` → MPI_Probe for completion
5. `callback()` → MPI_Recv, moves worker to free pool

**Worker side (worker_routines.cpp):**
1. `WorkerRecvWait` → MPI_Recv for next task
2. `WorkerTaskDispatch` → Execute task (compute_acceleration_irr, etc.)
3. `WorkerSendComplete` → MPI_Isend completion signal

### Existing Profiling Infrastructure

**From Phase 15 (reuse):**
- `OnlineStats` class — Welford's algorithm for running mean/variance
- `NeighborHistogram` — Distribution tracking
- `PROFILE_NEIGHBOR_TIME(count, time_ns)` macro

**From Phase 16 (reuse):**
- `QueueDepthTracker` — Queue depth sampling
- `PROFILE_STARVATION_EVENT()`, `PROFILE_DISPATCH_LATENCY()` macros

**Existing worker timers (Phase 8):**
- `TimerID::WorkerRecvWait` — Time waiting for next task
- `TimerID::WorkerTaskDispatch` — Task execution time
- `TimerID::WorkerSendComplete` — Completion notification time
- `TimerID::WorkerIdleTime` — Cumulative idle time

### Key Insight: What's Missing

Currently tracked:
- Per-interval aggregate time per timer (all workers combined on root output)
- Per-particle neighbor count and compute time

NOT tracked:
- Per-worker particle count
- Per-worker total compute time
- Which worker processed which particle
- Heavy particle identification with worker assignment
- Load balance ratio (max/avg worker time)

## Requirements Analysis

### WRKR-01: Particles per Worker per Interval

**Goal:** Count how many particles each worker processes.

**Implementation:**
- Root tracks assignment count per worker in `assignQueueAuto()`
- Use array indexed by worker rank: `int particles_per_worker[num_workers+1]`
- Reset at interval boundary

**Location:** `queue_scheduler.h` → `assignQueueAuto()` and `callback()`

### WRKR-02: Compute Time per Worker per Interval

**Goal:** Track total compute time per worker.

**Challenge:** Root doesn't know individual task duration on workers.

**Options:**
1. **Worker reports time in completion message** (adds MPI overhead)
2. **Root estimates via wall-clock between dispatch and callback** (includes MPI overhead)
3. **Aggregate via MPI_Reduce at interval end** (cleanest, minimal overhead)

**Recommended:** Option 3 — Aggregate worker-side `WorkerTaskDispatch` time at interval end via MPI_Reduce.

**Implementation:**
- Worker accumulates `WorkerTaskDispatch` time locally
- At output interval, root calls MPI_Reduce to get per-worker totals
- Store in `double worker_compute_time[num_workers+1]`

### WRKR-03: Heavy Particle Identification

**Goal:** Identify particles with >2σ compute time and their worker.

**Implementation:**
- Already have `PROFILE_NEIGHBOR_TIME(count, time_ns)` recording per-particle time
- Extend to track particle ID and assigned worker rank
- Need to associate particle with worker at dispatch time

**Data structure:**
```cpp
struct HeavyParticleInfo {
    int particle_id;
    int worker_rank;
    long long compute_time_ns;
    long long neighbor_count;
};
std::vector<HeavyParticleInfo> interval_heavy_particles_;
```

**Threshold:** Use `OnlineStats::isOutlier(time_ns, 2.0)` from Phase 15

### WRKR-04: Load Balance Ratio

**Goal:** Compute max_worker_time / avg_worker_time per interval.

**Implementation:**
- After MPI_Reduce gathers per-worker times, compute:
  - `avg = sum(worker_times) / num_workers`
  - `max = max(worker_times)`
  - `ratio = max / avg`
- Ratio > 1.5 indicates significant imbalance

**Location:** `Profiler::computeWorkerLoadBalance()` called after aggregation

## Integration Strategy

### New Data Structures in profiler.h

```cpp
// Per-worker statistics tracker (Phase 17)
class WorkerDistributionTracker {
public:
    void recordAssignment(int worker_rank);
    void recordCompletion(int worker_rank, long long compute_time_ns);
    void recordHeavyParticle(int particle_id, int worker_rank,
                            long long compute_time_ns, long long neighbor_count);
    void reset();

    // Accessors
    int getParticleCount(int worker_rank) const;
    double getComputeTime(int worker_rank) const;
    double getLoadBalanceRatio() const;
    const std::vector<HeavyParticleInfo>& getHeavyParticles() const;

private:
    std::vector<int> particles_per_worker_;
    std::vector<double> compute_time_per_worker_;  // seconds
    std::vector<HeavyParticleInfo> heavy_particles_;
    OnlineStats particle_time_stats_;  // for outlier detection
};
```

### New Macros

```cpp
#define PROFILE_WORKER_ASSIGNMENT(worker_rank) \
    Profiler::instance().recordWorkerAssignment(worker_rank)

#define PROFILE_WORKER_COMPLETION(worker_rank, particle_id, compute_ns) \
    Profiler::instance().recordWorkerCompletion(worker_rank, particle_id, compute_ns)
```

### Integration Points

1. **queue_scheduler.h:assignQueueAuto()**
   - Add `PROFILE_WORKER_ASSIGNMENT(worker->rank)` when assigning task

2. **queue_scheduler.h:callback()**
   - Track which worker completed task
   - Note: We don't have compute time here — need different approach

3. **profiler.h**
   - Add `WorkerDistributionTracker` class
   - Add methods and macros
   - Integrate into output (CSV, JSON, console)

4. **worker_routines.cpp**
   - Already timing `WorkerTaskDispatch`
   - Need to accumulate per-worker for later aggregation

### Heavy Particle Detection

**Challenge:** Root doesn't see individual particle compute times directly.

**Solution:** Track on worker side, aggregate at interval end:
1. Worker records heavy particles locally during execution
2. At interval output, MPI_Gather heavy particle info to root
3. Root combines and outputs

**Alternative (simpler):** Use existing `PROFILE_NEIGHBOR_TIME` data on root:
- Phase 15 already tracks per-particle compute time
- Just need to associate with worker rank at dispatch

**Implementation:**
- In `callback()`, when particle completes, check if it's an outlier
- Store particle_id → worker_rank mapping during assignment

## Output Format

### CSV Additions (per-interval)

```
WorkerParticleCount_min,WorkerParticleCount_max,WorkerParticleCount_mean,WorkerParticleCount_stddev,
WorkerComputeTime_min_s,WorkerComputeTime_max_s,WorkerComputeTime_mean_s,
LoadBalanceRatio,HeavyParticleCount
```

### Console Summary

```
--- Worker Distribution Statistics ---
Particles per worker: min=45, max=52, mean=48.2, stddev=2.1
Compute time per worker: min=1.2s, max=1.8s, mean=1.5s
Load balance ratio: 1.20 (max on rank 7)
Heavy particles (>2σ): 3
  PID 1234 on rank 5: 2.3ms (neighbors: 15420)
  PID 5678 on rank 7: 1.8ms (neighbors: 12890)
  PID 9012 on rank 7: 1.9ms (neighbors: 13200)
WARNING: Load imbalance detected (ratio > 1.5)
```

### JSON Additions

```json
"worker_distribution": {
  "particles_per_worker": {
    "min": 45, "max": 52, "mean": 48.2, "stddev": 2.1,
    "per_worker": [48, 45, 52, 47, ...]
  },
  "compute_time_per_worker": {
    "min_s": 1.2, "max_s": 1.8, "mean_s": 1.5,
    "per_worker_s": [1.5, 1.2, 1.8, 1.4, ...]
  },
  "load_balance_ratio": 1.20,
  "max_time_worker": 7,
  "heavy_particles": [
    {"pid": 1234, "worker": 5, "time_ns": 2300000, "neighbors": 15420}
  ]
}
```

## Plan Structure

| Plan | Focus | Dependencies |
|------|-------|--------------|
| 17-01 | WorkerDistributionTracker class and macros | None |
| 17-02 | Per-worker particle count instrumentation | 17-01 |
| 17-03 | Per-worker compute time aggregation | 17-01 |
| 17-04 | Heavy particle tracking + load balance ratio + output | 17-01, 17-02, 17-03 |

## Overhead Considerations

- Per-assignment tracking: O(1) counter increment, negligible
- Worker rank lookup: O(1), already available
- Heavy particle check: Uses existing OnlineStats::isOutlier(), O(1)
- MPI_Reduce for aggregation: One collective per interval, minimal
- Heavy particle gather: Only for outliers, bounded count

**Expected overhead:** <1% additional runtime.

## Risk Areas

1. **Worker rank availability** — Queue scheduler has worker pointer, can access rank
2. **Per-particle timing on root** — Phase 15 infrastructure handles this
3. **MPI aggregation overhead** — Single MPI_Reduce per interval is acceptable
4. **Heavy particle list size** — Cap at reasonable limit (e.g., 100 per interval)

## Success Criteria

1. [ ] Per-worker particle count tracked and output to CSV
2. [ ] Per-worker compute time tracked (via MPI aggregation)
3. [ ] Heavy particles (>2σ time) identified with worker assignment
4. [ ] Load balance ratio computed and output
5. [ ] Warning when ratio > 1.5
6. [ ] Overhead < 5% of runtime
