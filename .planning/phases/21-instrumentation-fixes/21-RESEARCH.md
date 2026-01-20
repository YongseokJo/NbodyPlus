# Phase 21: Instrumentation Fixes - Research

**Researched:** 2026-01-19
**Status:** Complete

## Executive Summary

The profiling instrumentation has a **fundamental architecture issue**: data is collected on WORKERS (where `compute_acceleration_irr()` executes) but profiling output is generated exclusively from ROOT's Profiler instance. There's no MPI aggregation of profiling data from workers to root.

All four issues share a common root cause: the distributed nature of the simulation vs. the centralized profiling output.

---

## Issue 1: Neighbor Count = 0

### Observation
```
neighbor_profiling.count = 0  (expected: ~281,463 based on IrregularForce calls)
```

### Code Path Analysis

**Instrumentation location:** `src/Particle/compute_acceleration.cpp:297`
```cpp
PROFILE_NEIGHBOR_TIME(total_neighbors, particle_compute_ns);
```

**Macro definition:** `src/profiler.h:2037`
```cpp
#define PROFILE_NEIGHBOR_TIME(count, time_ns) Profiler::instance().recordNeighborWithTime(count, time_ns)
```

**Recording method:** `src/profiler.h:1455-1460`
```cpp
void recordNeighborWithTime(long long neighbor_count, long long compute_time_ns) {
    recordNeighborCount(neighbor_count);
    neighbor_time_correlation_.update(...);
    interval_neighbor_time_correlation_.update(...);
}
```

### Root Cause

**Workers record, ROOT outputs:**
- `compute_acceleration_irr()` runs on WORKER processes (MPI ranks 1-N)
- Each worker has its own `Profiler::instance()` singleton
- Workers call `recordNeighborWithTime()` which updates their LOCAL profiler stats
- Profiling output (`dumpToCSV`, `dumpToJSON`) runs only on ROOT (rank 0)
- ROOT's profiler never received any neighbor data → count = 0

Confirmed by comment in `profiler.h:1815`:
```cpp
os << "\n--- Neighbor Count Statistics (Root rank) ---\n";
```

### Fix Strategy

Add MPI aggregation of neighbor statistics from workers to root before output.

---

## Issue 2: Worker Compute Time = 0

### Observation
```
worker_compute_time = 0.000s  (expected: ~120s distributed across 15 workers)
```

### Code Path Analysis

**Instrumentation location:** `src/Particle/compute_acceleration.cpp:301-319`
```cpp
#ifdef PERFORMANCETRACE
{
    auto& profiler = Profiler::instance();
    int worker_rank = profiler.getCurrentParticleWorkerRank();
    if (worker_rank > 0) {
        profiler.recordWorkerComputeTime(worker_rank, particle_compute_ns);
        ...
    }
}
#endif
```

**getCurrentParticleWorkerRank() definition:** `src/profiler.h:1558`
```cpp
int getCurrentParticleWorkerRank() const { return current_particle_worker_rank_; }
```

**Where it's set:** `src/profiler.h:1539-1542`
```cpp
void recordWorkerAssignment(int worker_rank) {
    worker_distribution_.recordAssignment(worker_rank);
    interval_worker_distribution_.recordAssignment(worker_rank);
    current_particle_worker_rank_ = worker_rank;  // Set on ROOT only!
}
```

**Called from:** `src/queue_scheduler.h:63`
```cpp
PROFILE_WORKER_ASSIGNMENT((*worker)->MyRank);  // Called on ROOT
```

### Root Cause

**Two-layer problem:**

1. **Worker rank context never set on workers:**
   - `recordWorkerAssignment()` is called on ROOT (in `queue_scheduler.h`)
   - Workers never call this, so their `current_particle_worker_rank_ = 0`
   - The check `if (worker_rank > 0)` always fails on workers
   - `recordWorkerComputeTime()` is never called

2. **Even if called, data would still be on workers:**
   - Same aggregation issue as Issue 1 — data on workers, output on root

### Fix Strategy

1. Workers should use `my_rank` directly (global variable) instead of `getCurrentParticleWorkerRank()`
2. Add MPI aggregation of worker distribution stats

---

## Issue 3: CM Particle Count = 0

### Observation
```
cm_count = 0  (may be correct if no binaries in test)
```

### Code Path Analysis

**Instrumentation location:** `src/Particle/compute_acceleration.cpp:317`
```cpp
// Phase 18: Track CM vs regular particle
profiler.recordParticleType(this->is_cm_particle, particle_compute_ns);
```

Note: This is OUTSIDE the `if (worker_rank > 0)` block, so it should be called.

### Root Cause

**Two possible causes:**

1. **MPI aggregation issue (same as Issues 1 & 2):**
   - `recordParticleType()` IS called on workers
   - But data stays in worker's local Profiler instance
   - ROOT never aggregates this data

2. **Test case may have no CM particles:**
   - `test4` configuration may not include primordial binaries
   - If no CM particles exist, count = 0 is correct

### Verification Needed

Check if `test4` includes binary systems:
```bash
grep -i "binary\|primordial" tests/test4/config.toml
```

### Fix Strategy

1. Add MPI aggregation for particle type stats
2. If test case has no binaries, create or use a binary-containing test case for verification

---

## Issue 4: Cache Counters Unavailable

### Observation
```
counters_available = false
```

### Code Path Analysis

**Initialization:** `src/profiler.h` (PROFILE_CACHE_INIT)
```cpp
void initCacheCounters() {
    cache_counters_available_ = false;  // Default to unavailable
#ifdef __linux__
    // Try to open perf_event file descriptors
    struct perf_event_attr pe;
    // ...
    fd = syscall(SYS_perf_event_open, &pe, 0, -1, -1, 0);
    if (fd == -1) {
        // Failed - counters not available
        return;
    }
    // ...
#endif
}
```

### Root Cause

**Permission restrictions on HPC systems:**
- `perf_event_open()` syscall requires specific permissions
- HPC clusters typically restrict this syscall for security
- Default kernel setting: `/proc/sys/kernel/perf_event_paranoid` >= 2
- Without root or CAP_PERFMON capability, syscall fails

### Fix Strategy Options

**Option A: Silent skip (current behavior is almost this)**
- Already returns gracefully when unavailable
- Just ensure output clearly indicates "unavailable" vs "0 misses"

**Option B: Estimation via timing**
- Estimate memory bandwidth from force loop timing
- Calculate operational intensity from known FLOPs per particle
- Classify as memory-bound or compute-bound based on roofline

**Option C: Sample-based estimation**
- Use RDTSC or high-resolution timing around small code regions
- Infer cache behavior from timing variance

**Recommendation:** Option A (silent skip with clear messaging) is sufficient. Cache profiling is "nice to have" — the critical metrics are neighbor count, worker compute time, and CM tracking.

---

## Fix Dependencies

```
Issue 1 (Neighbor count)  ─┐
Issue 2 (Worker compute)  ─┼── All need MPI aggregation
Issue 3 (CM count)        ─┘
Issue 4 (Cache counters)  ─── Independent (fallback mechanism)
```

Recommended fix order:
1. **Fix the worker rank context issue** (Issue 2) — quick fix, enables local recording
2. **Add MPI aggregation** — fixes Issues 1, 2, 3 together
3. **Improve cache fallback messaging** (Issue 4) — minor enhancement

---

## Implementation Plan

### Plan 1: Fix Worker Rank Context

**Location:** `src/Particle/compute_acceleration.cpp:301-319`

**Change:**
```cpp
// Before:
int worker_rank = profiler.getCurrentParticleWorkerRank();
if (worker_rank > 0) {

// After:
int worker_rank = my_rank;  // Use global MPI rank directly
if (worker_rank > 0) {  // Root is rank 0, workers are 1-N
```

This fixes the immediate issue of workers not recording data. However, data still only exists on workers.

### Plan 2: Add MPI Aggregation for Profiling Data

**Location:** `src/profiler.h` — new method `aggregateFromWorkers()`

**Design:**
1. At profiling output time (before `dumpToCSV`/`dumpToJSON`), call aggregation
2. Workers send their stats to root via `MPI_Reduce` or `MPI_Gather`
3. Root merges received data into its local stats
4. Output uses aggregated data

**Data to aggregate:**
- `OnlineStats` (neighbor counts) — use parallel merge formula
- `WorkerDistributionTracker` — aggregate per-worker arrays
- `ParticleTypeStats` — simple sums
- `Histogram` — element-wise sum of buckets

### Plan 3: Improve Cache Fallback

**Location:** `src/profiler.h` — cache statistics output

**Change:**
- When counters unavailable, output clear message: "Hardware counters not available (perf_event restricted)"
- Optionally estimate from timing if desired

---

## Verification Method

After implementing fixes:

1. **Build with profiling:**
   ```bash
   ./build.sh --slurm --test
   ```

2. **Run profiling job:**
   ```bash
   workflow/bin/submit.sh --test-dir tests/test4 --tag fix-instrumentation
   ```

3. **Expected results:**
   - `neighbor_profiling.count > 0` (should be ~281,463)
   - `worker_compute_time > 0` for each worker
   - `cm_count` matches actual CM particles in simulation
   - Cache stats show clear message if unavailable

---

## Key Files to Modify

| File | Changes |
|------|---------|
| `src/Particle/compute_acceleration.cpp` | Use `my_rank` instead of `getCurrentParticleWorkerRank()` |
| `src/profiler.h` | Add MPI aggregation method, modify output to use aggregated data |
| `src/irregular_routines.cpp` | Call aggregation before profiling output |

---

*Research completed: 2026-01-19*
*Issues diagnosed: 4/4*
