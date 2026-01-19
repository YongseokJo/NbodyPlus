# ABYSS Profiling Analysis Report

**Run:** profiling_20260119_013613
**Date:** 2026-01-19
**Branch:** MPI_opt_v2

## Simulation Configuration

| Parameter | Value |
|-----------|-------|
| Particles | 100,000 |
| MPI Ranks | 16 (1 root + 15 workers) |
| GPU | NVIDIA A100-PCIE-40GB |
| Simulation Time | 0.0 → 0.2 Myr |
| Output Intervals | 2 (at 0.1 and 0.2 Myr) |
| Wall-clock Time | 364 seconds total |

## Executive Summary

**Primary Bottleneck: MPI Dispatch Overhead**

The profiling reveals a **dispatch-bound** simulation where the root process cannot feed workers fast enough, despite excellent load balance. Key findings:

1. **Irregular force dominates (67%)** — expected for N-body
2. **Perfect worker balance (ratio=1.0)** — scheduling works well
3. **23M starvation events** — workers idle despite queued tasks
4. **103M MPI messages** — extreme message volume causes overhead
5. **Root spends 56% assigning** — dispatch is the bottleneck

**Recommended Action:** Implement MPI message batching (OPT-01) to reduce message count by 10-100x.

## Time Breakdown

### Interval 2 (0.1 → 0.2 Myr)

| Component | Time (s) | % | Analysis |
|-----------|----------|---|----------|
| **IrregularForce** | 119.8 | 67.1% | Neighbor force calculation — primary compute |
| MPISend | 42.1 | 23.6% | Communication overhead |
| RegularAdjust | 16.6 | 9.3% | Regular timestep adjustment |
| RegularGPU | 14.4 | 8.1% | GPU force calculation |
| RegularUpdate | 10.4 | 5.8% | Regular update |
| QueueWait | 7.2 | 4.0% | Worker queue waiting |
| MPIRecv | 6.0 | 3.3% | Message receive |
| IrregularUpdate | 4.2 | 2.3% | Irregular update |
| RegularSendToGPU | 3.6 | 2.0% | GPU data transfer |
| UpdateNextRegTime | 3.2 | 1.8% | Timestep management |
| SkipListCreate | 3.1 | 1.7% | Skip list construction |
| SkipListUpdate | 1.2 | 0.7% | Skip list updates |
| FewBodyTermination | 0.75 | 0.4% | SDAR termination |
| FewBodyInitialization | 0.72 | 0.4% | SDAR initialization |

### Compute vs Communication

| Category | Time (s) | % |
|----------|----------|---|
| Compute (Force + Update) | 167.0 | 93.5% |
| Communication (MPI) | 48.1 | 26.9% |
| Queue Management | 7.2 | 4.0% |

Note: Categories overlap (MPI happens during compute).

## Queue Dispatch Analysis

### Metrics

| Metric | Interval 1 | Interval 2 |
|--------|------------|------------|
| Queue depth (mean) | 19,092 | 19,175 |
| Queue depth (max) | 100,000 | 100,000 |
| Empty ratio | 0.6% | 0.6% |
| Starvation events | 22.4M | 23.6M |
| Assign time ratio | 9.8% | **56.2%** |
| Wait time ratio | 90.2% | 43.8% |

### Diagnosis

```
WARNING: Dispatch appears to be bottleneck (assign > 50%)
WARNING: 23553254 starvation events (workers waited with non-empty queue)
```

**Interpretation:**
- Queue rarely empties (0.6%) — there's always work available
- Root spends majority of time assigning tasks, not waiting
- Workers experience "starvation" — they wait even though tasks exist
- This indicates the **single-threaded root cannot dispatch fast enough**

## Worker Distribution Analysis

### Particle Distribution

| Metric | Value |
|--------|-------|
| Total particles processed | 104,755,134 |
| Workers | 15 |
| Min per worker | 6,942,246 |
| Max per worker | 7,019,173 |
| Mean per worker | 6,983,676 |
| Standard deviation | 21,449 (0.31%) |
| **Load balance ratio** | **1.00** |

**Conclusion:** Particle distribution is **near-perfect**. The scheduler successfully balances work across workers.

### Worker Compute Time

| Metric | Value |
|--------|-------|
| Min | 0.000s |
| Max | 0.000s |
| Mean | 0.000s |

**Issue:** Worker compute time not being recorded. See "Missing Data" section.

## MPI Communication Analysis

| Timer | Time (s) | Calls | Avg (μs) | Max (μs) |
|-------|----------|-------|----------|----------|
| MPISend | 42.1 | 104,755,260 | 0.40 | 14,143 |
| MPIRecv | 6.0 | 104,755,260 | 0.06 | 1,882 |

### Message Volume Analysis

- **Messages per interval:** ~104 million
- **Messages per particle:** ~1,048 per particle per interval
- **Effective rate:** 587,000 messages/second

The per-message latency is excellent (0.4μs avg), but the **volume** (104M messages) creates aggregate overhead of 48 seconds (27% of runtime).

### Batching Opportunity

| Batch Size | Messages | Est. Overhead |
|------------|----------|---------------|
| 1 (current) | 104M | 48s |
| 10 | 10.4M | ~5s |
| 100 | 1.04M | ~0.5s |

**Potential speedup from batching:** 10-20%

## Few-Body (SDAR) Timing

| Phase | Time (s) | % of Total |
|-------|----------|------------|
| Search | 0.0 | 0% |
| Initialization | 0.72 | 0.4% |
| Integration | 0.0 | 0% |
| Termination | 0.75 | 0.4% |
| **Total** | 1.47 | 0.8% |

**Conclusion:** Few-body overhead is negligible (<1%). No primordial binaries in this test case.

## Missing Data

Several Phase 15-19 metrics show zeros, indicating instrumentation issues:

### 1. Neighbor Count Profiling (Phase 15)

```json
"neighbor_profiling": {
  "count": 0,
  "min": 0, "max": 0, "mean": 0.00,
  "histogram": {"buckets":[0,0,0,...], "total":0}
}
```

**Issue:** Neighbor counting not being triggered. The `PROFILE_NEIGHBOR_COUNT` macro may not be placed in the active code path.

### 2. Worker Compute Time (Phase 17)

```json
"compute_time": {
  "min_s": 0.0000, "max_s": 0.0000, "mean_s": 0.0000
}
```

**Issue:** Per-worker compute time not being recorded. The `recordWorkerComputeTime()` calls may be missing or in wrong location.

### 3. Particle Type Breakdown (Phase 18)

```json
"particle_type_breakdown": {
  "regular_count": 0, "cm_count": 0,
  "regular_time_s": 0.0000, "cm_time_s": 0.0000
}
```

**Issue:** CM vs regular particle tracking not active. No CM particles in this simulation, or instrumentation not enabled.

### 4. Cache Statistics (Phase 19)

```json
"cache_statistics": {
  "counters_available": false,
  "l1d_misses": 0, "ll_misses": 0
}
```

**Issue:** Hardware performance counters not available on this system. Requires `perf_event_open` permissions or running as root.

## Conclusions

### What Works Well

1. **Load balancing is excellent** — ratio of 1.00, stddev < 0.5%
2. **GPU utilization appears reasonable** — 8% in RegularGPU
3. **Per-message MPI latency is good** — 0.4μs average
4. **Few-body overhead minimal** — <1%

### What Needs Improvement

1. **MPI message volume** — 104M messages is excessive
2. **Dispatch bottleneck** — Root can't keep up with workers
3. **Starvation events** — 23M events of workers waiting

### Root Cause

The simulation is **dispatch-bound**, not compute-bound or memory-bound. The single-threaded root process becomes a bottleneck when trying to feed 15 workers with individual particle tasks.

## Recommendations for v2.3

| Priority | Optimization | Expected Impact |
|----------|-------------|-----------------|
| **1** | MPI Message Batching | 10-20% speedup |
| **2** | Prefetch/Pipelining | 5-10% speedup |
| **3** | Work Stealing | Eliminate root bottleneck |
| **4** | Fix Missing Instrumentation | Enable data-driven decisions |

### Immediate Actions

1. **Fix neighbor count instrumentation** — Verify macro placement
2. **Fix worker compute time tracking** — Add missing `recordWorkerComputeTime()` calls
3. **Investigate CM tracking** — May need test case with binaries
4. **Request perf_event access** — For cache profiling on HPC system

---

*Report generated: 2026-01-19*
*Profiling infrastructure: ABYSS v2.2*
