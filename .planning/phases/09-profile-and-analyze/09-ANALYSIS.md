# Phase 9 Analysis Report

## Test Configuration

| Parameter | Value |
|-----------|-------|
| Test case | test1 |
| Simulation time | 0 → 1.0 Myr |
| Output intervals | 10 (every 0.1 Myr) |
| MPI ranks | 16 (per workflow config) |
| GPU | Enabled |
| Profiling | PERFORMANCETRACE enabled |
| Data source | `workflow/runs/run_20260117_225215/work/output/profiling.csv` |

## Time Breakdown

Averaged across all 10 output intervals (wall-clock ~130 seconds each):

| Timer | Total (s) | Percent | Calls/interval | Description |
|-------|-----------|---------|----------------|-------------|
| **IrregularForce** | **69.5** | **53.5%** | 139,934 | Irregular force coordination |
| RegularAdjust | 16.0 | 12.3% | 1,146 | Regular particle adjustment |
| RegularGPU | 14.8 | 11.4% | 1,146 | GPU regular force kernel |
| RegularUpdate | 11.1 | 8.5% | 1,146 | Regular particle update |
| IrregularUpdate | 4.4 | 3.4% | 139,934 | Irregular particle update |
| RegularSendToGPU | 4.0 | 3.1% | 1,146 | Data transfer to GPU |
| UpdateNextRegTime | 3.7 | 2.8% | 1,146 | Next regular time calc |
| SkipListCreate | 3.4 | 2.6% | 1,146 | Skip list creation |
| SkipListUpdate | 1.2 | 0.9% | 139,934 | Skip list update |
| FewBodyTermination | 0.78 | 0.6% | 139,934 | Few-body check |
| FewBodySearch | 0.77 | 0.6% | 139,934 | Few-body search |

**Accounted time:** 99.5% of wall clock

## Primary Bottleneck

### **IrregularForce: 53.5% of wall time**

The irregular force calculation dominates execution time, taking over half of total wall-clock time per output interval.

**Evidence:**
- Consistent across all 10 intervals (68-71 seconds per interval)
- 139,934 irregular force calls per interval vs 1,146 regular calls
- Irregular forces are evaluated ~122x more frequently than regular forces

**What this timer measures:**
- Time root rank spends coordinating irregular force calculations
- Includes time waiting for workers to complete force evaluations
- Includes MPI communication overhead for distributing work and collecting results

### Irregular Force Sub-Timers

The Phase 8 sub-timers (IrregularNeighborLoop, IrregularCMLoop, IrregularCorrection, IrregularPredict) show zero in the CSV because:

1. These timers run on **worker ranks**, not the root rank
2. The CSV is written by the root rank only
3. The aggregated console output (printAggregatedSummary) would show cross-rank data

**Implication:** To see the breakdown within irregular force, need to capture the console output which includes the MPI-aggregated statistics.

## Queue/MPI Analysis

The queue scheduler and MPI timers show significant activity:

| Timer | Total (s) | Calls/interval | Notes |
|-------|-----------|----------------|-------|
| QueueRun | 51.4 | 104,882,749 | Sending tasks to workers |
| MPISend | 43.5 | 104,882,749 | MPI send operations |
| QueueCallback | 10.6 | 80,342,599 | Handling worker completions |
| QueueAssign | 8.6 | 104,882,749 | Assigning work to queue |
| MPIRecv | 6.2 | 104,882,749 | MPI receive operations |
| QueueWait | varies | 24,439,984 | Waiting for workers |

**Note:** These timers are cumulative and may overlap with IrregularForce time, as they represent the mechanism by which irregular forces are distributed.

## Load Balance Analysis

The current CSV captures root-rank timing only. Load balance across workers requires the aggregated output from `printAggregatedSummary()`.

**Observable indicators:**
- QueueWait high variance (5-88 seconds between intervals) suggests variable load
- First interval has much higher QueueWait (88s) indicating initial load balancing overhead
- Subsequent intervals stabilize around 5-7 seconds

**Recommendation:** Capture console output in future runs to get per-rank min/avg/max and load balance ratios.

## Irregular vs Regular Force Ratio

| Category | Time (s) | Percent | Calls |
|----------|----------|---------|-------|
| Irregular (force + update) | 73.9 | 56.8% | 139,934 |
| Regular (GPU + adjust + update + send) | 45.9 | 35.3% | 1,146 |
| Data structures (skip list) | 4.6 | 3.5% | - |
| Few-body | 1.6 | 1.2% | 139,934 |

**Irregular force is 1.6x more expensive than regular force despite GPU acceleration for regular forces.**

## Work Distribution

Per interval:
- ~140,000 irregular force evaluations
- ~1,100 regular force evaluations
- ~105 million MPI messages (send + recv combined)

The high message count relative to force evaluations suggests fine-grained task distribution, which may contribute to communication overhead.

## Optimization Targets for Phase 10

### Priority 1: Irregular Force Optimization (53.5% of time)

**Current state:** Irregular force dominates wall time

**Potential optimizations:**
1. **Vectorization of neighbor loop** — SIMD instructions for force accumulation
2. **Cache optimization** — Improve data locality for neighbor accesses
3. **Reduce task granularity** — Batch multiple particles per MPI message to reduce overhead
4. **Overlap communication** — Use non-blocking MPI to hide latency

**Expected impact:** 10-30% reduction in irregular force time → 5-15% overall speedup

### Priority 2: Reduce MPI Message Count

**Current state:** ~105 million MPI messages per interval for ~140,000 force evaluations

**Potential optimizations:**
1. **Batch sends** — Group multiple particle updates into single messages
2. **Reduce synchronization points** — Fewer barriers between phases

**Expected impact:** Reduce communication overhead within IrregularForce

### Priority 3: GPU Utilization for Irregular Forces

**Current state:** Regular forces use GPU (14.8s), irregular forces do not

**Potential optimization:**
1. **GPU irregular force kernel** — Port irregular force to GPU

**Expected impact:** Significant if neighbor loop is compute-bound, but requires major code changes (defer to v2.1)

## Conclusion

**Primary bottleneck identified:** IrregularForce at 53.5% of wall time

**Root cause:** High frequency of irregular force evaluations (~122x more than regular) combined with per-particle MPI communication overhead.

**Recommended Phase 10 focus:**
1. Optimize irregular force loop (vectorization, cache locality)
2. Reduce MPI message count through batching
3. Profile worker-side sub-timers to pinpoint exact hotspot within irregular force

**Data gaps:**
- Worker-side sub-timer breakdown not captured in CSV (need console output)
- Per-rank load balance ratios not available from CSV alone

---
*Analysis completed: 2026-01-17*
*Data source: workflow/runs/run_20260117_225215/work/output/profiling.csv*
