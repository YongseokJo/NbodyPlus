# ABYSS Profiling Analysis Report

**Date:** 2026-01-19
**Branch:** MPI_opt_v2
**Milestone:** v2.3 Performance Optimization
**Author:** Claude Code analysis session

---

## Executive Summary

This document consolidates all profiling analysis performed during the v2.3 milestone, including:
- Phase 23 deep analysis results
- Neighbor distribution analysis
- Queue system effectiveness evaluation
- Hardware configuration analysis
- Optimization recommendations

**Key Finding:** The current dynamic dispatch system has **28% overhead** but the workload is **too uniform** to require it. MPI batching (Phase 24) is the recommended optimization path.

---

## 1. Hardware Configuration

### Cluster: ciera-gpu (Quest HPC)

| Node | GPU Type | GPUs/Node | CPUs | Memory |
|------|----------|-----------|------|--------|
| qgpu0602 | NVIDIA L40S | 2 | 64 | ~1TB |
| qgpu0210 | NVIDIA A100 | 1 | 52 | ~192GB |
| qgpu0701-0703 | NVIDIA A30 | 4 | 64 | ~1TB |
| qgpu0704 | NVIDIA L40S | 4 | 64 | ~1TB |

### Current Run Configuration

```bash
# workflow/config.sh
PARTITION="ciera-gpu"
NODES="1"
NTASKS="16"        # 1 root + 15 workers
CPUS_PER_TASK="1"
GPUS="1"           # Single GPU for RegularGPU kernel
```

### Architecture
- **Root process (rank 0):** Dispatches particle tasks, manages queue
- **Worker processes (ranks 1-15):** Execute irregular force calculations
- **GPU:** Handles regular force calculation (RegularGPU kernel)

---

## 2. Phase 23 Deep Analysis Results

### 2.1 Time Breakdown

Data from `analysis-100_20260119_125902` (100 intervals, 100K particles):

| Component | Time (s) | Percentage | Notes |
|-----------|----------|------------|-------|
| **IrregularForce** | 11.09 | 66.4% | Main compute (neighbor loop) |
| **MPI Communication** | 4.62 | 27.7% | Send + Recv overhead |
| → MPISend | 4.03 | 24.1% | 10.3M messages |
| → MPIRecv | 0.59 | 3.5% | Completion signals |
| **QueueRun** | 4.77 | 28.6% | Dispatch loop |
| **RegularGPU** | 1.45 | 8.7% | GPU kernel |
| **RegularAdjust** | 1.41 | 8.4% | Post-GPU correction |
| **RegularUpdate** | 1.12 | 6.7% | State update |
| **QueueCallback** | 1.03 | 6.1% | Completion handling |
| **QueueAssign** | 0.90 | 5.4% | Task assignment |
| **QueueWait** | 0.50 | 3.0% | Worker starvation |

**Wall time per interval:** 16.7 seconds

### 2.2 Load Balance

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Load balance ratio | 1.009 | Excellent (1.0 = perfect) |
| Particles per worker (mean) | 689,095 | Well distributed |
| Particles per worker (stddev) | 2,603 | Very low variance |
| Worker compute time range | 963-976s | Tight range |

### 2.3 Dispatch Analysis

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Queue depth (mean) | 19,269 | Plenty of work available |
| Queue empty ratio | 0.6% | Rarely starved |
| Starvation events | 2,438,522 | HIGH — dispatch bottleneck |
| MPI messages/interval | 10,336,438 | Very high volume |
| Is dispatch bottleneck | **true** | Confirmed |

### 2.4 Amdahl's Law Analysis

**MPI Batching (Phase 24):**
- MPI overhead: p = 0.28
- Expected reduction: 67% (s = 3)
- Speedup: S = 1 / (1 - 0.28 + 0.28/3) = **1.23 (23%)**

**Dispatch Pipelining (Phase 25):**
- QueueWait overhead: p = 0.03
- Expected reduction: 80% (s = 5)
- Speedup: S = **1.02 (2%)**

**Combined:** 18-30% speedup expected

---

## 3. Neighbor Distribution Analysis

### 3.1 Raw Statistics

From `profiling_99.json`:

```
neighbor_profiling:
  count: 794,935,952 (total evaluations)
  min: 0
  max: 513
  mean: 81.81
```

### 3.2 Distribution Histogram

20 buckets, width ≈ 25.7 neighbors each:

| Bucket | Range | Count | Percent | Cumulative |
|--------|-------|-------|---------|------------|
| 0 | 0-26 | 2,753,110 | 0.3% | 0.3% |
| 1 | 26-51 | 26,109,792 | 3.3% | 3.6% |
| 2 | 51-77 | 23,492,972 | 3.0% | 6.6% |
| 3 | 77-103 | 45,631,610 | 5.7% | 12.3% |
| 4 | 103-128 | 138,758,605 | **17.5%** | 29.8% |
| 5 | 128-154 | 268,561,269 | **33.8%** | 63.6% |
| 6 | 154-180 | 279,948,981 | **35.2%** | 98.8% |
| 7 | 180-205 | 9,666,386 | 1.2% | 100.0% |
| 8+ | 205+ | 13,227 | 0.0% | 100.0% |

**Key Finding:** 86.5% of particles have 103-180 neighbors.

### 3.3 Distribution Shape

```
Neighbors:   0    50   100   150   200   250   300
            |    |    |    |    |    |    |
Frequency:  █    ███  ███  ████████████████████  █
                           ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
                           Peak: 128-180 neighbors
```

### 3.4 Statistical Summary

| Metric | Value |
|--------|-------|
| Mode | ~154 neighbors |
| Mean (histogram-derived) | ~141 neighbors |
| Standard deviation | ~32 neighbors |
| Coefficient of Variation (CV) | **23%** |

**Interpretation:** Distribution is relatively tight with low variance.

---

## 4. Queue System Effectiveness Analysis

### 4.1 The Core Question

Is dynamic dispatch worth the 28% overhead, or would static partitioning (brute-force) work just as well?

### 4.2 Static Partitioning Estimate

If each worker gets a fixed set of particles (N/P):

| Parameter | Value |
|-----------|-------|
| Particles per worker | 689,095 |
| Neighbor CV | 23% |
| Expected worker CV (CLT) | 23% / √689 ≈ **0.9%** |
| Estimated max/mean ratio | 1 + 3×0.009 ≈ **1.03** |

**Result:** Static partitioning would only cause ~3% load imbalance.

### 4.3 Cost-Benefit Comparison

| Approach | Dispatch Overhead | Load Imbalance | Net Efficiency |
|----------|-------------------|----------------|----------------|
| Dynamic (current) | 28% | 0.9% | 72% |
| Static (estimated) | 0% | 3% | 97% |
| Batched (Phase 24) | ~3% | 0.9% | 96% |

### 4.4 Estimated Wall Times

| Scenario | Calculation | Time |
|----------|-------------|------|
| Current | 16.7s | 16.7s |
| Compute only | 16.7s × 0.72 | 12.0s |
| Static + imbalance | 12.0s × 1.03 | **12.4s** |
| Batched | 12.0s × 1.03 | **12.4s** |

### 4.5 Verdict

**The neighbor distribution is too uniform to justify 28% dispatch overhead.**

However, MPI batching (Phase 24) is still the recommended path because:
1. Reduces overhead from 28% to ~3%
2. Maintains dynamic dispatch flexibility
3. Handles edge cases (CM particles, few-body)
4. Incremental change, lower risk than restructuring

---

## 5. Why Not Static Partitioning?

### Arguments For Static:
- Simpler code (no queue scheduler)
- Zero dispatch overhead
- Workload is uniform enough

### Arguments Against Static:
- CM particles need specific worker assignment
- Few-body interactions cause variable work
- Future workloads may have higher variance
- Code restructuring risk

### Recommendation:
Implement Phase 24 (MPI batching) first. If successful, consider static partitioning as a future simplification.

---

## 6. Why Not More GPUs?

### Current GPU Usage
- RegularGPU: 8.7% of wall time (1.45s)
- Single A30/L40S GPU is sufficient

### Would 4 GPUs Help?
- IrregularForce (66%): CPU-bound, multi-GPU won't help
- RegularGPU (8.7%): Already fast, diminishing returns
- MPI overhead (28%): Not GPU-related

**Verdict:** More GPUs won't address the bottleneck. Focus on MPI optimization.

---

## 7. Optimization Roadmap

### Phase 24: MPI Message Batching (Recommended)

**Goal:** Reduce 10.3M messages → ~200K messages

**Approach:**
1. Create `BatchedQueue` struct (up to 100 particles)
2. Root accumulates tasks, sends batches
3. Worker processes batch, sends single completion
4. Expected speedup: 16-27%

**Plans created:**
- 24-01: Infrastructure (BatchedQueue, MPI type)
- 24-02: Root-side batch dispatch
- 24-03: Worker-side batch processing
- 24-04: Integration and tuning

### Phase 25: Dispatch Pipelining (Optional)

**Goal:** Prefetch next task while processing current

**Expected speedup:** 2-5% additional

**Recommendation:** Implement only if Phase 24 shows remaining starvation.

### Phase 26: Verification

**Goal:** Confirm optimizations work correctly

**Checks:**
- Energy conservation (|ΔE/E| < 1e-6)
- Performance improvement measured
- No correctness regressions

---

## 8. Data Sources

### Profiling Runs

| Run | Directory | Intervals | Purpose |
|-----|-----------|-----------|---------|
| analysis-100 | workflow/runs/analysis-100_20260119_125902/ | 100 | Primary analysis |
| variance-1 | workflow/runs/variance-1_20260119_130121/ | 200 | Variance check |
| variance-2 | workflow/runs/variance-2_20260119_130127/ | 200 | Variance check |

### Key Files

| File | Contents |
|------|----------|
| `.planning/ANALYSIS.md` | Phase 23 detailed analysis |
| `.planning/NEIGHBOR-ANALYSIS.md` | Neighbor distribution analysis |
| `.planning/v2.3-ROADMAP.md` | Milestone roadmap |
| `.planning/phases/24-mpi-batching/` | Phase 24 plans |
| `tools/analyze_profiling.py` | Profiling analysis tool |
| `tools/analyze_neighbor_distribution.py` | Neighbor analysis tool |

---

## 9. Conclusions

1. **Primary bottleneck:** MPI dispatch overhead (28%), not compute

2. **Root cause:** 10.3M individual MPI messages per interval

3. **Workload characteristic:** Uniform neighbor distribution (CV=23%)

4. **Current queue system:** Effective at load balancing (ratio=1.009) but expensive

5. **Recommended optimization:** MPI batching (Phase 24)
   - Expected speedup: 16-27%
   - Low risk, incremental change
   - Maintains flexibility for edge cases

6. **Future consideration:** Static partitioning could replace dynamic dispatch if batching proves the overhead hypothesis

7. **Hardware:** Single GPU sufficient; multi-GPU won't help this workload

---

## 10. Next Steps

1. **Execute Phase 24:** `/gsd:execute-phase 24`
2. **Verify correctness:** Energy conservation check
3. **Measure speedup:** Compare wall times
4. **Document results:** Update this report

---

*Report generated: 2026-01-19*
*Analysis session: Claude Code on MPI_opt_v2 branch*
