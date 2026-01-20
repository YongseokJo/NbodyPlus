# Neighbor Count Distribution Analysis

**Generated:** 2026-01-19
**Data source:** workflow/runs/analysis-100_20260119_125902/work/output/profiling_99.json

---

## Raw Data

From the profiling JSON:

```
neighbor_profiling:
  count: 794,935,952 (total neighbor evaluations)
  min: 0
  max: 513
  mean: 81.81
  histogram: [2753110, 26109792, 23492972, 45631610, 138758605, 268561269, 279948981, 9666386, 13227, 0, ...]
```

## Histogram Analysis

20 buckets covering 0-513 neighbors (bucket width ≈ 25.7):

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
| 8 | 205-231 | 13,227 | 0.0% | 100.0% |

**Key finding:** 86.5% of particles have 103-180 neighbors (buckets 4-6).

## Distribution Shape

```
Neighbors:   0    50   100   150   200   250   300   350   400   450   500
            |    |    |    |    |    |    |    |    |    |    |
Count:      █    ███  ███  █████████████████████████████  █
            └────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘
                           ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
                           Most particles here (103-180)
```

## Statistical Analysis

**From histogram:**
- Mode: ~154 neighbors (bucket 5-6 boundary)
- Mean: ~141 neighbors (histogram-derived)
- Standard deviation: ~32 neighbors (histogram-derived)
- Coefficient of Variation (CV): **23%** (stddev/mean)

**Interpretation:**
- Distribution is **relatively tight** (CV=23% is low)
- Most work is concentrated in a narrow range
- No extreme outliers (max 513 is rare: 13,227 samples out of 795M)

## Queue System Effectiveness

### Current System (Dynamic Dispatch)

From profiling:
- **Load balance ratio:** 1.009 (excellent, 1.0 = perfect)
- **Particles per worker:** mean=689,095, stddev=2,603 (CV=0.4%)
- **Starvation events:** 2,438,522 per interval
- **Dispatch overhead:** 28% of wall time

### Estimated Static Partitioning

If we used brute-force static partitioning (each worker gets N/P particles):

**Expected load imbalance:**
- With 689K particles per worker and CV=23% neighbor variance
- By Central Limit Theorem: worker CV ≈ 23% / √(689K/1000) ≈ 0.9%
- Expected max/mean ratio: 1 + 3×0.009 ≈ **1.03**

This means static partitioning would only cause ~3% load imbalance!

### Cost-Benefit Comparison

| Metric | Dynamic (Current) | Static (Estimated) |
|--------|-------------------|-------------------|
| Load balance ratio | 1.009 | ~1.03 |
| Dispatch overhead | 28% | 0% |
| Net efficiency | 72% | 97% |

**Estimated wall times (per interval):**
- Current: 16.7s
- Compute only: 16.7s × (1 - 0.28) = 12.0s
- Static with imbalance: 12.0s × 1.03 = **12.4s**
- Current with dispatch: 12.0s + 4.7s = **16.7s**

## Verdict

**⚠️ STATIC PARTITIONING MAY BE 25% FASTER**

- Dispatch overhead (28%) >> Load imbalance cost (3%)
- The neighbor distribution is too uniform to justify dynamic dispatch
- With 100K particles, random static assignment averages out variance

However, this analysis assumes neighbor count ∝ compute time. If there are other sources of variance (CM particles, few-body interactions), dynamic dispatch may still be valuable.

## Recommendations

### Option 1: Hybrid Approach (Recommended)
- Use **batched dynamic dispatch** (Phase 24)
- Batch 50-100 particles per message
- Get 90% of dispatch overhead reduction while keeping load balance flexibility
- Expected result: ~20% speedup

### Option 2: Static Partitioning Test
- Implement static work assignment as experiment
- Each worker gets particles `[i*N/P, (i+1)*N/P)`
- Compare wall time and energy conservation
- If faster with same correctness, simplify code permanently

### Option 3: Keep Current + Batch
- Current system provides safety margin for edge cases
- Batching (Phase 24) will reduce overhead from 28% to ~3%
- No risk of regression on unusual particle distributions

## Conclusion

The neighbor count distribution is **too uniform** to justify the 28% dispatch overhead. However:

1. **Phase 24 (MPI Batching) is still the right approach** because:
   - It reduces overhead to ~3% while keeping dynamic flexibility
   - No code restructuring needed (incremental change)
   - Handles edge cases (CM particles, few-body)

2. **Static partitioning could be a future simplification** after batching proves the overhead hypothesis.

---

*Analysis completed: 2026-01-19*
