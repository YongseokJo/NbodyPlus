# Phase 23: Deep Analysis - Research

**Researched:** 2026-01-19
**Status:** Complete

## Executive Summary

Phase 23 is an **analysis-only** phase — no code changes. The goal is to run extended profiling, analyze where time goes, and produce a priority matrix for optimizations. All tools and infrastructure are already in place from v2.2.

## 1. Existing Infrastructure

### Profiler Output (src/profiler.h)

The profiler produces two output formats:
- **JSON snapshots** (`profiling_step_NNNN.json`) — per-interval detailed metrics
- **CSV time series** (`profiling.csv`) — all intervals in tabular form

**JSON schema (v2.3):**
```json
{
  "schema_version": "2.3",
  "step": 10,
  "sim_time_myr": 10.0,
  "summary": {
    "wall_time_s": 36.0,
    "throughput_particles_per_s": 69000,
    "load_balance_ratio": 1.19,
    "primary_bottleneck": "MPI dispatch"
  },
  "timers": { /* non-zero timers with interval_ns, count, mean_ns */ },
  "worker_distribution": { /* per-worker compute times */ },
  "neighbor_profiling": { /* neighbor count stats */ },
  "particle_type_breakdown": { /* CM vs regular timing */ }
}
```

### Analysis Tool (tools/analyze_profiling.py)

Already supports:
- CSV loading with pandas
- JSON snapshot analysis
- Time breakdown by component (sorted by total time)
- Load balance analysis (Phase 15-19 metrics)
- Ranked findings with recommendations
- Markdown report generation (`--report` flag)
- Optional matplotlib plotting (`--plot` flag)

**Usage:**
```bash
python tools/analyze_profiling.py <output_dir>
python tools/analyze_profiling.py profiling.csv --report ANALYSIS.md
python tools/analyze_profiling.py profiling_*.json --plot
```

### Workflow Scripts (workflow/bin/*)

**Submit jobs:**
```bash
workflow/bin/submit.sh --scheduler slurm --tag analysis --profile
```

**Key options:**
- `--scheduler slurm` — Submit to SLURM cluster
- `--profile` — Enable profiling (`-DPERFORMANCETRACE`)
- `--tag <name>` — Name for run directory
- `--test-dir <path>` — Override test directory
- `--ntasks <N>` — MPI rank count

**Run directory structure:**
```
workflow/runs/<tag>_YYYYMMDD_HHMMSS/
├── job.sbatch       # Slurm job script
├── meta.txt         # Run metadata
├── config.sh        # Workflow config snapshot
├── output/          # Simulation output
│   ├── profiling.csv
│   ├── profiling_step_*.json
│   └── *.h5
└── summary.txt      # Analysis summary
```

### Test Configurations

**Available ICs:**
- `tests/test1/nbody.dat` — Small test case
- `tests/test4/nbody.dat` — Current default
- `tests/ICs/c1e4.dat` — 10K particles
- `tests/ICs/c1e5.dat` — 100K particles (for scaling)

**Config parameters (config.toml):**
```toml
StopTime = 3.0e5           # Simulation end time (years)
[output]
dtOutput = 1.0e5           # Output interval (years)
```

To run 100 steps: Set `StopTime` and `dtOutput` so `StopTime/dtOutput ≈ 100`.

## 2. Running Extended Profiling

### Main Run (100 steps)

1. Create config for 100 steps:
   ```toml
   StopTime = 1.0e7       # 10 Myr
   dtOutput = 1.0e5       # 0.1 Myr intervals → 100 outputs
   ```

2. Submit:
   ```bash
   workflow/bin/submit.sh --scheduler slurm --tag analysis-100 \
       --test-dir tests/test4 --profile
   ```

### Variance Runs (2 × 20 steps)

Submit two additional runs with different random seeds or ICs:
```bash
workflow/bin/submit.sh --tag variance-1 --test-dir tests/test4 --profile
workflow/bin/submit.sh --tag variance-2 --test-dir tests/test4 --profile
```

### Scaling Runs (100K particles)

User will provide 100K ICs. Submit:
```bash
workflow/bin/submit.sh --tag scale-100k-1 --test-dir <100k-ic-dir> --profile
workflow/bin/submit.sh --tag scale-100k-2 --test-dir <100k-ic-dir> --profile
```

## 3. Time Breakdown Analysis

### Categories to Report

**High-level (3-5 buckets):**
- **Compute:** IrregularForce + RegularGPU + FewBody*
- **Communication:** MPISend + MPIRecv + QueueAssign + QueueWait
- **Idle:** Worker starvation time, QueueWait beyond threshold
- **Data Structures:** SkipList* + UpdateNextRegTime
- **I/O:** FileWrite

**Detailed appendix:**
All 30+ timers sorted by percentage.

### Calculating Percentages

From JSON:
```python
whole_time = data['timers']['WholeRoutine']['interval_ns']
for timer, stats in data['timers'].items():
    pct = 100.0 * stats['interval_ns'] / whole_time
```

## 4. Amdahl's Law for Priority Matrix

### Formula

Speedup = 1 / (1 - p + p/s)

Where:
- p = fraction of time in the optimized component
- s = speedup factor for that component

### Example Calculations

**MPI Batching (Phase 24):**
- Current MPI overhead: 32% (from STATE.md)
- If batching reduces MPI time by 50%: p=0.32, s=2
- Speedup = 1 / (1 - 0.32 + 0.32/2) = 1 / 0.84 = **1.19 (19% faster)**

**Dispatch Pipelining (Phase 25):**
- Worker starvation events: 1.09M
- If pipelining eliminates 80% of starvation: additional 5-10% gain
- Combined with batching: **25-30% potential**

### Threshold Check

Minimum threshold: 15%+ speedup to pursue.
- MPI batching: 19% → **PASS**
- Dispatch pipelining: 5-10% alone, but synergistic → **CONDITIONAL**

## 5. Deliverables Checklist

| Deliverable | Location | Format |
|-------------|----------|--------|
| Extended profiling data | workflow/runs/analysis-*/ | JSON/CSV |
| Time breakdown summary | .planning/ANALYSIS.md | Markdown tables |
| Time breakdown appendix | .planning/ANALYSIS.md | Detailed timer list |
| Scaling analysis | .planning/ANALYSIS.md | Wall time vs particles |
| Priority matrix | .planning/ANALYSIS.md | Impact-ranked table |
| Charts (optional) | .planning/analysis-charts/ | PNG via matplotlib |

## 6. Key Metrics to Extract

From prior Phase 21.5 results:

| Metric | Value | Source |
|--------|-------|--------|
| Wall time | 36s | WholeRoutine |
| IrregularForce | 16.9s (47%) | Timers |
| MPI overhead | 11.4s (32%) | MPISend + MPIRecv + Queue* |
| Load balance ratio | 1.19 | worker_distribution |
| Starvation events | 1.09M | dispatch_analysis |
| MPI messages | 25M | dispatch_analysis |

**Goal for Phase 23:** Replicate with 100 steps, add variance bounds, confirm scaling.

## 7. Report Structure

```markdown
# ABYSS Performance Analysis Report

## Executive Summary
- Key findings (3 bullets)
- Recommended actions (go/no-go for Phases 24-25)

## Time Breakdown
### High-Level
| Category | Time (s) | Percentage |
### Detailed Appendix
[Full timer list]

## Load Balance Analysis
- Worker compute time distribution
- Starvation events analysis
- Load balance ratio over time

## Scaling Analysis
- Wall time vs particle count
- Scaling efficiency

## Optimization Priority Matrix
| Optimization | Expected Speedup | Threshold |
```

---

*Research complete: 2026-01-19*
