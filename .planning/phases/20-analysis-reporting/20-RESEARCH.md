# Phase 20 Research: Analysis & Reporting

## Requirements

From ROADMAP.md:
- **ANLYS-01**: Generate per-interval load balance summary (CSV)
- **ANLYS-02**: Generate per-worker histogram of compute times
- **ANLYS-03**: Generate neighbor count distribution histogram
- **ANLYS-04**: Produce analysis report with optimization recommendations

## Current State Analysis

### ANLYS-01: Per-interval Load Balance Summary (CSV)

**Status: ALREADY COMPLETE**

The `writeCSV()` method (profiler.h:1191-1290) already outputs comprehensive per-interval data including:
- All timer values (ns and count)
- Phase 15: Neighbor count stats (count, min, max, mean, stddev, outliers, correlation)
- Phase 16: Queue dispatch stats (depth, latency, starvation events, assign/wait ratio)
- Phase 17: Worker distribution stats (particle counts, compute times, load balance ratio, heavy particles)
- Phase 18: Particle type breakdown (regular/CM counts and times)
- Phase 19: Cache statistics (L1D/LL misses, bandwidth, operational intensity)

This requirement is satisfied by existing infrastructure.

### ANLYS-02: Per-worker Histogram of Compute Times

**Status: NOT IMPLEMENTED**

Current state:
- `WorkerDistributionTracker` (profiler.h:396-530) tracks per-worker compute times in `compute_time_per_worker_ns_[]` vector
- Has `getLoadBalanceRatio()`, `getMaxTimeWorker()`, `getMeanComputeTimeSeconds()` accessors
- Does NOT have histogram of individual worker compute times

What's needed:
- Add histogram of per-worker total compute times to understand distribution
- Can reuse existing `Histogram` class (profiler.h:123-217) which has log-scale buckets
- Add method to `WorkerDistributionTracker::getWorkerTimeHistogram()`
- Output histogram in JSON and console summary

### ANLYS-03: Neighbor Count Distribution Histogram

**Status: ALREADY COMPLETE**

The `NeighborHistogram` class (profiler.h:698-780) already exists with:
- 15 buckets covering 0 to 50K+ neighbor counts
- `print()` method for console output
- `toJSON()` method for JSON output
- Already output in JSON via `interval_neighbor_histogram_.toJSON()` (line 1330)

This requirement is satisfied by existing infrastructure.

### ANLYS-04: Analysis Report with Optimization Recommendations

**Status: PARTIALLY COMPLETE**

Current state:
- `tools/analyze_profiling.py` provides basic recommendations based on:
  - Top timer dominance
  - MPI overhead percentage
  - QueueWait percentage
  - SkipList overhead
- Recommendations are generic (batching, async comm, load balancing)

What's needed:
- Extended Python script with load balance-specific analysis
- Analysis of Phase 15-19 metrics to identify top imbalance sources
- Specific recommendations based on:
  - Neighbor count variance (high variance → work stealing)
  - Queue dispatch bottleneck (high latency → MPI batching)
  - Worker distribution imbalance (high ratio → better scheduling)
  - CM particle overhead (high CM ratio → CM-specific optimization)
  - Memory-bound classification (memory-bound → cache optimization)
- Final ANALYSIS.md report in `.planning/phases/20-analysis-reporting/`

## Implementation Plan

### Plan 20-01: Worker Compute Time Histogram

Add histogram capability to WorkerDistributionTracker:
1. Add `Histogram worker_time_histogram_` member to track per-worker times
2. Add `buildWorkerTimeHistogram()` method that populates histogram from `compute_time_per_worker_ns_[]`
3. Add `getWorkerTimeHistogram()` accessor
4. Output histogram in `printIntervalSummary()` and JSON output

### Plan 20-02: Enhanced Python Analysis Script

Extend `tools/analyze_profiling.py`:
1. Add functions to analyze new Phase 15-19 columns
2. Add `analyze_load_balance()` function that identifies top 3 imbalance sources
3. Add specific recommendations based on metrics thresholds
4. Generate markdown report output option

### Plan 20-03: Analysis Report Generation

Create final analysis report template:
1. Add `--report` flag to analyze_profiling.py
2. Generate `.planning/phases/20-analysis-reporting/ANALYSIS.md`
3. Include:
   - Summary of findings
   - Top 3 imbalance sources (ranked by impact)
   - Recommended optimization approach for v2.3
   - Supporting data tables

## Key Files

| File | Purpose |
|------|---------|
| `src/profiler.h` | Add worker time histogram |
| `tools/analyze_profiling.py` | Extend analysis capabilities |
| `.planning/phases/20-analysis-reporting/ANALYSIS.md` | Final report |

## Dependencies

Phase 20 depends on all previous phases (15-19) being complete, which they are.

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| No profiling data available | Use test4 simulation to generate sample data |
| Python pandas not available | Keep fallback mode without pandas |
| Histogram adds overhead | Only compute at output intervals |

## Estimated Plans

3 plans in 2 waves:
- **Wave 1**: Plan 20-01 (worker histogram), Plan 20-02 (Python analysis)
- **Wave 2**: Plan 20-03 (final report - depends on running simulation)

---
*Research completed: 2026-01-19*
