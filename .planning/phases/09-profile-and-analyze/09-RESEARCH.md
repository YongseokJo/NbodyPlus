# Research: Phase 9 — Profile & Analyze

## Objective

Research how to collect and analyze profiling data to identify the primary bottleneck with quantitative evidence.

## Phase 8 Infrastructure Available

The enhanced profiler (Phase 8) provides:

### Timers Available

| Category | TimerID | What it measures |
|----------|---------|------------------|
| **Main** | WholeRoutine | Total wall-clock time per interval |
| **Irregular** | IrregularForce | Total irregular force calculation |
| | IrregularNeighborLoop | Main neighbor loop time |
| | IrregularCMLoop | CM particle loop time |
| | IrregularCorrection | 4th order correction time |
| | IrregularPredict | Self/neighbor prediction time |
| | IrregularPairsEvaluated | Work counter (pairs) |
| **Regular** | RegularForce | Regular force on GPU |
| | RegularGPU | GPU kernel time |
| | RegularUpdate | Regular update |
| **FewBody** | FewBodyIntegration | SDAR integration |
| **Worker** | WorkerRecvWait | MPI_Recv wait time |
| | WorkerTaskDispatch | Task execution time |
| | WorkerSendComplete | MPI_Isend/Wait time |
| **Queue** | QueueAssign | assignQueueAuto() time |
| | QueueRun | runQueueAuto() time |
| | QueueWait | Queue wait time |
| **MPI** | MPISend, MPIRecv, MPIWait, MPIBarrier | Communication primitives |

### Aggregation Features

- `aggregateAcrossRanks(MPI_Comm)` — collects min/avg/max via MPI_Reduce
- `printAggregatedSummary()` — outputs formatted table with load balance ratio
- Load balance warning when ratio > 1.2
- Throughput metrics for neighbor pairs

### Output Formats

1. **Console output** — real-time summary at each output interval
2. **CSV file** — `profiling.csv` with all timer data per interval
3. **JSON file** — `profiling_*.json` snapshots

### Analysis Tool

`tools/analyze_profiling.py` can:
- Load CSV or JSON profiling data
- Calculate time breakdown by component
- Identify top timers and percentages
- Generate optimization recommendations
- Plot timeline charts (with matplotlib)

## Test Cases Available

| Test | Path | Particles | Purpose |
|------|------|-----------|---------|
| test1 | `tests/test1/` | Small | Quick smoke test |
| test2 | `tests/test2/` | Small | Alternative config |
| test3 | `tests/test3/` | Medium | Extended test |
| test4 | `tests/test4/` | Medium | Extended test |

Default test config (`tests/test1/config.toml`):
- StopTime: 1.0e7 years
- dtOutput: 1.0e6 years (10 output intervals)
- 100 neighbors per particle

## Execution Workflow

### Build with Profiling

```bash
# workflow/config.sh already has:
ENABLE_PROFILING="1"

# Build
./build.sh --slurm --test
```

### Run with Profiler

The profiler is enabled via `PERFORMANCETRACE` compile flag (already set via `ENABLE_PROFILING=1`).

```bash
# Submit job
cd tests/test1
sbatch ../../workflow/templates/slurm.sbatch.in

# Or run manually (16 MPI ranks)
srun -n 16 ../../src/abyss.exe -c config.toml
```

### Collect Output

After simulation completes:
1. Console output contains aggregated summary
2. `output/profiling.csv` contains per-interval data
3. Slurm log file captures all console output

### Analyze Results

```bash
python tools/analyze_profiling.py tests/test1/output/
# or
python tools/analyze_profiling.py tests/test1/output/profiling.csv --plot
```

## Key Questions to Answer

### ANLZ-01: Profile Data Collection

1. Which test case to use? → test1 (quick) vs larger test for realistic workload
2. How many MPI ranks? → Match typical production (8-32 ranks)
3. How long to run? → Enough output intervals for statistical significance

### ANLZ-02: Bottleneck Identification

Primary analysis targets:
1. **Irregular vs Regular force ratio** — which dominates?
2. **Irregular force breakdown** — NeighborLoop vs CMLoop vs Correction vs Predict
3. **MPI overhead** — communication time as % of total
4. **Queue wait time** — are workers idle?

Bottleneck thresholds:
- Component >50% of wall time = primary bottleneck
- Component >30% = secondary target
- MPI overhead >15% = communication-bound

### ANLZ-03: Load Imbalance

Key metrics:
- Load balance ratio (max/avg) for each timer
- Ratio >1.2 indicates imbalance worth addressing
- Which ranks are slowest (max_rank in output)?

## Analysis Document Structure

The analysis output should include:

1. **Test Configuration**
   - Test case used
   - MPI ranks, GPU usage
   - Simulation parameters

2. **Time Breakdown**
   - Top 10 timers by wall time
   - Percentage of total for each
   - Primary bottleneck identified

3. **Irregular Force Analysis**
   - Sub-timer breakdown
   - Pairs/second throughput
   - Which sub-component dominates?

4. **Load Balance Analysis**
   - Per-timer imbalance ratios
   - Worst-offending ranks
   - Correlation with work distribution

5. **Optimization Targets**
   - Prioritized list for Phase 10
   - Expected impact estimates

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Profiler overhead skews results | Check overhead < 5% via WholeRoutine comparison |
| Small test not representative | Run additional larger test if needed |
| GPU time dominates | Document; GPU optimization deferred to v2.1 |
| Balanced profile (no clear bottleneck) | Document; pick highest-impact area |

## Dependencies

- Phase 8 complete (profiler infrastructure) ✓
- Build system with PERFORMANCETRACE enabled ✓
- Analysis tools available ✓
- Test cases available ✓

## Recommended Plan Structure

**Plan 09-01: Run Profiled Simulation**
- Build with profiling enabled
- Run test1 (or selected test case)
- Collect output and logs

**Plan 09-02: Analyze and Document**
- Parse profiling output
- Create ANALYSIS.md with findings
- Identify optimization targets for Phase 10

---
*Research completed: 2026-01-17*
