# Summary: Plan 09-01

## What Was Built

Successfully ran a profiled simulation to collect timing data:

- **Test case**: test1 configuration
- **Simulation time**: 0 → 1.0 Myr (10 output intervals)
- **Output intervals**: 10 (every 0.1 Myr)
- **Profiling enabled**: PERFORMANCETRACE active

## Data Collected

**Location**: `/gpfs/home/vjl4366/pkg/ABYSS/workflow/runs/run_20260117_225215/work/output/profiling.csv`

**Contents**:
- Per-interval timing for all profiler categories
- 10 data points across simulation
- All Phase 8 timers captured (IrregularForce, QueueWait, MPI timers, etc.)

## Key Observations (Raw Data)

Per interval (~130 seconds wall time each):
- WholeRoutine: ~129-133 seconds per interval
- IrregularForce: ~68-71 seconds (dominant)
- QueueWait: 5-88 seconds (high variance, first interval highest)
- MPISend: ~43 seconds cumulative
- QueueRun: ~51 seconds
- RegularGPU: ~14-15 seconds
- RegularAdjust: ~15-19 seconds

## Commits

None (human execution step - no code changes)

## Files Changed

- `workflow/runs/run_20260117_225215/` — simulation run directory (created)
- `workflow/runs/run_20260117_225215/work/output/profiling.csv` — profiling data

## Deviations

None.

## Requirements Addressed

- ANLZ-01: Profile data collected from representative simulation runs ✓
