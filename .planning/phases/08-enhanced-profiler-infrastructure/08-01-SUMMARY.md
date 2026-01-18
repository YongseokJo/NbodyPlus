# Summary: Plan 08-01

## What Was Built

Added per-rank MPI statistics and aggregation to the profiler:

- **AggregatedStats struct**: Holds min/avg/max seconds, counts, min/max rank, load balance ratio, and throughput
- **aggregateAcrossRanks(MPI_Comm)**: Uses MPI_Reduce with MPI_MINLOC/MPI_MAXLOC to collect stats across all ranks
- **getLoadBalanceRatio(TimerID)**: Returns max/avg ratio (>1.2 indicates imbalance)
- **printAggregatedSummary()**: Outputs comprehensive cross-rank timing breakdown with load balance warnings
- **Work unit tracking**: Added work_units and interval_work_units fields to TimerStats
- **recordWork(TimerID, units)**: Method to track throughput
- **PROFILE_WORK macro**: Convenience macro for work unit recording

## Commits

- `9c9d673` feat(08-01,08-02): add MPI aggregation, histograms, and new timers

## Files Changed

- `src/profiler.h` — Added AggregatedStats, aggregation methods, work tracking

## Deviations

None. Combined with 08-02 implementation since both modify profiler.h.

## Requirements Addressed

- PROF-01: Per-rank MPI statistics (min/avg/max)
- PROF-04: Load balance metrics
- PROF-05: Throughput tracking
