# Summary: Plan 22-03 JSON Output Cleanup

**Status:** Completed
**Executed:** 2026-01-19

## What Was Done

### Task 1: Add schema_version and summary to JSON ✓

Added to `Profiler::writeJSON()`:
- `schema_version: "2.3"` at top of JSON output
- `summary` section with key metrics:
  - `wall_time_s`: from WholeRoutine timer
  - `throughput_particles_per_s`: computed from particle counts / wall time
  - `load_balance_ratio`: from worker distribution (aggregated or local)
  - `primary_bottleneck`: computed by `identifyPrimaryBottleneck()` helper

Added helper methods:
- `computeAggregatedLoadBalance()` — computes load balance from aggregated worker times
- `identifyPrimaryBottleneck()` — identifies bottleneck based on metrics (dispatch starvation, MPI communication, queue scheduling, load imbalance, or compute-bound)

### Task 2: Verbosity flag — SKIPPED

Skipped as it adds complexity without clear immediate benefit. The simplified output (Task 3) provides enough information.

### Task 3: Filter zero timers from output ✓

Modified timer output loop to skip timers with zero interval data:
```cpp
if (s.interval_total_ns == 0 && s.interval_count == 0) {
    continue;
}
```

Also simplified timer output to include only essential fields:
- `interval_ns`
- `count`
- `mean_ns`

### Task 4: Update analyze_profiling.py ✓

Updated `analyze_json()` function to:
1. Read and display `schema_version`
2. Display `summary` section when available
3. Handle both old format (full timer fields) and new format (simplified, filtered)
4. Gracefully handle missing keys for backwards compatibility

## Files Modified

| File | Changes |
|------|---------|
| `src/profiler.h` | Added schema_version, summary section, zero-timer filtering, helper methods |
| `tools/analyze_profiling.py` | Updated JSON parsing for new format |

## JSON Output Changes

**Before (all timers, verbose):**
```json
{
  "step": 10,
  "timers": {
    "WholeRoutine": {"total_ns": ..., "interval_ns": ..., "count": ..., ...},
    "IrregularTotal": {"total_ns": 0, "interval_ns": 0, ...},  // zero-value
    ...
  }
}
```

**After (filtered, with summary):**
```json
{
  "schema_version": "2.3",
  "step": 10,
  "summary": {
    "wall_time_s": 36.0,
    "throughput_particles_per_s": 69000,
    "load_balance_ratio": 1.19,
    "primary_bottleneck": "dispatch starvation"
  },
  "timers": {
    "WholeRoutine": {"interval_ns": ..., "count": ..., "mean_ns": ...},
    // zero-value timers omitted
    ...
  }
}
```

## Commits

- (pending) feat(22-03): add schema version and summary to JSON output

## Verification

Build and runtime testing required to verify:
1. JSON output includes schema_version
2. Summary section populated correctly
3. Zero-value timers filtered
4. analyze_profiling.py handles both formats
