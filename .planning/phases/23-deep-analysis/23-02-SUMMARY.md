# Plan 23-02 Summary: Collect Profiling Results

**Completed:** 2026-01-19
**Status:** Complete

## Deliverables

| Deliverable | Status | Location |
|-------------|--------|----------|
| Main run data (100 steps) | ✓ | workflow/runs/analysis-100_20260119_125902/ |
| Variance run 1 (200 steps) | ✓ | workflow/runs/variance-1_20260119_130121/ |
| Variance run 2 (200 steps) | ✓ | workflow/runs/variance-2_20260119_130127/ |
| Data validated | ✓ | All profiling.csv files present |

## Data Collected

| Run | Intervals | profiling.csv | JSON snapshots |
|-----|-----------|---------------|----------------|
| analysis-100 | 100 | ✓ | 99 files |
| variance-1 | 200 | ✓ | present |
| variance-2 | 200 | ✓ | present |

## Validation

All runs completed successfully:
- Job IDs: 6376323, 6376486, 6376487
- Exit status: COMPLETED
- profiling.csv files contain expected data
- Non-zero metrics confirmed for key timers

## Notes

- User ran variance runs with 200 steps instead of 20 (more data is better)
- All data suitable for analysis in Plan 23-03

## Requirements Addressed

- ANLYS-01: Run extended profiling (collection phase complete)

---
*Plan completed: 2026-01-19*
