# Plan 23-01 Summary: Submit Profiling Jobs

**Completed:** 2026-01-19
**Status:** Complete

## Deliverables

| Deliverable | Status | Location |
|-------------|--------|----------|
| Config for 100-step run | ✓ | tests/test4/config.toml |
| Main profiling job (100 steps) | ✓ Submitted | Job ID: 6376323 |
| Variance run 1 (20 steps) | ✓ Submitted | Job ID: 6376486 |
| Variance run 2 (20 steps) | ✓ Submitted | Job ID: 6376487 |
| Job tracking file | ✓ | .planning/phases/23-deep-analysis/JOBS.md |

## Configuration

Modified `tests/test4/config.toml`:
- StopTime: 1.0e6 years (100 Myr for ~100 outputs)
- dtOutput: 1.0e4 years (10 kyr intervals)
- Profiling enabled via `--profile` flag

## Run Directories

| Run | Directory |
|-----|-----------|
| analysis-100 | workflow/runs/analysis-100_20260119_125902/ |
| variance-1 | workflow/runs/variance-1_20260119_130121/ |
| variance-2 | workflow/runs/variance-2_20260119_130127/ |

## Notes

- Jobs submitted to SLURM cluster with `--profile` flag
- Jobs run asynchronously; Plan 23-02 will collect results after completion
- Config uses adjusted parameters (StopTime=1.0e6, dtOutput=1.0e4) for ~100 output intervals

## Requirements Addressed

- ANLYS-01: Run extended profiling (submission phase complete)

---
*Plan completed: 2026-01-19*
