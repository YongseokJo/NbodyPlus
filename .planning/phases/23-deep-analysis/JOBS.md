# Phase 23: Profiling Jobs Tracking

**Submitted:** 2026-01-19
**Status:** Running

## Jobs

| Tag | Job ID | Run Directory | Steps | Status |
|-----|--------|---------------|-------|--------|
| analysis-100 | 6376323 | workflow/runs/analysis-100_20260119_125902/ | 100 | COMPLETE |
| variance-1 | 6376486 | workflow/runs/variance-1_20260119_130121/ | 200 | COMPLETE |
| variance-2 | 6376487 | workflow/runs/variance-2_20260119_130127/ | 200 | COMPLETE |

## Configuration

**Main run (analysis-100):**
- StopTime = 1.0e6 years
- dtOutput = 1.0e4 years
- Expected outputs: 100

**Variance runs:**
- Same config as main run (20 outputs each based on shorter StopTime)

## Monitoring

Check job status:
```bash
squeue -j 6376323,6376486,6376487
```

Check completion:
```bash
sacct -j 6376323,6376486,6376487 --format=JobID,State,ExitCode
```

## Expected Outputs

After completion, each run directory should contain:
- `work/output/profiling.csv` — time series data
- `work/output/profiling_*.json` — per-interval snapshots
- `summary.txt` — analysis summary

---
*Tracking file created: 2026-01-19*
