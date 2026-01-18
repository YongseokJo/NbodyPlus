# Phase 7: Validation — Research

## Workflow Infrastructure

### submit.sh Workflow
The validation runs via `workflow/bin/submit.sh` which:
1. Creates timestamped run directory under `workflow/runs/`
2. Captures git info (commit, branch, tag) in `meta.txt`
3. Executes compile → run → analyze pipeline
4. Supports `--tag <name>` for identifying runs
5. Supports `--scheduler slurm|pbs|local`

### Analyze Pipeline
`workflow/bin/analyze.sh` runs after simulation:
1. Generates `summary.txt` with run metadata
2. Runs tools from `ANALYZE_TOOLS` array (configured in `workflow/config.sh`)
3. Calls `tools/summarize_run.py` to extract metrics
4. Appends row to `summary_runs.tsv` (stacked summary file)

### Key Tools
- `tools/analyze_energy.py` — computes dE/E0 from HDF5 output
- `tools/summarize_run.py` — extracts metrics, generates TSV rows
- Both use code units from `src/def.h`

## summary_runs.tsv Format

Tab-separated file with columns:
| Column | Description |
|--------|-------------|
| tag | Run identifier |
| test_name | Test directory (e.g., tests/test1) |
| simulation_duration_myr | Simulated time |
| total_wall_s | Wall clock time |
| dE_over_E0_mean | Mean energy conservation error |
| dE_over_E0_std | Std dev of dE/E0 |
| scheduler | slurm/pbs/local |
| cpu_arch | CPU model |
| gpu_arch | GPU model |
| nodes | Number of nodes |
| ntasks | MPI ranks |
| gpus | Number of GPUs |
| timestamp | Run timestamp |
| run_dir | Full path to run directory |
| git_commit | Short commit hash |
| git_commit_long | Full commit hash |
| git_branch | Branch name |
| git_tag | Tag or short commit |

## Existing Test Runs

From `summary_runs.tsv`:
| Tag | dE/E0 mean | Wall time | Git commit | Branch |
|-----|------------|-----------|------------|--------|
| baseline | 3.20665e-05 | 28.33s | 7e90c74 | stable_candidate |
| mpi_test | 1.82796e-05 | 28.01s | 665a8ad | AoS_to_SoA |
| phase4 | 2.41778e-05 | 27.96s | 2b10d70 | AoS_to_SoA |
| phase5 | 2.63974e-05 | 28.14s | c5c0967 | AoS_to_SoA |
| phase6 | 2.51783e-05 | 28.40s | b4a6116 | AoS_to_SoA |

**Key observation:** All SoA runs maintain dE/E0 < baseline (3.2e-5), meaning energy conservation is already validated. Wall times are comparable (~28s).

## Test Configurations

Default test: `tests/test1/`
- `config.toml` — 10 Myr simulation, eta=0.01, output every 1 Myr
- `nbody.dat` — Initial conditions file

Other tests available:
- `tests/test2/`, `tests/test3/`, `tests/test4/` — each with config.toml and nbody.dat

## Acceptance Criteria

| REQ-ID | Requirement | How to Verify |
|--------|-------------|---------------|
| VAL-01 | Baseline captured | ✓ Already exists in summary_runs.tsv (baseline_20260116_234947) |
| VAL-02 | Energy conservation ≤ baseline | Compare dE_over_E0_mean column |
| VAL-03 | Simulation completes | Check run.log for "Simulation Done" |
| VAL-04 | Performance improvement | Compare total_wall_s column |

## Validation Strategy

Given that intermediate phase runs (mpi_test, phase4, phase5, phase6) all pass energy conservation checks, the final validation should:

1. **Run final test** with tag `soa_final` to capture complete SoA implementation
2. **Compare against baseline** directly using summary_runs.tsv
3. **Document results** including any performance changes

The workflow already handles:
- Compilation with current codebase
- Running simulation
- Analyzing energy conservation
- Recording metrics in summary_runs.tsv

## Hardware Considerations

Baseline was run on:
- CPU: Intel Xeon Gold 6338 @ 2.00GHz
- GPU: NVIDIA A30
- 16 MPI ranks, 1 GPU

For valid comparison, final validation should use same hardware configuration.

## Research Summary

The validation infrastructure is mature and already working. Phase 7 can be executed with a single test run using `--tag soa_final`, followed by comparison against baseline metrics. All intermediate phase runs show acceptable energy conservation.
