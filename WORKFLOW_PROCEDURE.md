# ABYSS Workflow (compile → run → analyze)

This repository includes a small workflow harness under `workflow/` to make it easy to:

1) select an environment (local / Slurm / PBS),
2) compile ABYSS,
3) run a chosen test case,
4) analyze logs and summarize results.

## Files

- `workflow/config.sh`
  - Main configuration (scheduler choice, resources, modules, test selection).
  - You can override locally by creating `workflow/config.local.sh` (not tracked).

- `workflow/templates/slurm.sbatch.in`
- `workflow/templates/pbs.pbs.in`
  - Job templates used for submission.
  - These are rendered into a concrete job script per run.

- `workflow/bin/submit.sh`
  - Renders a job script and submits it via `sbatch` or `qsub`.

- `workflow/bin/compile.sh`
  - Runs `make` in `src/` and writes logs into the run directory.

- `workflow/bin/run.sh`
  - Copies `src/abyss.exe` into the test directory and runs it.
  - Uses `srun` inside a Slurm allocation (recommended) or `mpirun` if configured.

- `workflow/bin/analyze.sh`
  - Greps for common failure markers and produces `summary.txt`.

## One-time setup

Make scripts executable:

```bash
chmod +x workflow/bin/*.sh
```

If you need local overrides (highly recommended), create:

```bash
cp workflow/config.sh workflow/config.local.sh
```

Then edit `workflow/config.local.sh`.

## Configuration

Edit `workflow/config.sh` (or better: `workflow/config.local.sh`). Key settings:

- `SCHEDULER`: `slurm` | `pbs` | `local`
- `ACCOUNT`, `PARTITION` (Slurm), `QUEUE` (PBS)
- `NODES`, `NTASKS`, `CPUS_PER_TASK`, `GPUS`, `WALLTIME`
- `TEST_DIR`: default `test/test1`
- `TEST_CONFIG`: default `config.txt`
- `MODULES`: list of module names to load (optional)
- `CXX`: MPI compiler wrapper, usually `mpicxx`
- `USE_CUDA`: `auto` | `1` | `0`

## Submit a job (Slurm/PBS)

```bash
bash workflow/bin/submit.sh run
```

Artifacts are placed under:

- `workflow/runs/<timestamp>/` (compile logs, run logs, summary)
- `logs/` (scheduler stdout/stderr)

## What to look at after a run

In the run directory:

- `compile_stdout.txt`, `compile_stderr.txt`
- `run_stdout.txt`, `run_stderr.txt`
- `summary.txt`

If the compile fails with `mpi.h: No such file or directory`, it means the job environment does not have MPI headers available. Fix by:

- setting `MODULES=(...)` to load the correct OpenMPI/IntelMPI module, and/or
- setting `CXX=mpicxx` (and ensuring it exists in PATH inside the job).

## Running locally (no scheduler)

This harness focuses on scheduler runs. For local runs you can still use:

```bash
SCHEDULER=local bash workflow/bin/compile.sh workflow/runs/local
SCHEDULER=local bash workflow/bin/run.sh     workflow/runs/local
SCHEDULER=local bash workflow/bin/analyze.sh workflow/runs/local
```
