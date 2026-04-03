# ABYSS workflow — tutorial

This tutorial explains the lightweight workflow in `workflow/` for building, running, and
analyzing ABYSS runs. The workflow is intentionally small and scheduler-agnostic: it supports
`local` and `slurm` (and templates for `pbs`) out of the box and is designed to be extended.

Contents
- Overview
- Prerequisites
- Quickstart (local)
- Submitting to Slurm
- Configuration and customization
- Run directory layout and logs
- Common troubleshooting
- Advanced: adding a scheduler template

Overview
--------

The workflow provides three main step scripts and a submit wrapper:

- `workflow/bin/compile.sh` — compile ABYSS in a reproducible per-run workspace
- `workflow/bin/run.sh` — run the compiled `abyss.exe` inside the run workspace
- `workflow/bin/analyze.sh` — collect basic artifacts and produce a summary
- `workflow/bin/submit.sh` — orchestrates the above (local or scheduler submission)

Each invocation creates a run directory under `workflow/runs/` named like:

- `TAG_YYYYMMDD_HHMMSS`

It contains logs and a `work/` area with the staged test inputs and outputs.

Prerequisites
-------------

Ensure the following are available on the machine you use to run the workflow:

- `bash` (POSIX shell compatible)
- `make` and usual build toolchain (C/C++ compiler, `mpicxx`)
- Optional for GPU builds: `nvcc` and a CUDA-capable driver; `nvidia-smi` helps auto-detection
- `sbatch` / Slurm on clusters when using the `slurm` scheduler

If your cluster uses environment modules, prefer adding a `workflow/config.local.sh` to
automatically load cluster modules on `submit.sh` calls (example below).

Quickstart — local smoke run
----------------------------

1. Run a local smoke run (stages `test/test1` into a run workspace and runs one task):

```bash
workflow/bin/submit.sh --scheduler local --tag smoke_local --test-dir test/test1 --config config.toml
```

2. What happens:

- A run dir is created at `workflow/runs/smoke_local_<timestamp>`
- `compile.sh` builds `src/abyss.exe` (logs to `build.log`)
- `run.sh` stages inputs into `work/`, copies the `abyss.exe` and runs it (logs to `run.log`)
- `analyze.sh` produces a `summary.txt` in the run dir

3. Inspect run results:

```bash
ls -l workflow/runs
less workflow/runs/smoke_local_<timestamp>/build.log
less workflow/runs/smoke_local_<timestamp>/run.log
less workflow/runs/smoke_local_<timestamp>/summary.txt
```

Note: Local runs auto-disable CUDA when no GPU is detected (this avoids crashes on login nodes).

Submitting to Slurm
-------------------

To run on Slurm, call `submit.sh` with `--scheduler slurm`. Example:

```bash
workflow/bin/submit.sh --scheduler slurm --tag my_experiment --test-dir test/test1 --config config.toml
```

This renders `workflow/templates/slurm.sbatch.in` into `workflow/runs/my_experiment_<ts>/job.sbatch` and
submits it with `sbatch`. The job script calls the same `compile.sh`, `run.sh`, and `analyze.sh` in the
run directory so the pipeline is identical across schedulers.

How test cases are staged (what `run.sh` does)
---------------------------------------------

Yes — today `workflow/bin/run.sh` stages the *entire* `TEST_DIR` into the per-run work directory:

- Source: `TEST_DIR` (e.g. `test/test1/`)
- Destination: `RUN_DIR/work/`

This is done for reproducibility and safety:

- ABYSS reads multiple inputs besides the TOML config (e.g. `nbody.dat`), and the config may reference files
  by relative path.
- ABYSS writes outputs relative to the working directory (e.g. `output/`), so running in `RUN_DIR/work/`
  keeps outputs isolated per run.

In principle, the *minimum* you need is:

- the executable (`abyss.exe`)
- the config file (e.g. `config.toml`)
- every input file referenced by the config (commonly `nbody.dat`)

But copying the whole test directory avoids missing an input file and prevents writing outputs back into
the source test directory.

Configuration and customization
-------------------------------

- `workflow/config.sh` contains default settings (scheduler, Slurm account/partition, NTASKS, default test dir).
- You can create `workflow/config.local.sh` (gitignored) to override defaults for your cluster. Example `config.local.sh`:

```bash
# workflow/config.local.sh (example)
SCHEDULER=slurm
ACCOUNT=b1094
PARTITION=ciera-gpu
NTASKS=32
MPI_CANDIDATES=(/path/to/openmpi /another/path)
CUDA_CANDIDATES=(/usr/local/cuda)
```

- `submit.sh` accepts the following flags:
  - `--scheduler` — `local` or `slurm` (overrides `config.sh`)
  - `--tag` — short name appended to the run directory
  - `--test-dir` — relative path to the test case to stage
  - `--config` — the run config filename inside `TEST_DIR` (default `config.toml`)
  - `--ntasks` — override `NTASKS` for this run
  - `--summary-file` — set the summary file name (relative to run dir unless absolute)
  - `--summary-stack-file` — set the stacked TSV file (relative to repo root unless absolute)

The script exports namespaced `WF_*` overrides so the step scripts inherit the resolved settings even when run
in new processes.

Run directory layout and logs
-----------------------------

Each run directory looks like:

- `meta.txt` — resolved settings and timestamp
- `config.sh`, `config.local.sh` — snapshots of workflow-level config
- `job.sbatch` (for scheduler runs) — rendered job script
- `submit.log` — sbatch/qsub output when applicable
- `build.log` — output from `compile.sh` (stdout+stderr)
- `run.log` — output from `run.sh` (stdout+stderr)
- `work/` — staged test inputs and produced outputs
- `summary.txt` — produced by `analyze.sh` (name configurable)
- `tools_analysis.log` — tool execution log for Python analysis
- `summary_runs.tsv` — project-level stacked summary (one TSV row per run)

Use `build.log` first if a run fails to compile, then `run.log` for runtime issues.

Common troubleshooting
----------------------

- Build fails because `mpicxx` or `nvcc` not found: set `MPI_CANDIDATES`/`CUDA_CANDIDATES` in `config.local.sh` or load modules in `config.local.sh`.
- `nvcc` host-compiler errors like `cc1plus: unrecognized option '-std=c++11'`: set `CUDAHOSTCXX` to a modern `g++` (done in the workflow via `workflow/bin/common.sh` and validated at runtime).
- Missing third-party headers (e.g., `sevn.h`): edit `workflow/config.local.sh` to point `USE_SEVN=1` and provide the `SEVN` include path or disable `USE_SEVN`.
- If the program reports `Failed to get CUDA context` on local runs, either run on a GPU node (Slurm) or let the workflow disable CUDA for local runs (it does this automatically when no GPU is found).

How to rerun specific steps
---------------------------

- Re-run only the compile step for an existing run dir (useful for iterative development):

```bash
WF_SCHEDULER_OVERRIDE=local WF_TEST_DIR_OVERRIDE=test/test1 WF_RUN_CONFIG_OVERRIDE=config.toml WF_USE_CUDA_OVERRIDE=0 workflow/bin/compile.sh workflow/runs/<timestamp>_tag
```

- Re-run only the run step (after a successful build):

```bash
workflow/bin/run.sh workflow/runs/<timestamp>_tag
```

- Run analysis only:

```bash
workflow/bin/analyze.sh workflow/runs/<timestamp>_tag
```

Running analysis tools in `tools/` with your Python
---------------------------------------------------

`workflow/bin/analyze.sh` does two things:

1. Always writes `summary.txt` (log tails + quick status)
2. Optionally runs Python tools from `tools/` and writes:
  - `tools_analysis.log`
  - per-tool outputs like `analyze_profiling.out`, `analyze_energy.out`

By default, the workflow uses `python3`, but you can point it to your venv Python.
For your setup:

```bash
workflow/bin/submit.sh --scheduler local --tag smoke_local --python /gpfs/home/vjl4366/pyenv/venv/bin/python
```

Or, set it permanently in `workflow/config.local.sh`:

```bash
PYTHON=/gpfs/home/vjl4366/pyenv/venv/bin/python
```

Note: some tools need extra Python deps (e.g. `h5py`, `matplotlib`, `pandas`). If a tool fails, check
the corresponding `*.err` next to `tools_analysis.log`.

Concise performance + energy summary
------------------------------------

`analyze.sh` appends a concise summary (performance timer totals + energy conservation stats) to the
summary file using `tools/summarize_run.py`. You can configure:

- `SUMMARY_FILE` in `workflow/config.sh` (or `workflow/config.local.sh`)
- `--summary-file` flag on `submit.sh`

The summary includes architecture metadata (CPU/GPU mode, nodes, tasks, GPUs) pulled from `meta.txt`.

Stacked summary (project-level TSV)
----------------------------------

By default, `analyze.sh` also appends a row to a project-level TSV file so multiple runs can be compared
or plotted easily. The default file is:

- `summary_runs.tsv` in the repo root

You can override it via:

- `SUMMARY_STACK_FILE` in `workflow/config.sh` (or `workflow/config.local.sh`)
- `--summary-stack-file` on `submit.sh`

The TSV includes (in this order):

- `tag`
- `simulation_duration_myr`
- `total_wall_s`
- `dE_over_E0_mean`
- `dE_over_E0_std`
- `scheduler`
- `cpu_arch`
- `gpu_arch`
- `nodes`
- `ntasks`
- `gpus`
- `timestamp`
- `run_dir`

Advanced: adding a scheduler template
------------------------------------

Scheduler job scripts are template files found in `workflow/templates/`. The Slurm template is `slurm.sbatch.in` and calls the step scripts in the run directory.

To add a new scheduler:

1. Create `workflow/templates/<your-sched>.in`.
2. Call the step scripts relative to the `RUN_DIR` (the renderer will substitute `RUN_DIR` and `REPO_ROOT`).
3. Update `submit.sh` to handle the new scheduler and submission command.

Troubleshooting tips for templates

- Always render your template to a run dir and inspect the resulting job script before submitting:

```bash
workflow/bin/render.sh workflow/templates/slurm.sbatch.in /tmp/job.sbatch ACCOUNT b1094 PARTITION ciera-gpu RUN_DIR /path/to/run
less /tmp/job.sbatch
```

McLuster IC generation
----------------------

The workflow supports McLuster integration for automatic initial condition (IC) generation. When a config file contains a `[mcluster]` section, the workflow automatically handles IC generation.

**Build requirements:**

- `gfortran` compiler (for McLuster's Fortran components)
- `gcc` (for SSE/BSE stellar evolution libraries)

**How it works (two-job workflow for Slurm):**

For Slurm submissions with `[mcluster]` configs, the workflow submits **two separate jobs**:

1. **McLuster job** — runs on CPU partition with OpenMP parallelization
2. **ABYSS job** — runs on GPU partition, depends on McLuster job completion

This separation is necessary because:
- McLuster uses OpenMP (single process, multiple threads)
- ABYSS uses MPI (multiple processes)
- Running them in the same job would limit McLuster to 1 CPU core

The workflow automatically:
- Submits McLuster job first (`mcl_<tag>`)
- Submits ABYSS job with `--dependency=afterok:<mcluster_job_id>`
- ABYSS detects pre-generated IC file and skips McLuster execution

**McLuster SLURM configuration:**

Configure in `workflow/config.sh` or `workflow/config.local.sh`:

```bash
# McLuster job settings (CPU-only, OpenMP)
MCLUSTER_PARTITION="ciera-std"    # CPU partition (no GPU needed)
MCLUSTER_WALLTIME="02:00:00"      # IC generation time limit
MCLUSTER_CPUS="16"                # OpenMP threads
MCLUSTER_MEM="32G"                # Memory for large N
```

**Enabling McLuster:**

By default, `USE_MCLUSTER=1` in `workflow/config.sh`. The workflow will:

1. Build McLuster alongside ABYSS using the root Makefile
2. For Slurm: submit separate McLuster and ABYSS jobs with dependency
3. For local: stage McLuster binary and let ABYSS invoke it

**Disabling McLuster:**

```bash
# Via flag
workflow/bin/submit.sh --no-mcluster --scheduler local --tag my_run

# Or in config.local.sh
USE_MCLUSTER=0
```

**Example config with McLuster:**

```toml
# config.toml
Filename = "nbody.dat"
StopTime = 1e7
OutputDirectory = "output"

[mcluster]
N = 10000          # 10,000 stars
P = 0              # Plummer profile
R = 0.8            # Half-mass radius in pc
f = 1              # Kroupa IMF
Z = 0.02           # Solar metallicity
```

**Large N simulations (1M+ particles):**

For large simulations, the two-job workflow is essential:

```toml
[mcluster]
N = 1000000        # 1 million stars
P = 0              # Plummer profile
R = 0.8
f = 1
Z = 0.02
```

With 16 OpenMP threads, generating 1M particles takes ~6-7 minutes. Without proper
parallelization, it would take 60+ minutes.

**Generate IC only (no simulation):**

```toml
[mcluster]
N = 100000
generate_only = true   # Exit after IC generation
```

**Run directory with McLuster:**

When using the two-job workflow, the run directory contains additional files:

- `mcluster.sbatch` — McLuster job script
- `mcluster.log` — McLuster execution log
- `mcluster_stdout.log`, `mcluster_stderr.log` — SLURM output for McLuster job
- `work/mcluster_ic.txt` — raw McLuster output
- `work/mcluster_abyss.dat` — transformed IC file for ABYSS

**Troubleshooting:**

- If `gfortran` not found, McLuster build is skipped with a warning
- Config using `[mcluster]` without McLuster binary will show a warning at runtime
- Set `GFORTRAN_CANDIDATES` in `config.local.sh` if gfortran is in a non-standard location
- Check `mcluster.log` for IC generation progress and errors
- If ABYSS job shows `DependencyNeverSatisfied`, the McLuster job failed — check `mcluster_stderr.log`

Final notes and tips
--------------------

- The workflow intentionally stages each test into `work/` for reproducibility. If your tests include large artifacts, consider editing `run.sh` to use a symlink instead.
- If you're developing on a laptop or head node, prefer `--scheduler local --ntasks 1` for fast iteration.
- If you want the workflow to set up environment modules automatically, put the necessary `module load` commands in `workflow/config.local.sh`.

If you want, I can:

- add step-by-step screenshots or example outputs from a successful run;
- add a `workflow/README.quickstart.md` with a tiny end-to-end example that you can run with zero edits;
- or add `--use-cuda/--no-cuda` flags to `submit.sh` for explicit control.

— end
