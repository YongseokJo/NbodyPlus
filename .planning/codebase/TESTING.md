# ABYSS Testing

## Test Infrastructure

### Test Locations

| Path | Description |
|------|-------------|
| `tests/test1/` | Primary test case with TOML config |
| `tests/test_1/`, `tests/test_2/` | Legacy test cases |
| `tests/test_10/`, `tests/test_20/` | Larger particle count tests |
| `tests/ICs/` | Initial condition files |

### Test Configuration Files

**TOML Format** (`tests/test1/config.toml`):
```toml
Filename = "nbody.dat"
StopTime = 1.0e7
OutputDirectory = "output"

[numerics]
eta = 0.01
FixNumNeighbor = 100
InitialRadius = 0.2

[output]
dtOutput = 1.0e6
Compression = true
```

**Legacy Format** (`tests/test_1/config.txt`):
```
Filename = nbody.dat
eta = 0.01
FixNumNeighbor = 100
```

## Running Tests

### Manual Test Execution

```bash
cd tests/test1
../../src/abyss.exe -c config.toml
```

### Workflow-Based Testing

```bash
# Configure workflow
source workflow/config.sh

# Build and run
workflow/bin/build.sh
workflow/bin/run.sh
```

### SLURM Execution

Tests are typically run via SLURM on HPC clusters:

```bash
sbatch workflow/bin/run.sh
```

## Analysis Tools

### Energy Conservation (`tools/analyze_energy.py`)

Verifies energy conservation (dE/E0) across simulation:

```bash
python tools/analyze_energy.py output/output.h5 --plot energy.png
```

**Outputs**:
- Table of kinetic, potential, total energy per timestep
- Max/mean |dE/E0| statistics
- Energy evolution plot

### Profiling Analysis (`tools/analyze_profiling.py`)

Analyzes performance when `PERFORMANCETRACE` is enabled.

### Run Summary (`tools/summarize_run.py`)

Generates summary statistics for a simulation run.

## Test Output Structure

### HDF5 Output Format

```
output.h5
├── Step_0/
│   ├── Time_Myr (attribute)
│   ├── Mass_Msun (dataset)
│   ├── X_pc, Y_pc, Z_pc (datasets)
│   ├── Vx_km_s, Vy_km_s, Vz_km_s (datasets)
│   └── E_binary, E_merger, E_PN (attributes)
├── Step_1/
│   └── ...
```

### Additional Output Files

| File | Description |
|------|-------------|
| `binary_output.txt` | Binary formation/evolution events |
| `merger_output.txt` | Merger events with detailed properties |
| `worker_output_N.txt` | Per-worker debug output |
| `SEVN_output.txt` | Stellar evolution events (when USE_SEVN=1) |

## Validation Criteria

### Energy Conservation

Primary validation metric is relative energy error:

```
|dE/E0| = |(E_total - E_initial) / E_initial|
```

Typical acceptable thresholds:
- Direct N-body: |dE/E0| < 10^-6 to 10^-8
- With binaries/mergers: Higher tolerance due to energy tracking

### Few-Body Events

Tracked in `binary_output.txt` and `merger_output.txt`:
- Binary formations
- Exchanges
- GW-driven mergers
- Tidal disruptions
- Stellar collisions

## Current Test Cases

### test1 (Primary)
- Small particle count
- TOML configuration
- Quick smoke test

### test_10
- 10^4 particles (from `tests/ICs/c1e4.dat`)
- Production-scale testing
- GPU performance validation

### test_20
- Larger test case
- Extended runtime

## Continuous Integration

**No formal CI pipeline currently exists.**

Typical validation workflow:
1. Build with `make` in `src/`
2. Run test case
3. Check energy conservation with `analyze_energy.py`
4. Review output files for errors

## Debug Modes

### Compile-Time Debug

```bash
# In src/Makefile, add:
CXXFLAGS += -DDEBUG
```

### CUDA Debugging

```bash
# Enable CUDA memory checking
cuda-memcheck ./abyss.exe -c config.toml

# Or cuda-gdb for interactive debugging
cuda-gdb ./abyss.exe
```

### Profiling

```bash
# Enable in workflow/config.sh
ENABLE_PROFILING="1"

# Analyze with
python tools/analyze_profiling.py
```
