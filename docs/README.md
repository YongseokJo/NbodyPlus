# ABYSS Documentation

ABYSS (Astrophysical Bridges for Young Star Simulations) is a high-performance N-body simulation code for star cluster dynamics with integrated initial condition generation.

## Quick Start

```bash
# Build ABYSS and McLuster
make

# Run with a configuration file
mpirun -np 4 ./ABYSS config.toml
```

## Documentation

| Document | Description |
|----------|-------------|
| [Building](building.md) | Compilation instructions and requirements |
| [Configuration](configuration.md) | All configuration options and parameters |
| [Input/Output](io.md) | File formats for initial conditions and output |

## Features

- **N-body dynamics**: Direct gravitational N-body simulation
- **MPI parallelization**: Scales across multiple compute nodes
- **Integrated IC generation**: McLuster integration for automatic initial condition creation
- **HDF5 output**: Compressed, portable output format
- **Checkpoint/restart**: Resume interrupted simulations

## Basic Workflow

### 1. Using External Initial Conditions

```toml
# config.toml
Filename = "nbody.dat"
StopTime = 1e7           # 10 Myr in years
OutputDirectory = "output"
```

```bash
mpirun -np 4 ./ABYSS config.toml
```

### 2. Using McLuster IC Generation

```toml
# config.toml
Filename = "nbody.dat"
StopTime = 1e7
OutputDirectory = "output"

[mcluster]
N = 10000          # 10,000 stars
P = 0              # Plummer profile
R = 0.8            # Half-mass radius in pc
```

```bash
mpirun -np 4 ./ABYSS config.toml
# McLuster generates IC automatically, then simulation runs
```

### 3. Generate IC Only (No Simulation)

```toml
[mcluster]
N = 100000
P = 0
generate_only = true   # Exit after IC generation
```

## Requirements

- C++11 or later
- MPI implementation (OpenMPI, MPICH, Intel MPI)
- HDF5 library (optional but recommended)
- gfortran (only for McLuster IC generation)

## License

See LICENSE file in repository root.
