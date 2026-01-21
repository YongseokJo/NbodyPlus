# Building ABYSS

## Prerequisites

### Required

- **C++ compiler**: g++ 7+ or clang++ with C++11 support
- **MPI**: OpenMPI, MPICH, or Intel MPI
- **Make**: GNU Make

### Optional

- **HDF5**: For compressed output (highly recommended)
- **gfortran**: Required only for McLuster IC generation

## Quick Build

```bash
# Build everything (ABYSS + McLuster if gfortran available)
make

# Build only ABYSS
make abyss

# Build only McLuster
make mcluster
```

## Build Targets

| Target | Description |
|--------|-------------|
| `make` or `make all` | Build ABYSS and McLuster (default) |
| `make abyss` | Build only ABYSS binary |
| `make mcluster` | Build only McLuster binary |
| `make clean` | Remove all build artifacts |
| `make mcluster-clean` | Remove only McLuster artifacts |
| `make mcluster-rebuild` | Clean and rebuild McLuster |

## Environment Variables

| Variable | Description |
|----------|-------------|
| `DISABLE_MCLUSTER=1` | Skip McLuster build even if gfortran available |
| `QUIET=1` | Suppress verbose build output |

### Examples

```bash
# Build without McLuster
DISABLE_MCLUSTER=1 make

# Quiet build
QUIET=1 make

# Both
DISABLE_MCLUSTER=1 QUIET=1 make
```

## Build Output

After successful build:

```
./ABYSS              # Main simulation binary
./mcluster/mcluster_sse    # McLuster IC generator (if built)
src/mcluster -> ../mcluster/mcluster_sse  # Symlink for ABYSS access
```

## Compiler Configuration

The build system uses the following defaults in `src/Makefile`:

```makefile
CXX = mpicxx
CXXFLAGS = -O3 -std=c++11
LDFLAGS = -lhdf5
```

To customize, edit `src/Makefile` or set environment variables:

```bash
CXX=mpiicpc make    # Use Intel MPI compiler
```

## HDF5 Setup

ABYSS uses HDF5 for output. On most systems:

```bash
# Ubuntu/Debian
sudo apt install libhdf5-dev

# RHEL/CentOS
sudo yum install hdf5-devel

# macOS (Homebrew)
brew install hdf5

# HPC systems (module)
module load hdf5
```

If HDF5 is in a non-standard location:

```bash
export HDF5_DIR=/path/to/hdf5
export CPATH=$HDF5_DIR/include:$CPATH
export LIBRARY_PATH=$HDF5_DIR/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH
```

## McLuster Build

McLuster is built with gfortran and gcc for the SSE/BSE stellar evolution libraries.

### Requirements

- gfortran (Fortran compiler)
- gcc (C compiler for SSE/BSE)

### Checking gfortran

```bash
which gfortran
gfortran --version
```

If gfortran is missing, the build will skip McLuster with a warning:

```
warning: gfortran not found, skipping mcluster build
```

### Installing gfortran

```bash
# Ubuntu/Debian
sudo apt install gfortran

# RHEL/CentOS
sudo yum install gcc-gfortran

# macOS (Homebrew)
brew install gcc
```

## Verification

After building, verify the installation:

```bash
# Check ABYSS binary
./ABYSS --help 2>&1 || echo "ABYSS built (no --help flag)"

# Check McLuster (if built)
ls -la src/mcluster
./mcluster/mcluster_sse -h 2>&1 | head -5
```

## Troubleshooting

### "mcluster not found" at runtime

If ABYSS reports McLuster not found when using `[mcluster]` config:

1. Rebuild with gfortran available: `make mcluster`
2. Verify symlink exists: `ls -la src/mcluster`
3. If symlink broken, recreate: `cd src && ln -sf ../mcluster/mcluster_sse mcluster`

### HDF5 linking errors

```
undefined reference to `H5Fcreate'
```

Ensure HDF5 is installed and in library path:

```bash
ldconfig -p | grep hdf5
export LD_LIBRARY_PATH=/path/to/hdf5/lib:$LD_LIBRARY_PATH
```

### MPI compiler not found

```
mpicxx: command not found
```

Load MPI module or install MPI:

```bash
module load openmpi   # HPC systems
# or
sudo apt install libopenmpi-dev   # Ubuntu
```

## Clean Rebuild

For a complete clean rebuild:

```bash
make clean
make
```
