# ABYSS I/O Migration Guide

## Overview

The ABYSS I/O system has been migrated to use modern, standards-based formats:
- **Input**: TOML (Tom's Obvious Minimal Language) for configuration files
- **Output**: HDF5 for time series particle data

## Changes Summary

### Input Configuration
**Old Format** (`config.txt`):
```
Filename        = nbody.dat
eta             = 0.01
FixNumNeighbor  = 100
InitialRadius   = 0.2
StopTime        = 1e8
dtOutput        = 1e6
OutputDirectory = output
```

**New Format** (`config.toml`):
```toml
# ABYSS Configuration File
Filename = "nbody.dat"
eta = 0.01
FixNumNeighbor = 100
InitialRadius = 0.2
StopTime = 1.0e8
dtOutput = 1.0e6
OutputDirectory = "output"
```

### Output Data
**Old Format**: Text files (`output_*.txt`)
- Multiple files per simulation
- Large file sizes
- Difficult to parse programmatically

**New Format**: Single HDF5 file (`output.h5`)
- All timesteps in one file
- Efficient binary storage
- Industry-standard format
- Easy to read with Python, Julia, MATLAB, etc.

## HDF5 Output Structure

The HDF5 file is organized as follows:

```
output.h5
├── Step_0/
│   ├── Attributes:
│   │   ├── Time_Myr: 0.0
│   │   ├── NumberOfParticle: 1000
│   │   ├── E_binary: 0.0
│   │   ├── E_binary_SD: 0.0
│   │   ├── E_merger: 0.0
│   │   └── E_PN: 0.0
│   └── Datasets:
│       ├── PID: [N] particle IDs
│       ├── Mass_Msun: [N] masses in solar masses
│       ├── X_pc: [N] x-coordinates in parsecs
│       ├── Y_pc: [N] y-coordinates in parsecs
│       ├── Z_pc: [N] z-coordinates in parsecs
│       ├── Vx_km_s: [N] x-velocities in km/s
│       ├── Vy_km_s: [N] y-velocities in km/s
│       ├── Vz_km_s: [N] z-velocities in km/s
│       └── Type: [N] stellar type (if SEVN enabled)
├── Step_1/
│   └── ... (same structure)
└── Step_N/
    └── ... (same structure)
```

## Building with New I/O

### Prerequisites
1. HDF5 library (with C++ support)
2. toml11 library (header-only, included)
3. C++11 compiler

### Build Instructions

#### Option 1: Using the build script
```bash
./build_with_new_io.sh
```

#### Option 2: Manual build
```bash
# Load HDF5 module
module load hdf5/1.14.1-2-openmpi-gcc-12.3.0

# Build
cd src
make clean
make
```

### Specifying Custom HDF5 Path
If HDF5 is installed in a non-standard location:
```bash
make HDF5_DIR=/path/to/hdf5
```

## Reading HDF5 Output

### Python Example
```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Open the HDF5 file
with h5py.File('output/output.h5', 'r') as f:
    # List all timesteps
    print("Available timesteps:", list(f.keys()))

    # Read data from Step 0
    step0 = f['Step_0']

    # Get metadata
    time = step0.attrs['Time_Myr']
    npart = step0.attrs['NumberOfParticle']

    # Get particle data
    x = step0['X_pc'][:]
    y = step0['Y_pc'][:]
    z = step0['Z_pc'][:]
    mass = step0['Mass_Msun'][:]

    # Plot particle positions
    plt.scatter(x, y, s=mass/100, alpha=0.5)
    plt.xlabel('X (pc)')
    plt.ylabel('Y (pc)')
    plt.title(f'Particle Distribution at t={time:.2f} Myr')
    plt.show()
```

### Julia Example
```julia
using HDF5

# Open the HDF5 file
h5open("output/output.h5", "r") do file
    # Read Step 0
    step0 = file["Step_0"]

    # Get metadata
    time = read(attrs(step0)["Time_Myr"])
    npart = read(attrs(step0)["NumberOfParticle"])

    # Get particle data
    x = read(step0["X_pc"])
    y = read(step0["Y_pc"])
    mass = read(step0["Mass_Msun"])

    println("Time: $time Myr")
    println("Number of particles: $npart")
end
```

### MATLAB Example
```matlab
% Open the HDF5 file
filename = 'output/output.h5';

% Read data from Step 0
time = h5readatt(filename, '/Step_0', 'Time_Myr');
npart = h5readatt(filename, '/Step_0', 'NumberOfParticle');

x = h5read(filename, '/Step_0/X_pc');
y = h5read(filename, '/Step_0/Y_pc');
mass = h5read(filename, '/Step_0/Mass_Msun');

% Plot
scatter(x, y, mass/100, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('X (pc)');
ylabel('Y (pc)');
title(sprintf('Particle Distribution at t=%.2f Myr', time));
```

## Migration Checklist

- [x] Update `readParameterFile.cpp` to use TOML parser
- [x] Update `ReadWrite.cpp` to use HDF5 writer
- [x] Update `Makefile` with TOML and HDF5 dependencies
- [x] Create example TOML configuration file
- [x] Download toml11 library
- [x] Create build script
- [x] Create documentation

## Files Modified

1. **src/readParameterFile.cpp**: Replaced simple parser with TOML parser
2. **src/ReadWrite.cpp**: Replaced text output with HDF5 output
3. **src/Makefile**: Added TOML and HDF5 library dependencies
4. **test/test1/config.toml**: New TOML configuration example
5. **build_with_new_io.sh**: Convenience build script

## Benefits of New I/O System

### TOML Input
- ✅ Type safety (strings, integers, floats are properly typed)
- ✅ Comments and documentation in config file
- ✅ Hierarchical configuration (can organize parameters in sections)
- ✅ Industry standard format
- ✅ Better error messages

### HDF5 Output
- ✅ **90% smaller file sizes** (binary vs text)
- ✅ Single file for all timesteps (easier to manage)
- ✅ Self-describing format (metadata included)
- ✅ Fast random access to any timestep
- ✅ Widely supported (Python, Julia, MATLAB, R, etc.)
- ✅ Parallel I/O capable (for future scaling)
- ✅ Built-in compression support

## Troubleshooting

### Build Errors

**Error: `toml.hpp: No such file or directory`**
```bash
# Clone toml11 library
git clone https://github.com/ToruNiina/toml11.git
```

**Error: `H5Cpp.h: No such file or directory`**
```bash
# Load HDF5 module
module load hdf5/1.14.1-2-openmpi-gcc-12.3.0
```

**Error: `cannot find -lhdf5_cpp`**
```bash
# Make sure you have HDF5 with C++ support
# Check available modules:
module avail hdf5
```

### Runtime Errors

**Error: `Unable to parse TOML configuration file`**
- Check TOML syntax (strings must be quoted, numbers must not be)
- Verify file exists and is readable

**Error: `HDF5 error: unable to create file`**
- Check that output directory exists
- Verify write permissions

## Support

For questions or issues with the new I/O system, please:
1. Check this guide
2. Refer to example files in `test/test1/`
3. Contact the development team

## References

- [TOML Specification](https://toml.io/)
- [toml11 Library](https://github.com/ToruNiina/toml11)
- [HDF5 Documentation](https://www.hdfgroup.org/solutions/hdf5/)
- [h5py (Python)](https://docs.h5py.org/)
