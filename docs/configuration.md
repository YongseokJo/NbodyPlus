# Configuration Reference

ABYSS uses TOML configuration files. This document describes all available parameters.

## Minimal Configuration

```toml
Filename = "nbody.dat"
StopTime = 1e7              # End time in years
OutputDirectory = "output"
```

## Complete Configuration Example

```toml
# === Required Parameters ===
Filename = "nbody.dat"           # Initial conditions file
StopTime = 1e7                   # End time in years (10 Myr)
OutputDirectory = "output"       # Output directory

# === Numerics ===
[numerics]
eta = 0.01                       # Timestep parameter (smaller = more accurate)
fixed_num_neighbors = 100        # Number of neighbors for force calculation
InitialRadius = 0.2              # Initial neighbor search radius [pc]
r_search = 2.5e-4                # Few-body search radius [pc]
t_search = 1e-6                  # Few-body search time [Myr]

# === Output ===
[output]
dtOutput = 1e6                   # Output interval [years]
Compression = true               # Enable HDF5 compression
compression_level = 6            # Compression level (1-9)

# === Restart ===
[restart]
Enabled = false                  # Enable restart from checkpoint
checkpoint_file = ""             # Checkpoint file path

# === McLuster IC Generation (optional) ===
[mcluster]
N = 10000                        # Number of stars
P = 0                            # Density profile (0=Plummer)
R = 0.8                          # Half-mass radius [pc]
f = 1                            # IMF (1=Kroupa)
Z = 0.02                         # Metallicity
b = 0.0                          # Binary fraction
e = 0.0                          # Stellar evolution epoch [Myr]
generate_only = false            # If true, exit after IC generation
```

---

## Required Parameters

### `Filename`
**Type:** String
**Required:** Yes

Path to initial conditions file. If `[mcluster]` section is present, this file will be created/overwritten by McLuster.

```toml
Filename = "nbody.dat"
Filename = "/path/to/initial_conditions.dat"
```

### `StopTime`
**Type:** Float
**Required:** Yes
**Units:** Years

Simulation end time.

```toml
StopTime = 1e7      # 10 Myr
StopTime = 1e8      # 100 Myr
StopTime = 1e9      # 1 Gyr
```

### `OutputDirectory`
**Type:** String
**Default:** `"output"`

Directory for output files. Created if it doesn't exist.

```toml
OutputDirectory = "output"
OutputDirectory = "/scratch/user/run001"
```

---

## Numerics Section

### `eta`
**Type:** Float
**Default:** `0.01`
**Range:** > 0

Timestep accuracy parameter. Smaller values give more accurate but slower simulations.

```toml
eta = 0.01     # Default, good balance
eta = 0.001    # Higher accuracy
eta = 0.1      # Faster but less accurate
```

### `fixed_num_neighbors`
**Type:** Integer
**Default:** `100`
**Range:** 10 to MAX_NUM_NEIGHBOR

Number of neighbors used in force calculations.

```toml
fixed_num_neighbors = 100    # Default
fixed_num_neighbors = 200    # More neighbors, higher accuracy
```

### `InitialRadius`
**Type:** Float
**Default:** `0.2`
**Units:** Parsecs

Initial neighbor search radius.

```toml
InitialRadius = 0.2    # Default
InitialRadius = 0.5    # Larger initial search
```

### `r_search`
**Type:** Float
**Default:** `2.5e-4`
**Units:** Parsecs

Search radius for few-body interactions.

```toml
r_search = 2.5e-4    # Default
r_search = 1e-3      # Larger search radius
```

### `t_search`
**Type:** Float
**Default:** `1e-6`
**Units:** Myr

Time criterion for few-body search.

```toml
t_search = 1e-6    # Default
```

---

## Output Section

### `dtOutput`
**Type:** Float
**Default:** `StopTime / 10`
**Units:** Years

Output interval for snapshots.

```toml
[output]
dtOutput = 1e6      # Output every 1 Myr
dtOutput = 1e5      # Output every 0.1 Myr (more frequent)
```

### `Compression`
**Type:** Boolean
**Default:** `true`

Enable HDF5 compression for output files.

```toml
[output]
Compression = true     # Compressed output (smaller files)
Compression = false    # Uncompressed (faster I/O)
```

### `compression_level`
**Type:** Integer
**Default:** `6`
**Range:** 1-9

HDF5 compression level. Higher = smaller files, slower writes.

```toml
[output]
compression_level = 1    # Fast, less compression
compression_level = 6    # Default, good balance
compression_level = 9    # Maximum compression, slowest
```

---

## Restart Section

### `Enabled`
**Type:** Boolean
**Default:** `false`

Enable checkpoint restart.

```toml
[restart]
Enabled = true
checkpoint_file = "checkpoint.h5"
```

### `checkpoint_file`
**Type:** String
**Default:** `""`

Path to checkpoint file for restart.

```toml
[restart]
Enabled = true
checkpoint_file = "output/checkpoint_001.h5"
```

---

## McLuster Section

The `[mcluster]` section enables automatic initial condition generation using the McLuster cluster generator. If this section is present, ABYSS will:

1. Run McLuster to generate initial conditions
2. Transform output to ABYSS format (unit conversion)
3. Continue to simulation (or exit if `generate_only = true`)

### `N`
**Type:** Integer
**Default:** `0`

Number of stars. Must specify either `N` or `M` (not both).

```toml
[mcluster]
N = 1000        # 1,000 stars
N = 100000      # 100,000 stars
N = 1000000     # 1 million stars
```

**Minimum:** 3 (N-body simulation requires at least 3 particles)

### `M`
**Type:** Float
**Default:** `0.0`
**Units:** Solar masses (Msun)

Total cluster mass. Alternative to `N`. If both specified, `M` takes precedence.

```toml
[mcluster]
M = 5000.0      # 5,000 solar mass cluster
```

### `P`
**Type:** Integer
**Default:** `0`
**Range:** -1 to 3

Density profile type.

| Value | Profile | Description |
|-------|---------|-------------|
| -1 | None | No density gradient |
| 0 | Plummer | Plummer sphere (default) |
| 1 | King | King model |
| 2 | Subr | Subr et al. model |
| 3 | EFF/Nuker | EFF or Nuker profile |

```toml
[mcluster]
P = 0    # Plummer (most common)
P = 1    # King model
```

### `R`
**Type:** Float
**Default:** `0.8`
**Units:** Parsecs

Half-mass radius of the cluster.

```toml
[mcluster]
R = 0.8     # Default, typical open cluster
R = 2.0     # Larger cluster
R = 0.3     # Compact cluster
```

### `f`
**Type:** Integer
**Default:** `1`
**Range:** 0-2

Initial Mass Function (IMF) selection.

| Value | IMF | Description |
|-------|-----|-------------|
| 0 | Single mass | All stars have same mass |
| 1 | Kroupa | Kroupa (2001) IMF (default) |
| 2 | User-defined | Custom IMF from file |

```toml
[mcluster]
f = 1    # Kroupa IMF (recommended)
f = 0    # Equal mass (for testing)
```

### `Z`
**Type:** Float
**Default:** `0.02`
**Range:** 0.0001 to 0.03

Metallicity (mass fraction of metals).

| Value | Description |
|-------|-------------|
| 0.02 | Solar metallicity (default) |
| 0.001 | Low metallicity (halo stars) |
| 0.03 | High metallicity |

```toml
[mcluster]
Z = 0.02      # Solar (default)
Z = 0.001     # Metal-poor
Z = 0.03      # Metal-rich
```

### `b`
**Type:** Float
**Default:** `0.0`
**Range:** 0.0 to 1.0

Binary fraction (0 = no binaries, 1 = all binaries).

```toml
[mcluster]
b = 0.0     # No binaries (default)
b = 0.5     # 50% binaries
b = 1.0     # All stars in binaries
```

### `e`
**Type:** Float
**Default:** `0.0`
**Units:** Myr

Stellar evolution epoch. Evolve stars before simulation starts.

```toml
[mcluster]
e = 0.0     # No pre-evolution (default)
e = 10.0    # Pre-evolve for 10 Myr
e = 100.0   # Pre-evolve for 100 Myr
```

### `generate_only`
**Type:** Boolean
**Default:** `false`

If true, exit after generating initial conditions (no simulation).

```toml
[mcluster]
N = 100000
generate_only = true    # Generate IC only, then exit
```

---

## Example Configurations

### Simple Simulation with External IC

```toml
Filename = "cluster.dat"
StopTime = 5e7
OutputDirectory = "run_001"
```

### Large Cluster with McLuster

```toml
Filename = "nbody.dat"
StopTime = 1e8
OutputDirectory = "large_cluster"

[numerics]
eta = 0.005
fixed_num_neighbors = 150

[output]
dtOutput = 5e5
Compression = true
compression_level = 6

[mcluster]
N = 100000
P = 0
R = 1.0
f = 1
Z = 0.02
b = 0.3
```

### Generate IC Only

```toml
Filename = "plummer_1M.dat"
StopTime = 1e6
OutputDirectory = "ic_output"

[mcluster]
N = 1000000
P = 0
R = 0.8
generate_only = true
```

### Metal-Poor Globular Cluster

```toml
Filename = "nbody.dat"
StopTime = 1e10
OutputDirectory = "globular"

[mcluster]
N = 50000
P = 1         # King model
R = 3.0
Z = 0.001     # Low metallicity
b = 0.1
e = 1000.0    # Pre-evolve 1 Gyr
```

---

## Configuration Summary Printout

When ABYSS runs, it prints a configuration summary:

```
========== ABYSS Configuration ==========
Input file:        nbody.dat
Output directory:  output

--- Time ---
End time:          10 Myr
Output interval:   1 Myr

--- Numerics ---
eta:               0.01
fixed_num_neighbors:    100
InitialRadius:     0.2 pc
r_search:          0.00025 pc
t_search:          1e-06 Myr

--- Output ---
Compression:       enabled
Compression level: 6

--- Restart ---
Restart enabled:   no

--- McLuster IC Generation ---
N (star count):    10000
P (profile):       0 (Plummer)
R (half-mass):     0.8 pc
f (IMF):           1 (Kroupa)
Z (metallicity):   0.02
b (binary frac):   0
e (epoch):         0 Myr
generate_only:     no
==========================================
```
