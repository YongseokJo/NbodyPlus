# Input/Output Formats

This document describes the file formats used by ABYSS for initial conditions and output.

## Initial Conditions (Input)

### Format: ASCII Text (nbody.dat)

ABYSS reads initial conditions from a whitespace-delimited ASCII file with 7 columns:

```
x y z vx vy vz mass
```

| Column | Quantity | Units |
|--------|----------|-------|
| 1 | x position | kpc |
| 2 | y position | kpc |
| 3 | z position | kpc |
| 4 | x velocity | km/s |
| 5 | y velocity | km/s |
| 6 | z velocity | km/s |
| 7 | mass | 10^-9 Msun |

### Example IC File

```
# nbody.dat (comments not supported, this is for illustration)
0.000123  0.000045 -0.000012  1.234 -0.567  0.890  0.95
0.000456 -0.000078  0.000034  0.123  2.345 -1.234  1.05
-0.000234  0.000089  0.000067 -1.456  0.789  0.456  0.85
...
```

### Unit Conversion

ABYSS uses internal code units. The input units are converted as follows:

| Input | Code Unit Conversion | Constant |
|-------|---------------------|----------|
| Position (kpc) | `pos_code = pos_kpc * 1000 / POSITION_UNIT` | POSITION_UNIT = 4.0 |
| Velocity (km/s) | `vel_code = vel_kms / VELOCITY_UNIT` | VELOCITY_UNIT = 4e-10 |
| Mass (10^-9 Msun) | `mass_code = mass * 1e9 / MASS_UNIT` | MASS_UNIT = 0.0001424198 |

### McLuster-Generated ICs

When using the `[mcluster]` section, ABYSS automatically:

1. Runs McLuster with specified parameters
2. Transforms McLuster output (Msun, pc, km/s) to ABYSS input format (10^-9 Msun, kpc, km/s)
3. Writes the transformed IC file

**McLuster output format:**
```
mass x y z vx vy vz [extras...]
```
Units: Msun, pc, km/s

**Transformed to ABYSS input format:**
```
x y z vx vy vz mass
```
Units: kpc, km/s, 10^-9 Msun

**Conversion applied:**
- Position: pc → kpc (divide by 1000)
- Mass: Msun → 10^-9 Msun units (divide by 10^9)
- Velocity: unchanged (km/s)
- Column order: reordered from mass-first to mass-last

---

## Output Files

ABYSS writes output in HDF5 format.

### Directory Structure

```
output/
├── snapshot_000.h5
├── snapshot_001.h5
├── snapshot_002.h5
├── ...
└── checkpoint_NNN.h5  (if restart enabled)
```

### Snapshot Files

Each snapshot contains the full system state at a given time.

**File naming:** `snapshot_NNN.h5` where NNN is zero-padded sequence number.

**Contents:**

| Dataset | Shape | Type | Description |
|---------|-------|------|-------------|
| `/x` | (N,) | float64 | x positions (code units) |
| `/y` | (N,) | float64 | y positions (code units) |
| `/z` | (N,) | float64 | z positions (code units) |
| `/vx` | (N,) | float64 | x velocities (code units) |
| `/vy` | (N,) | float64 | y velocities (code units) |
| `/vz` | (N,) | float64 | z velocities (code units) |
| `/mass` | (N,) | float64 | masses (code units) |
| `/id` | (N,) | int64 | particle IDs |

**Attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `time` | float64 | Simulation time (code units) |
| `time_myr` | float64 | Simulation time (Myr) |
| `n_particles` | int64 | Number of particles |
| `snapshot_number` | int64 | Snapshot sequence number |

### Reading Snapshots with Python

```python
import h5py
import numpy as np

with h5py.File('output/snapshot_010.h5', 'r') as f:
    # Read particle data
    x = f['x'][:]
    y = f['y'][:]
    z = f['z'][:]
    vx = f['vx'][:]
    vy = f['vy'][:]
    vz = f['vz'][:]
    mass = f['mass'][:]

    # Read metadata
    time_myr = f.attrs['time_myr']
    n_particles = f.attrs['n_particles']

    print(f"Snapshot at t = {time_myr:.2f} Myr")
    print(f"Particles: {n_particles}")
```

### Reading Snapshots with h5dump

```bash
# View structure
h5dump -H output/snapshot_000.h5

# View attributes
h5dump -A output/snapshot_000.h5

# View specific dataset
h5dump -d /x output/snapshot_000.h5 | head -20
```

---

## Checkpoint Files

Checkpoints enable restart of interrupted simulations.

**File naming:** `checkpoint_NNN.h5`

**Contents:** Same as snapshots, plus additional state information for exact restart.

### Restarting from Checkpoint

```toml
[restart]
Enabled = true
checkpoint_file = "output/checkpoint_005.h5"
```

---

## Code Units

ABYSS uses internal code units for computation. Key conversion constants:

| Constant | Value | Meaning |
|----------|-------|---------|
| POSITION_UNIT | 4.0 | 1 code length = 4 pc |
| VELOCITY_UNIT | 4e-10 | Velocity scaling factor |
| MASS_UNIT | 0.0001424198 | Mass scaling factor |
| G | 1.0 | Gravitational constant in code units |

### Converting Output to Physical Units

```python
# Code unit conversion factors (from src/def.h)
POSITION_UNIT = 4.0      # pc
MASS_UNIT = 0.0001424198 # Msun scaling
VELOCITY_UNIT = 4e-10    # km/s scaling

# Convert code units to physical units
def to_physical(x_code, v_code, m_code):
    x_pc = x_code * POSITION_UNIT
    v_kms = v_code * VELOCITY_UNIT
    m_msun = m_code * MASS_UNIT
    return x_pc, v_kms, m_msun
```

---

## Energy and Diagnostics

### Virial Ratio

The virial ratio Q = 2K/|U| indicates dynamical equilibrium:

- Q = 1.0: Virial equilibrium
- Q < 1.0: Contracting (cold)
- Q > 1.0: Expanding (hot)

McLuster generates clusters in approximate virial equilibrium (Q ≈ 0.5 for bound systems, though the ratio 2K/|U| = 1.0).

### Computing Energy from Snapshots

```python
import numpy as np

def compute_energy(x, y, z, vx, vy, vz, mass):
    """Compute kinetic and potential energy from snapshot data."""
    n = len(mass)

    # Kinetic energy: K = 0.5 * sum(m * v^2)
    v2 = vx**2 + vy**2 + vz**2
    K = 0.5 * np.sum(mass * v2)

    # Potential energy: U = -G * sum_i<j(m_i * m_j / r_ij)
    U = 0.0
    for i in range(n):
        for j in range(i+1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dz = z[i] - z[j]
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            if r > 0:
                U -= mass[i] * mass[j] / r

    # Total energy
    E = K + U

    # Virial ratio
    Q = 2.0 * K / abs(U) if U != 0 else 0

    return K, U, E, Q
```

---

## Testing File Formats

### Verify IC File Format

```bash
# Check column count
head -1 nbody.dat | awk '{print NF}'  # Should print 7

# Check line count (number of particles)
wc -l nbody.dat

# Check for invalid values
awk '{for(i=1;i<=NF;i++) if($i !~ /^-?[0-9]/) print NR": invalid"}' nbody.dat
```

### Verify HDF5 Output

```bash
# Check file structure
h5ls -r output/snapshot_000.h5

# Check dataset shapes
h5ls -v output/snapshot_000.h5

# Verify compression
h5dump -p output/snapshot_000.h5 | grep -A5 "FILTERS"
```

---

## File Size Estimates

### Input IC Files

```
N particles × 7 columns × ~20 bytes/value ≈ 140 × N bytes

Examples:
  N = 1,000:      ~140 KB
  N = 10,000:     ~1.4 MB
  N = 100,000:    ~14 MB
  N = 1,000,000:  ~140 MB
```

### Output HDF5 Files (compressed)

Compression ratio depends on data patterns. Typical ratio: 2-5x.

```
Per snapshot: ~50-100 bytes/particle (compressed)

Examples (100 snapshots):
  N = 1,000:      ~5-10 MB total
  N = 10,000:     ~50-100 MB total
  N = 100,000:    ~500 MB - 1 GB total
  N = 1,000,000:  ~5-10 GB total
```
