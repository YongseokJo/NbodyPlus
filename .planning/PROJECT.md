# ABYSS SoA Conversion

## What This Is

A performance optimization project to convert ABYSS's core data structures from Array of Structures (AoS) to Structure of Arrays (SoA) layout. This affects the `Particle` struct, GPU data types (`j_particle_t`, `i_particle_t`), and MPI shared memory windows. The goal is improved CPU cache efficiency and GPU memory coalescing while maintaining simulation correctness.

## Core Value

**Improve simulation performance through optimized memory layout while maintaining physics correctness (energy conservation, test reproducibility).**

## Requirements

### Validated

<!-- Existing capabilities that must continue working -->

- ✓ 4th-order Hermite integration with block timesteps — existing
- ✓ GPU-accelerated force calculations (CUDA) — existing
- ✓ MPI root/worker parallelism with shared memory — existing
- ✓ Few-body dynamics via SDAR (binaries, multiples) — existing
- ✓ Post-Newtonian corrections (PN1.0, PN2.0, PN2.5) — existing
- ✓ Merger handling (GW-driven, TDE, stellar collisions) — existing
- ✓ HDF5 output with compression — existing
- ✓ TOML configuration parsing — existing
- ✓ Optional SEVN stellar evolution integration — existing
- ✓ Energy conservation within tolerance — existing

### Active

<!-- Current scope being built toward -->

- [ ] Convert `Particle` struct to SoA container with accessor functions
- [ ] Convert GPU types (`j_particle_t`, `i_particle_t`) to SoA
- [ ] Convert MPI shared memory to multiple `MPI_Win` windows (one per array)
- [ ] Update all particle data access to use accessor functions
- [ ] Update CUDA kernels to work with SoA layout
- [ ] Update CPU force calculations to work with SoA layout
- [ ] Maintain SDAR/Group compatibility with new data layout
- [ ] Validate energy conservation matches pre-conversion baseline
- [ ] Benchmark performance improvement vs AoS baseline

### Out of Scope

- SDAR `Group` struct conversion — complex integration with SDAR library, defer
- Hybrid AoS/SoA approach — decided on full conversion for cleaner architecture
- Algorithm changes — this is purely a data layout optimization
- New physics features — focus on conversion only

## Context

**Current Architecture:**
- `Particle` struct in `src/particle.h` contains ~120 fields including position[3], velocity[3], acceleration arrays [3][4], neighbor info, timestep data, and SDAR/SEVN pointers
- GPU types defined in `src/cuda/cuda_defs.h`
- MPI uses `MPI_Win_allocate_shared` for particle array in `src/mpi_routines.cpp`
- Force calculations in `src/cuda/cuda_kernels.cu` and `src/Particle/compute_acceleration.cpp`

**Key Files to Modify:**
- `src/particle.h` — main struct → SoA container
- `src/cuda/cuda_defs.h` — GPU types
- `src/mpi_routines.cpp` — shared memory allocation
- `src/cuda/cuda_kernels.cu` — kernel access patterns
- `src/cuda/cuda_acceleration.cu` — GPU dispatch
- All files accessing `particles[i].field`

**Existing Tests:**
- `tests/test1/` with TOML config
- `tools/analyze_energy.py` for energy conservation validation
- Performance profiling via `PERFORMANCETRACE` flag

## Constraints

- **Tech stack**: C++11, CUDA 12.x, MPI (OpenMPI), HDF5
- **Compatibility**: Must maintain SDAR library integration (uses macro-based field aliases)
- **Correctness**: Energy conservation must remain within existing tolerances
- **Build**: Must work with existing Makefile and workflow scripts

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Full SoA (not hybrid) | Cleaner architecture, no data format conversions at runtime | — Pending |
| Accessor functions | Encapsulation, easier future changes, explicit intent | — Pending |
| SoA for MPI shared memory | Avoids AoS↔SoA conversion overhead on every access | — Pending |
| Keep SDAR Group as-is | Complex library integration, conversion risk outweighs benefit | — Pending |

---
*Last updated: 2026-01-17 after initialization*
