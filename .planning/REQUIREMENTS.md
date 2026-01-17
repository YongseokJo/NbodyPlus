# ABYSS SoA Conversion — Requirements

## Core Container

| REQ-ID | Requirement | Acceptance Criteria | Status |
|--------|-------------|---------------------|--------|
| CONT-01 | Create `ParticleData` SoA container class | Class exists with separate arrays for all Particle fields | Complete ✓ |
| CONT-02 | Implement accessor functions for all fields | `get_X(i)`, `set_X(i, val)`, and `X()` (pointer) for each field | Complete ✓ |
| CONT-03 | Implement memory allocation/deallocation | `allocate(N)`, `deallocate()` methods work correctly | Complete ✓ |

## MPI Integration

| REQ-ID | Requirement | Acceptance Criteria | Status |
|--------|-------------|---------------------|--------|
| MPI-01 | Allocate each SoA array as separate `MPI_Win` | Multiple `MPI_Win_allocate_shared` calls succeed | Complete ✓ |
| MPI-02 | Workers can access all particle data arrays | Read access works across all MPI ranks | Complete ✓ |
| MPI-03 | Root can update particle data arrays | Write access works from root rank | Complete ✓ |

## GPU Integration

| REQ-ID | Requirement | Acceptance Criteria | Status |
|--------|-------------|---------------------|--------|
| GPU-01 | Create `ParticleDataGPU` device container | Struct with device pointers for all GPU-needed fields | Complete ✓ |
| GPU-02 | Implement host-to-device transfer | `copy_from_host()` transfers all arrays | Complete ✓ |
| GPU-03 | Implement device-to-host transfer | `copy_to_host()` transfers results back | Complete ✓ |
| GPU-04 | Update `compute_forces_kernel` for SoA | Kernel accepts separate array pointers | Complete ✓ |
| GPU-05 | Update `predict_particle_kernel` for SoA | Kernel accepts separate array pointers | N/A (no kernel found) |

## CPU Routines

| REQ-ID | Requirement | Acceptance Criteria |
|--------|-------------|---------------------|
| CPU-01 | Update regular force calculation for SoA | `calculate_regular_acceleration()` uses accessors |
| CPU-02 | Update irregular force calculation for SoA | `calculate_irregular_acceleration()` uses accessors |
| CPU-03 | Update prediction routines for SoA | All prediction functions use accessors |
| CPU-04 | Update correction routines for SoA | All correction functions use accessors |
| CPU-05 | Update timestep routines for SoA | Timestep calculations use accessors |

## SDAR Compatibility

| REQ-ID | Requirement | Acceptance Criteria |
|--------|-------------|---------------------|
| SDAR-01 | Create `ParticleProxy` adapter class | Proxy provides SDAR-compatible interface |
| SDAR-02 | Implement `sync_from_soa()` | Data copies from SoA to proxy buffers |
| SDAR-03 | Implement `sync_to_soa()` | Modified data copies back to SoA |

## I/O

| REQ-ID | Requirement | Acceptance Criteria |
|--------|-------------|---------------------|
| IO-01 | Update HDF5 output for SoA | `output.h5` contains correct particle data |
| IO-02 | Update particle initialization from file | Reading from input files populates SoA correctly |

## Validation

| REQ-ID | Requirement | Acceptance Criteria |
|--------|-------------|---------------------|
| VAL-01 | Baseline captured before conversion | `summary_runs.tsv` contains baseline row with `dE_over_E0_mean`, `total_wall_s` |
| VAL-02 | Energy conservation matches baseline | Run `workflow/bin/submit.sh --tag soa_test`, verify `dE_over_E0_mean` in `summary_runs.tsv` is ≤ baseline (≈3.2e-5) |
| VAL-03 | All tests pass | `workflow/bin/submit.sh` completes with `status=ok (Simulation Done found)` |
| VAL-04 | Performance improvement measurable | Compare `total_wall_s` in `summary_runs.tsv` between baseline and SoA runs |

---

## Summary

| Category | Count |
|----------|-------|
| Core Container | 3 |
| MPI Integration | 3 |
| GPU Integration | 5 |
| CPU Routines | 5 |
| SDAR Compatibility | 3 |
| I/O | 2 |
| Validation | 4 |
| **Total** | **25** |

---
*Generated: 2026-01-17*
