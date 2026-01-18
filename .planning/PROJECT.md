# ABYSS SoA Conversion

## What This Is

ABYSS N-body simulation code with Structure of Arrays (SoA) data layout. The v1.0 conversion replaced the original Array of Structures (AoS) `Particle` struct with a `ParticleData` SoA container (66 arrays), `ParticleDataMPI` for shared memory (66 MPI_Win handles), and `ParticleDataGPU` for device data. SDAR few-body integration maintained via exit-only sync pattern.

## Core Value

**Physics correctness (energy conservation) with clean SoA architecture for future optimizations.**

## Requirements

### Validated

<!-- Shipped in v1.0 -->

- ✓ ParticleData SoA container with 66 arrays — v1.0
- ✓ Accessor functions (get/set/pointer) for all fields — v1.0
- ✓ ParticleDataMPI with 66 MPI_Win shared memory windows — v1.0
- ✓ ParticleDataGPU device container with async transfers — v1.0
- ✓ SoA-compatible CUDA force kernels — v1.0
- ✓ CPU routines with SoA helper functions — v1.0
- ✓ SDAR compatibility via exit-only sync pattern — v1.0
- ✓ I/O sync (initialization and checkpoint restore) — v1.0
- ✓ Energy conservation validated (2.52e-5 vs 3.21e-5 baseline) — v1.0

<!-- Existing capabilities maintained -->

- ✓ 4th-order Hermite integration with block timesteps — existing
- ✓ GPU-accelerated force calculations (CUDA) — existing
- ✓ MPI root/worker parallelism with shared memory — existing
- ✓ Few-body dynamics via SDAR (binaries, multiples) — existing
- ✓ Post-Newtonian corrections (PN1.0, PN2.0, PN2.5) — existing
- ✓ Merger handling (GW-driven, TDE, stellar collisions) — existing
- ✓ HDF5 output with compression — existing
- ✓ TOML configuration parsing — existing
- ✓ Optional SEVN stellar evolution integration — existing

### Active

<!-- v2.0 Performance Profiling & Optimization -->

- [ ] Enhanced profiler with per-rank statistics (min/avg/max across MPI ranks)
- [ ] Sub-timers for irregular force breakdown (NeighborLoop, CMLoop, Correction)
- [ ] Work counters (particles processed, neighbor pairs evaluated)
- [ ] Histogram support for call-time distributions
- [ ] Load balance metrics and analysis output
- [ ] Worker-side timing (MPI_Recv wait, task dispatch, send completion)
- [ ] Queue scheduler timing (assign, run operations)
- [ ] Profile runs to identify actual bottleneck
- [ ] Optimization based on profiling findings (vectorization, MPI, load balancing)

### Out of Scope

- SDAR `Group` struct conversion — kept as-is with exit-only sync
- Performance improvement — v1.0 showed no change; bottleneck elsewhere
- Algorithm changes — SoA is data layout only

## Context

**Current State (v1.0 shipped):**

- `ParticleData` SoA container in `src/particle_data.h` — 66 arrays with accessors
- `ParticleDataMPI` in `src/particle_data_mpi.h` — 66 MPI_Win handles
- `ParticleDataGPU` in `src/particle_data_gpu.h` — device container
- Original `Particle` struct retained for SDAR compatibility (exit-only sync)
- 57 commits on AoS_to_SoA branch

**Validation Results:**

| Metric | Baseline | v1.0 SoA | Result |
|--------|----------|----------|--------|
| dE/E0 mean | 3.21e-05 | 2.52e-05 | ✓ Better |
| Wall time | 28.33s | 28.40s | No change |

**Key Files:**

- `src/particle_data.h/.cpp` — SoA container
- `src/particle_data_mpi.h/.cpp` — MPI shared memory
- `src/particle_data_gpu.h/.cu` — GPU container
- `src/particle.h` — original AoS struct (kept for SDAR)

**Validation Workflow:**

- Run tests: `workflow/bin/submit.sh --tag <name> --scheduler slurm`
- Review results: `summary_runs.tsv`

## Constraints

- **Tech stack**: C++11, CUDA 12.x, MPI (OpenMPI), HDF5
- **Compatibility**: Must maintain SDAR library integration (uses macro-based field aliases)
- **Correctness**: Energy conservation must remain within existing tolerances
- **Build**: Must work with existing Makefile and workflow scripts

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Full SoA (not hybrid) | Cleaner architecture, no data format conversions at runtime | ✓ Good — simpler code |
| Accessor functions | Encapsulation, easier future changes, explicit intent | ✓ Good — clean API |
| SoA for MPI shared memory | Avoids AoS↔SoA conversion overhead on every access | ✓ Good — 66 windows work fine |
| Keep SDAR Group as-is | Complex library integration, conversion risk outweighs benefit | ✓ Good — exit-only sync works |
| Exit-only sync pattern | SDAR uses AoS, sync to SoA only on exit | ✓ Good — simpler than proxy |

## Current Milestone: v2.0 Performance Profiling & Optimization

**Goal:** Build detailed profiler infrastructure to identify irregular force bottlenecks, then optimize based on findings.

**Target features:**
- Enhanced profiler with per-rank MPI statistics and sub-timers
- Work counters and histogram support for performance analysis
- Bottleneck identification through systematic profiling
- Targeted optimizations (vectorization, MPI improvements, load balancing)

---
*Last updated: 2026-01-17 after v2.0 milestone initialization*
