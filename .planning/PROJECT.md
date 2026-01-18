# ABYSS N-body Simulation

## What This Is

ABYSS N-body simulation code with Structure of Arrays (SoA) data layout and enhanced profiling infrastructure. The v1.0 conversion replaced the original Array of Structures (AoS) with `ParticleData` SoA containers, `ParticleDataMPI` for shared memory, and `ParticleDataGPU` for device data. v2.0 added comprehensive profiler infrastructure with per-rank statistics, sub-timers, work counters, and AVX-512 vectorized irregular force calculation.

## Core Value

**Physics correctness (energy conservation) with clean architecture for targeted optimizations.**

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

<!-- Shipped in v2.0 -->

- ✓ Enhanced profiler with per-rank MPI statistics (min/avg/max) — v2.0
- ✓ Sub-timers for irregular force breakdown (NeighborLoop, CMLoop, Correction) — v2.0
- ✓ Work counters (particles processed, neighbor pairs evaluated) — v2.0
- ✓ Histogram support for call-time distributions — v2.0
- ✓ Load balance metrics and analysis output — v2.0
- ✓ Worker-side timing (MPI_Recv wait, task dispatch, send completion) — v2.0
- ✓ Queue scheduler timing (assign, run operations) — v2.0
- ✓ Bottleneck identified: IrregularForce at 53.5% of wall time — v2.0
- ✓ AVX-512 vectorization of irregular force neighbor loop — v2.0
- ✓ Performance improvement: 3.2% (69.5s → 67.3s IrregularForce) — v2.0
- ✓ Energy conservation maintained (dE/E0 = 1.34e-6) — v2.0

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

<!-- v2.1 MPI Communication Optimization -->

**Async MPI Implementation:**
- [ ] Convert root-side MPI to async operations (MPI_Isend/MPI_Irecv)
- [ ] Replace blocking waits with MPI_Waitany/MPI_Testany
- [ ] Overlap communication with computation on root
- [ ] Evaluate worker-side async (if beneficial)
- [ ] Measure MPI wait time reduction
- [ ] Validate energy conservation

**Batching Research:**
- [ ] Research batching strategies for CM particle dependencies
- [ ] Research load balancing for variable work per particle
- [ ] Research callback mechanism changes for batch completion
- [ ] Document recommended approach for v2.2+

### Out of Scope

- SDAR `Group` struct conversion — kept as-is with exit-only sync
- Algorithm changes — infrastructure optimization only
- GPU irregular forces — deferred to v2.2+
- Implementing batching — research only for this milestone, implementation deferred

## Context

**Current State (v2.0 shipped):**

- `src/profiler.h` — Enhanced profiler with MPI aggregation, histograms, work counters
- `src/simd_force.h/.cpp` — AVX-512 vectorized force kernel with pre-gather pattern
- `src/Particle/compute_acceleration.cpp` — Integrated vectorized neighbor loop
- IrregularForce bottleneck identified and addressed (53.5% → still dominant but 3.2% faster)
- Energy conservation: ✓ Pass (dE/E0 = 1.34e-6)

**Optimization Results:**

| Metric | v2.0 Baseline | v2.0 Optimized | Result |
|--------|---------------|----------------|--------|
| IrregularForce | 69.5s | 67.3s | -3.2% |
| Wall time | 130s | 126.4s | -2.8% |
| dE/E0 | — | 1.34e-6 | ✓ Pass |

**Key Files:**

- `src/particle_data.h/.cpp` — SoA container (v1.0)
- `src/particle_data_mpi.h/.cpp` — MPI shared memory (v1.0)
- `src/particle_data_gpu.h/.cu` — GPU container (v1.0)
- `src/profiler.h` — Enhanced profiler (v2.0)
- `src/simd_force.h/.cpp` — AVX-512 vectorization (v2.0)

**Validation Workflow:**

- Run tests: `workflow/bin/submit.sh --tag <name> --scheduler slurm`
- Review results: `summary_runs.tsv`

## Constraints

- **Tech stack**: C++11, CUDA 12.x, MPI (OpenMPI), HDF5
- **Compatibility**: Must maintain SDAR library integration (uses macro-based field aliases)
- **Correctness**: Energy conservation must remain within existing tolerances
- **Build**: Must work with existing Makefile and workflow scripts
- **AVX-512**: Requires Skylake-AVX512 or later CPU architecture

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Full SoA (not hybrid) | Cleaner architecture, no data format conversions at runtime | ✓ Good — simpler code |
| Accessor functions | Encapsulation, easier future changes, explicit intent | ✓ Good — clean API |
| SoA for MPI shared memory | Avoids AoS↔SoA conversion overhead on every access | ✓ Good — 66 windows work fine |
| Keep SDAR Group as-is | Complex library integration, conversion risk outweighs benefit | ✓ Good — exit-only sync works |
| Exit-only sync pattern | SDAR uses AoS, sync to SoA only on exit | ✓ Good — simpler than proxy |
| Pre-gather for SIMD | Align scattered neighbor data before vectorized compute | ⚠ Revisit — gather overhead limits gains |
| AVX-512 vectorization | Target IrregularForce bottleneck with SIMD | ⚠ Limited — 3.2% gain vs 20% target |

## Current Milestone: v2.1 MPI Communication Optimization

**Goal:** Reduce MPI overhead in irregular force communication through async operations, and research batching strategies for future work.

**Key insight from v2.0 profiling:**
- ~105M MPI messages per interval for ~140K particle evaluations
- Current pattern: 1 particle per message, blocking send/recv
- Root and workers spend significant time waiting on MPI operations

**Approach:**
1. Convert to async MPI (MPI_Isend/MPI_Irecv) starting with root-side
2. Replace blocking waits with MPI_Waitany/MPI_Testany
3. Research batching strategies (complex due to CM dependencies, variable work, callback mechanism)

## Recommendations for v2.2+

1. **Implement batching** — Based on v2.1 research findings
2. **SIMD gather intrinsics** — Use `_mm512_i64gather_pd` to avoid pre-gather copies
3. **GPU irregular forces** — Port irregular force kernel to GPU

---
*Last updated: 2026-01-18 after v2.1 milestone start*
