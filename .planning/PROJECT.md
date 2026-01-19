# ABYSS N-body Simulation

## What This Is

ABYSS N-body simulation code with Structure of Arrays (SoA) data layout and enhanced profiling infrastructure. The v1.0 conversion replaced the original Array of Structures (AoS) with `ParticleData` SoA containers, `ParticleDataMPI` for shared memory, and `ParticleDataGPU` for device data. v2.0 added comprehensive profiler infrastructure with per-rank statistics, sub-timers, work counters, and AVX-512 vectorized irregular force calculation.

## Core Value

**Physics correctness (energy conservation) with clean architecture for targeted optimizations.**

## Current Milestone: v2.2 Load Balance Profiling

**Goal:** Instrument and analyze load balance across MPI workers to understand where imbalance comes from before committing to optimization approach (OpenMP, batching, work stealing, etc.).

**Target measurements:**
- Neighbor count variance (per-particle counts, distribution, correlation with compute time)
- Queue dispatch overhead (dispatch latency, root overhead, queue depth over time)
- Worker distribution (particles per worker, compute time per worker, heavy particle detection)
- CM particle breakdown (time comparison: CM vs regular particles)
- Few-body overhead (search/init/integration time breakdown)
- Memory access patterns (cache miss rates in force loops)

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

<!-- v2.2 Load Balance Profiling -->

- [ ] Neighbor count profiling (per-particle, distribution, outlier detection)
- [ ] Queue dispatch latency measurement (time between completion and next task)
- [ ] Queue depth tracking over time (worker starvation detection)
- [ ] Per-worker particle distribution analysis
- [ ] Per-worker compute time breakdown
- [ ] CM particle timing comparison (CM vs regular particles)
- [ ] Few-body timing breakdown (search, init, integration separately)
- [ ] Cache miss rate measurement in force loops
- [ ] Load balance analysis report with recommendations

### Out of Scope

- SDAR `Group` struct conversion — kept as-is with exit-only sync
- Algorithm changes — profiling work is measurement only
- Optimization implementation — v2.2 is analysis; optimization deferred to v2.3+
- Async MPI — archived in v2.1, single-message async doesn't help

## Context

**Current State (v2.0 shipped, v2.1 archived):**

- `src/profiler.h` — Enhanced profiler with MPI aggregation, histograms, work counters
- `src/simd_force.h/.cpp` — AVX-512 vectorized force kernel with pre-gather pattern
- `src/Particle/compute_acceleration.cpp` — Integrated vectorized neighbor loop
- IrregularForce bottleneck identified (53.5% of wall time)
- v2.1 async MPI archived — single-message async adds overhead, doesn't help

**v2.1 Archive Learnings:**

- Async MPI at single-message granularity is slower than blocking (126.4B vs 125.0B ns)
- ~105M MPI messages/interval creates unavoidable overhead at current granularity
- Future optimization must reduce message count (batching) or change parallelization model

**Key Files:**

- `src/particle_data.h/.cpp` — SoA container (v1.0)
- `src/particle_data_mpi.h/.cpp` — MPI shared memory (v1.0)
- `src/particle_data_gpu.h/.cu` — GPU container (v1.0)
- `src/profiler.h` — Enhanced profiler (v2.0)
- `src/queue_scheduler.h` — Worker dispatch (target for load balance profiling)
- `src/irregular_routines.cpp` — Main irregular force loop (target for profiling)

**Validation Workflow:**

- Run tests: `workflow/bin/submit.sh --tag <name> --scheduler slurm`
- Review results: `summary_runs.tsv`

## Constraints

- **Tech stack**: C++11, CUDA 12.x, MPI (OpenMPI), HDF5
- **Compatibility**: Must maintain SDAR library integration (uses macro-based field aliases)
- **Correctness**: Energy conservation must remain within existing tolerances
- **Build**: Must work with existing Makefile and workflow scripts
- **AVX-512**: Requires Skylake-AVX512 or later CPU architecture
- **Profiling overhead**: New instrumentation should not significantly impact runtime (< 5% overhead)

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
| Archive async MPI (v2.1) | Single-message async adds overhead, doesn't help | ✓ Good — data-driven decision |
| Profile before optimizing | Understand load imbalance source before committing to fix | — Pending |

## Recommendations for v2.3+

After v2.2 analysis, likely directions:
1. **MPI batching** — Reduce message count if queue dispatch is bottleneck
2. **OpenMP hybrid** — Add thread parallelism if workers are underutilized
3. **Work stealing** — Rebalance if particle distribution is uneven
4. **GPU irregular forces** — Port if compute is the bottleneck, not communication

---
*Last updated: 2026-01-18 after v2.2 milestone start*
