# ABYSS N-body Simulation

## What This Is

ABYSS N-body simulation code with Structure of Arrays (SoA) data layout and comprehensive profiling infrastructure. The v1.0 conversion replaced the original Array of Structures (AoS) with `ParticleData` SoA containers, `ParticleDataMPI` for shared memory, and `ParticleDataGPU` for device data. v2.0-v2.3 added profiling infrastructure with MPI aggregation, bottleneck analysis identifying dispatch starvation as the primary optimization target.

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

<!-- Shipped in v2.2 -->

- ✓ Per-particle neighbor count tracking and statistics — v2.2
- ✓ Queue dispatch profiling (latency, depth, starvation detection) — v2.2
- ✓ Worker distribution tracking (particles/worker, compute time) — v2.2
- ✓ CM particle vs regular particle timing breakdown — v2.2
- ✓ Cache miss rate measurement via perf_event — v2.2
- ✓ Comprehensive analysis report with optimization recommendations — v2.2

<!-- Shipped in v2.3 -->

- ✓ Fixed MPI aggregation (worker data flows to root correctly) — v2.3
- ✓ Descriptive cache statistics status messages — v2.3
- ✓ JSON schema versioning and summary section — v2.3
- ✓ Bottleneck analysis: 66% compute, 28% MPI, 1.009 load balance ratio — v2.3
- ✓ Amdahl's Law estimates: 16-27% speedup from MPI batching — v2.3

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

(Ready for next milestone — see `/gsd:new-milestone`)

### Out of Scope

- SDAR `Group` struct conversion — kept as-is with exit-only sync
- Algorithm changes — profiling work is measurement only
- Async MPI — archived in v2.1, single-message async doesn't help
- OpenMP hybrid — load balance excellent (1.009), not needed
- Work stealing — particle distribution already balanced
- GPU irregular forces — compute not the bottleneck (MPI is)

## Context

**Current State (v2.3 shipped):**

- `src/profiler.h` — Enhanced profiler with MPI aggregation, histograms, JSON schema v2.3
- `src/simd_force.h/.cpp` — AVX-512 vectorized force kernel with pre-gather pattern
- `src/Particle/compute_acceleration.cpp` — Integrated vectorized neighbor loop
- `.planning/ANALYSIS.md` — Comprehensive bottleneck analysis report

**v2.3 Analysis Findings:**

- IrregularForce: 66% of wall time (primary compute)
- MPI overhead: 28% of wall time (10.3M messages/interval)
- Load balance ratio: 1.009 (excellent — no optimization needed)
- Dispatch starvation: 2.4M events/interval (primary bottleneck)
- Neighbor distribution: 86.5% have 103-180 neighbors (tight, low variance)

**Optimization Roadmap (deferred to v2.4):**

1. **MPI batching** — Reduce 10.3M messages to 100K-1M (16-27% expected speedup)
2. **Dispatch pipelining** — Hide dispatch latency (2-5% additional)
3. **Verification** — Baseline vs optimized benchmarking

**Key Files:**

- `src/particle_data.h/.cpp` — SoA container (v1.0)
- `src/particle_data_mpi.h/.cpp` — MPI shared memory (v1.0)
- `src/particle_data_gpu.h/.cu` — GPU container (v1.0)
- `src/profiler.h` — Enhanced profiler (v2.0-v2.3)
- `src/queue_scheduler.h` — Worker dispatch (target for v2.4 batching)
- `src/queue.h` — Queue structures (BatchedQueue added for v2.4)

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
| Profile before optimizing | Understand bottleneck source before committing to fix | ✓ Good — identified MPI dispatch |
| MPI batching over alternatives | Incremental change, handles CM particles, preserves dynamic balance | — Pending (v2.4) |
| Rescope v2.3 | Ship analysis, defer optimization to v2.4 | ✓ Good — clean milestone boundary |

---
*Last updated: 2026-01-20 after v2.3 milestone completion*
