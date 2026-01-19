# Project Milestones: ABYSS Performance Optimization

## v2.1 MPI Communication Optimization (Archived: 2026-01-18)

**Status:** Archived — async MPI at single-message granularity did not achieve performance improvement

**Phases attempted:** 11-14

**Key findings:**

- Async MPI (MPI_Isend/Irecv/Waitany) adds overhead that exceeds overlap benefits
- Best async configuration: MPI_Recv(ANY_SOURCE) + sync MPI_Send = 126.4B ns
- Baseline blocking: 125.0B ns (1% faster than best async)
- Root cause: ~105M MPI messages/interval creates unavoidable async overhead

**Attempted optimizations:**

1. Pre-posted receives + MPI_Waitany — too much overhead
2. Pre-allocated arrays to avoid vector allocations — insufficient
3. MPI_Recv(ANY_SOURCE) instead of MPI_Waitany — best but still slower
4. MPI_Probe + specific source recv — worse performance
5. Sequential worker waiting — catastrophic (4x slower)

**Conclusion:** Async MPI at single-message granularity does not provide speedup. Future optimization must reduce message count (batching) or change parallelization model.

**Archive location:** Tag `archive-phase13-async-attempts` on branch `MPI_async`

---

## v2.0 Performance Profiling & Optimization (Shipped: 2026-01-18)

**Delivered:** Comprehensive profiling infrastructure to identify bottlenecks, plus targeted AVX-512 vectorization of the irregular force calculation (primary bottleneck at 53.5% of wall time).

**Phases completed:** 8-10 (11 plans total)

**Key accomplishments:**

- Built enhanced profiler with per-rank MPI statistics (min/avg/max across ranks)
- Added sub-timers for irregular force breakdown (NeighborLoop, CMLoop, Correction)
- Implemented work counters, histogram support, and load balance metrics
- Added worker-side timing and queue scheduler profiling
- Identified IrregularForce as primary bottleneck (53.5% of wall time, 122x more calls than regular forces)
- Implemented AVX-512 vectorization of neighbor loop with pre-gather pattern
- Validated optimization (3.2% improvement, energy conservation maintained)

**Stats:**

- 3 phases, 11 plans, 22 requirements
- 19 commits on IrrForce_optimization branch
- 40 files changed, +4,638 / -116 lines
- Performance: 3.2% improvement (69.5s → 67.3s IrregularForce)
- Energy conservation: ✓ Pass (dE/E0 = 1.34e-6)

**Git range:** `dbbbfcc` → `03dd8dc`

**Validation results:**

| Metric | Baseline | Optimized | Result |
|--------|----------|-----------|--------|
| IrregularForce | 69.5s | 67.3s | -3.2% |
| Wall time | 130s | 126.4s | -2.8% |
| dE/E0 | — | 1.34e-6 | ✓ Pass |

**Notes:** AVX-512 vectorization achieved modest gain limited by gather overhead and memory-bound workload. Profiling infrastructure provides foundation for future optimization work.

---

## v1.0 AoS to SoA Conversion (Shipped: 2026-01-17)

**Delivered:** Complete conversion of ABYSS particle data structures from Array of Structures (AoS) to Structure of Arrays (SoA) layout, maintaining full simulation correctness.

**Phases completed:** 1-7 (28 plans total)

**Key accomplishments:**

- Created ParticleData SoA container with 66 array pointers and accessor functions
- Implemented ParticleDataMPI with 66 MPI_Win handles for shared memory
- Created ParticleDataGPU device container with async transfer support
- Added SoA helper functions for AoS↔SoA bidirectional sync
- Maintained SDAR few-body compatibility via exit-only sync pattern
- Updated I/O to sync SoA after initialization and checkpoint restore
- Validated energy conservation matches baseline (2.52e-5 vs 3.21e-5)

**Stats:**

- 7 phases, 28 plans
- 57 commits on AoS_to_SoA branch
- Energy conservation: ✓ Pass (better than baseline)
- Performance: No measurable change (28.4s vs 28.3s baseline)

**Git range:** `feat(01-01)` → `feat(06-02)`

**Validation results:**

| Metric | Baseline | SoA | Result |
|--------|----------|-----|--------|
| dE/E0 mean | 3.21e-05 | 2.52e-05 | ✓ Pass |
| Wall time | 28.33s | 28.40s | No change |

**Notes:** No performance improvement observed. Bottleneck likely in GPU compute or MPI communication rather than CPU memory access patterns. SoA layout provides foundation for future SIMD optimizations.

---
