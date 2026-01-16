# ABYSS Optimization Summary & Plan (Codex)

## Executive Summary
This plan targets performance, scalability, and code quality improvements across MPI, OpenMP, and CUDA while **not touching SDAR/SEVN**. The biggest wins are in GPU kernel efficiency, GPU/CPU overlap, and MPI scheduling. Refactors are scoped to non-SDAR/SEVN code paths only, primarily in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu), [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu), [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp), [src/MPIRoutines.cpp](src/MPIRoutines.cpp), [src/IrregularRoutines.cpp](src/IrregularRoutines.cpp), and [src/RegularRoutines.cpp](src/RegularRoutines.cpp).

### Primary Goals Mapped
- **Performance / latency / memory:** GPU kernel tuning, overlap H2D/D2H with compute, eliminate repeated allocations, improve neighbor reduction.
- **Code quality / architecture:** consolidate repeated MPI send logic, introduce reusable GPU buffers, add clear data-flow staging.
- **Scalability:** non-blocking MPI task distribution, better GPU multi-stream utilization, OpenMP in CPU update loops.
- **MPI/OpenMP/CUDA:** explicit actions listed below.

---

## Observations (Current Bottlenecks)
1. **Blocking MPI task assignment** in [src/MPIRoutines.cpp](src/MPIRoutines.cpp) (`InitialAssignmentOfTasks` overloads) limits overlap and adds latency at scale.
2. **GPU work is chunked**, but host-side merging and neighbor list reductions in [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu) are mostly serial per target; compute/transfer overlap is limited.
3. **CUDA kernel `compute_forces`** in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu) uses fixed `BatchSize`, `GridDimY`, and unroll factor; Ampere/A100 has room for higher occupancy and shared memory reuse.
4. **`reduce_forces_kernel`** sums sequentially over `n` for each column, leaving reduction parallelism untapped.
5. **sendAllParticlesToGPU** in [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp) allocates vectors each call; this increases allocation overhead and pressure on the allocator.
6. **CPU update loops** in irregular/regular phases are serial and can be safely parallelized with OpenMP in most cases.

---

## Specific Optimization Recommendations (with rationale)

### A) MPI Optimizations
1. **Non-blocking task assignment** in [src/MPIRoutines.cpp](src/MPIRoutines.cpp):
   - Convert `MPI_Send` to `MPI_Isend` with `MPI_Waitall` or `MPI_Test` for overlap.
   - Expected impact: **5–15%** in high-processor runs with large task counts.

2. **Batch queue sends** (when multiple tasks):
   - Group tasks into a single send for each worker to reduce per-message overhead.
   - Expected impact: **3–8%** depending on task granularity.

### B) OpenMP Optimizations
1. **Parallelize particle update loops** in irregular and regular updates:
   - Example: irregular update loop in [src/IrregularRoutines.cpp](src/IrregularRoutines.cpp) and any per-particle update loops in [src/RegularRoutines.cpp](src/RegularRoutines.cpp).
   - Expected impact: **20–40%** in CPU-bound phases on multi-core nodes.

2. **Thread-safe data structures**: where shared data is modified, use per-thread buffers or reductions.

### C) CUDA Optimizations
1. **Tune kernel launch parameters** in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu):
   - Increase `BatchSize` (e.g., 128) and `GridDimY` (e.g., 64) on A100.
   - Increase unroll factor for the inner loop.
   - Expected impact: **20–40%** on GPU compute kernels.

2. **Parallel reduction of `reduce_forces_kernel`**:
   - Replace the serial for-loop over `n` with a block-level reduction and/or use CUB reduction.
   - Expected impact: **10–25%** in regular acceleration reduction phase.

3. **Persistent GPU buffers and pinned host memory** in [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu):
   - Allocate once per max capacity (`J_capacity`, `I_capacity`) and reuse.
   - Use pinned host buffers for `h_result`, `h_neighbor`, `h_neighbor_count`.
   - Expected impact: **10–20%** in H2D/D2H overhead.

4. **Overlap transfers with compute**:
   - Use per-GPU streams and async copies (`cudaMemcpyAsync`) to overlap H2D/D2H with kernel execution.
   - Expected impact: **10–20%** depending on data size.

### D) Memory / Data Structure Optimizations
1. **Reuse `std::vector` capacity** in [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp):
   - Avoid repeated allocations for `Jparticles`, `Iparticles`, and index arrays.
   - Expected impact: **5–10%** on large particle counts.

2. **Prefer contiguous arrays for frequent access** (e.g., convert `std::unordered_set` to `std::vector` when stable ordering is acceptable).

---

## Refactored Code for Critical Improvements (examples)

### 1) MPI: Non-blocking task assignment
Target: [src/MPIRoutines.cpp](src/MPIRoutines.cpp)

```cpp
// Before (blocking):
for (int i = 0; i < NumberOfWorker; i++) {
    if (i >= NumTask) break;
    MPI_Send(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD);
}

// After (non-blocking):
std::vector<MPI_Request> reqs;
reqs.reserve(std::min(NumberOfWorker, NumTask));
for (int i = 0; i < NumberOfWorker; i++) {
    if (i >= NumTask) break;
    MPI_Request req;
    MPI_Isend(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &req);
    reqs.push_back(req);
}
if (!reqs.empty()) {
    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
}
```

### 2) OpenMP: Parallel irregular update loop
Target: [src/IrregularRoutines.cpp](src/IrregularRoutines.cpp)

```cpp
// Before (serial):
for (int ptcl_id : ThisLevelNode->ParticleList) {
    Particle* ptcl = &particles[ptcl_id];
    if (ptcl->NumberOfNeighbor != 0)
        ptcl->updateParticle();
    ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
    ptcl->CurrentTimeIrr = ptcl->CurrentBlockIrr * time_step;
}

// After (parallel):
#pragma omp parallel for schedule(static)
for (size_t i = 0; i < ThisLevelNode->ParticleList.size(); i++) {
    int ptcl_id = ThisLevelNode->ParticleList[i];
    Particle* ptcl = &particles[ptcl_id];
    if (ptcl->NumberOfNeighbor != 0)
        ptcl->updateParticle();
    ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
    ptcl->CurrentTimeIrr = ptcl->CurrentBlockIrr * time_step;
}
```

### 3) CUDA: Parallel reduction in `reduce_forces_kernel`
Target: [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu)

```cpp
// Example sketch: block-level reduction over n
__global__ void reduce_forces_kernel_opt(const CUDA_REAL* diff, CUDA_REAL* result, int n, int m) {
    int col = blockIdx.x;        // one block per column
    int comp = threadIdx.y;      // component (0..5)
    int tid  = threadIdx.x;      // reduction thread
    extern __shared__ CUDA_REAL s[];

    if (col >= m || comp >= 6) return;

    CUDA_REAL sum = 0.0;
    for (int i = tid; i < n; i += blockDim.x) {
        sum += diff[comp * m * n + (col * n) + i];
    }
    s[comp * blockDim.x + tid] = sum;
    __syncthreads();

    // parallel reduction
    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (tid < offset) {
            s[comp * blockDim.x + tid] += s[comp * blockDim.x + tid + offset];
        }
        __syncthreads();
    }
    if (tid == 0) {
        result[col * 6 + comp] = s[comp * blockDim.x];
    }
}
```

### 4) CUDA: Reuse buffers and pinned host memory
Target: [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu)

```cpp
// Allocate once (init path)
if (!gpu[i].h_result) {
    cudaHostAlloc(&gpu[i].h_result, _six * I_capacity * sizeof(CUDA_REAL), cudaHostAllocDefault);
    cudaHostAlloc(&gpu[i].h_neighbor, I_capacity * MaxNumNeighbor * sizeof(int), cudaHostAllocDefault);
    cudaHostAlloc(&gpu[i].h_neighbor_count, I_capacity * sizeof(int), cudaHostAllocDefault);
}

// Reuse in each call
toHost(gpu[i].h_result, gpu[i].d_result, _six * NumTarget, gpu[i].stream);
```

---

## Performance Impact Estimates (high-level)
- **MPI non-blocking sends:** 5–15% latency reduction for high-rank runs.
- **OpenMP in update loops:** 20–40% speedup in CPU-bound segments.
- **CUDA kernel tuning (BatchSize/GridDimY/unroll):** 20–40% GPU kernel speedup.
- **Parallel reduction kernel:** 10–25% GPU reduction speedup.
- **Pinned memory + overlap:** 10–20% end-to-end improvement in GPU phases.

These are **conservative** and depend on particle count, GPU occupancy, and interconnect characteristics.

---

## Implementation Plan (Phased)

### Phase 0 (Baseline & Safety)
- Enable detailed profiling (already in place: `PERFORMANCETRACE`, NSIGHT markers).
- Record baseline timings from existing `profiling.csv`.

### Phase 1 (Low Risk / High Reward)
1. Non-blocking MPI assignment in [src/MPIRoutines.cpp](src/MPIRoutines.cpp).
2. Reuse allocations in [src/cuda/CalculateRegularAcceleration.cpp](src/cuda/CalculateRegularAcceleration.cpp).
3. Pinned host memory in [src/cuda/cuda_my_acceleration.cu](src/cuda/cuda_my_acceleration.cu).

### Phase 2 (OpenMP + GPU Kernel Tuning)
1. Add OpenMP pragma to irregular and regular update loops in [src/IrregularRoutines.cpp](src/IrregularRoutines.cpp) and [src/RegularRoutines.cpp](src/RegularRoutines.cpp).
2. Tune `BatchSize`, `GridDimY`, and unroll in [src/cuda/cuda_kernels.cu](src/cuda/cuda_kernels.cu).

### Phase 3 (Deep CUDA Optimizations)
1. Replace `reduce_forces_kernel` with a parallel reduction or CUB.
2. Stream overlap for H2D/D2H and kernel launch pipelines.
3. Optional: add persistent L2 cache window on A100 (advanced).

### Phase 4 (Scalability & Architecture)
1. Investigate data structure choices for `RegularList` / `RegularMap` to reduce overhead at large N.
2. Evaluate task batching to reduce MPI message counts.

---

## Notes / Constraints
- **SDAR/SEVN code is untouched** by all recommendations above.
- Avoid `-ffast-math` or any non-IEEE math transformations to preserve numerical precision.

---

## Next Actions (If You Want Me to Implement)
- Pick Phase 1 items to implement first; I can apply safe changes and validate with the existing profiling pipeline.
- Confirm your target GPU architecture and MPI stack for more precise tuning.
