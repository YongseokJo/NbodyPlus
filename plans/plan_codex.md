# ABYSS Optimization Plan (Codex)

## Scope and Constraints
- Do not modify any code under `SDAR/` or any `#ifdef SEVN` / `SEVN` integration paths (for example, `src/StellarEvolution.cpp` or SEVN blocks in `src/particle.h`).
- All recommendations and refactors below target non-SDAR/SEVN paths only.

## High-level Findings (Non-SDAR/SEVN)
1. MPI task assignment and queue distribution are largely blocking and synchronous (`src/MPIRoutines.cpp`, `src/QueueScheduler.h`).
2. GPU path spends significant time in host-side gather and serial reduction (`src/cuda/CalculateRegularAcceleration.cpp`, `src/cuda/cuda_kernels.cu`, `src/cuda/cuda_my_acceleration.cu`).
3. Per-step allocations and hash-based neighbor sets are pervasive in CPU routines (`src/Particle/ComputeAcceleration.cpp`, `src/RegularRoutines.cpp`).
4. CPU update loops are serial even when data-parallel (`src/IrregularRoutines.cpp`, `src/WorkerRoutines.cpp`).
5. Build flags are conservative for performance (release flags and OpenMP disabled) (`src/Makefile`).

## Optimization Recommendations

### MPI / Scheduling
- Replace blocking `MPI_Send` in `InitialAssignmentOfTasks` with `MPI_Isend` + `MPI_Waitall` to reduce root stalls and improve overlap.
  - Expected impact: 5-15% on high-rank runs with many short tasks.
- Batch task queues per worker (send arrays of queue entries) to reduce per-message overhead in `QueueScheduler` and `WorkerRoutines`.
  - Expected impact: 3-8%.
- Use `MPI_Igatherv` in `sendAllParticlesToGPU` to overlap gather with CPU work, then complete before `SendToDevice`.
  - Expected impact: 5-10% during regular GPU steps.

### CUDA / GPU
- Parallelize `reduce_forces_kernel` (currently serial over n). Use block-level reduction or CUB `DeviceReduce`.
  - Expected impact: 10-25% on GPU reduction phase.
- Tune kernel launch parameters (`BatchSize`, `GridDimY`) in `src/def.h` for target GPU; add `__launch_bounds__` in `compute_forces` to control register pressure.
  - Expected impact: 15-35% in compute kernel.
- Reuse and pin host buffers for `Jparticles`, `Iparticles`, and `RegularListIndices` in `sendAllParticlesToGPU` to avoid per-step allocations.
  - Expected impact: 5-15% in host-side overhead.
- Overlap H2D/D2H copies with compute via streams and events; avoid `cudaStreamSynchronize` until the chunk is fully staged.
  - Expected impact: 10-20% in GPU pipeline.

### OpenMP / CPU
- Parallelize pure per-particle loops in irregular update and regular update stages (`src/IrregularRoutines.cpp`, `src/WorkerRoutines.cpp`) with `#pragma omp parallel for`.
  - Expected impact: 20-40% on CPU-heavy phases.
- Replace `std::unordered_set` neighbor caches in `computeAccelerationReg` and `updateRegularParticleCuda` with a reusable `std::vector<int>` plus a scratch bitset (per thread) to reduce allocations.
  - Expected impact: 5-15% for CPU regular/neighbor updates.

### Memory / Data Layout
- Consider a SoA layout for `Iparticle`/`Jparticle` on GPU (separate arrays for pos/vel/mass). This improves coalescing and allows vectorized loads.
  - Expected impact: 10-30% on memory-bound kernels (medium effort).
- Pre-allocate and reuse `RegularList` containers to reduce churn across steps.

### Build / Tooling
- Provide a release build profile: `-O3 -DNDEBUG -fopenmp -flto` (and remove `-g` in release builds).
  - Expected impact: 10-20% on CPU-heavy runs.
- For CUDA builds, set `-O3 -arch=sm_80` (or target arch) and keep PTXAS verbosity for register tuning.

## Refactored Code (Critical Improvements)

### 1) MPI non-blocking task assignment
Target: `src/MPIRoutines.cpp`

```cpp
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG) {
    int count = std::min(NumberOfWorker, NumTask);
    std::vector<MPI_Request> reqs(count);
    for (int i = 0; i < count; i++) {
        MPI_Isend(&data[i], 1, MPI_INT, i + 1, TAG, MPI_COMM_WORLD, &reqs[i]);
    }
    if (count > 0) {
        MPI_Waitall(count, reqs.data(), MPI_STATUSES_IGNORE);
    }
}
```

### 2) CUDA reduction kernel parallelization
Target: `src/cuda/cuda_kernels.cu`

```cpp
__global__ void reduce_forces_kernel_opt(const CUDA_REAL* diff, CUDA_REAL* result, int n, int m) {
    int col = blockIdx.x;
    int comp = threadIdx.y; // 0..5
    int tid = threadIdx.x;
    extern __shared__ CUDA_REAL sdata[];

    if (col >= m || comp >= 6) return;

    CUDA_REAL sum = 0.0;
    int base = comp * m * n + col * n;
    for (int i = tid; i < n; i += blockDim.x) {
        sum += diff[base + i];
    }
    sdata[comp * blockDim.x + tid] = sum;
    __syncthreads();

    for (int offset = blockDim.x / 2; offset > 0; offset >>= 1) {
        if (tid < offset) {
            sdata[comp * blockDim.x + tid] += sdata[comp * blockDim.x + tid + offset];
        }
        __syncthreads();
    }
    if (tid == 0) {
        result[col * 6 + comp] = sdata[comp * blockDim.x];
    }
}
```

### 3) OpenMP irregular update loop
Target: `src/IrregularRoutines.cpp`

```cpp
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

### 4) Reuse gather buffers for GPU staging
Target: `src/cuda/CalculateRegularAcceleration.cpp`

```cpp
// Sketch: move these to static or a dedicated staging struct, reuse capacity
static std::vector<Jparticle> Jparticles;
static std::vector<Iparticle> Iparticles;
static std::vector<int> RegularListIndices;
static std::vector<int> counts, Jcounts, Icounts, Jdispls, Idispls;

Jparticles.resize(NumberOfParticle);
Iparticles.resize(RegularListSize);
RegularListIndices.resize(RegularListSize);
counts.resize(NumberOfProcessor * 2);
Jcounts.resize(NumberOfProcessor);
Icounts.resize(NumberOfProcessor);
Jdispls.resize(NumberOfProcessor);
Idispls.resize(NumberOfProcessor);
```

## Performance Impact Estimates (Conservative)
- MPI non-blocking + batching: 5-15% on multi-rank runs with many short tasks.
- CUDA reduction + kernel tuning: 20-45% in GPU compute and reduction phases.
- Buffer reuse + overlap: 10-20% in GPU pipeline wall time.
- OpenMP updates: 20-40% in CPU-heavy phases.
- Build flags (release): 10-20% on CPU-only or mixed runs.

## Phased Execution Plan
1. Phase 0: Profiling baseline using existing `PERFORMANCETRACE` and NVTX markers; capture wall time and GPU kernel stats.
2. Phase 1 (low risk): MPI non-blocking sends, buffer reuse for GPU staging, release build flags.
3. Phase 2 (medium): OpenMP on update loops, CUDA reduction kernel, stream overlap.
4. Phase 3 (high impact): SoA data layout for GPU, queue scheduler refactors, GPU resident particle arrays.
