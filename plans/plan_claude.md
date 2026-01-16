# ABYSS Codebase Optimization Plan

## Executive Summary
ABYSS is an N-body astrophysical simulation code using a Hermite integrator with block time-steps. It employs an MPI master-worker model with optional CUDA GPU acceleration for regular force calculations. This plan outlines optimizations across MPI, CUDA, compiler settings, and code quality while respecting the constraint that **SDAR and SEVN code should not be touched**.

### User Preferences (Clarified)
- **Priority:** CUDA/GPU performance optimization
- **Target GPU:** NVIDIA Ampere (A100, A30) - Compute Capability 8.0+
- **Precision:** Strict numerical precision required (no `-ffast-math`)

---

## Quick Reference: Key Optimizations

| Priority | Optimization | Expected Impact | Effort |
|----------|-------------|-----------------|--------|
| **1** | **CUDA kernel tuning for A100** | **30-50% on GPU path** | **Medium** |
| **2** | **CUDA memory/stream optimization** | **15-25% on GPU path** | **Medium** |
| **3** | **NVCC flags for Ampere sm_80** | **5-15%** | **Low** |
| 4 | Compiler flags (-O3, -fopenmp, -flto) | 10-20% | Low |
| 5 | Enable MULTIMAP option | 10-20% for large N | Low |
| 6 | OpenMP parallelization of update loops | 20-40% on multi-core | Medium |
| 7 | MPI non-blocking sends | 5-15% | Medium |
| 8 | Memory optimization | 5-10% | High |

---

## 1. Compiler Optimization Enhancements

### 1.1 Makefile Improvements (`src/Makefile`)

**Current:** `-O2 -march=native -ftree-vectorize`

**Recommended Changes:**
- Upgrade to `-O3` for more aggressive optimizations
- Add `-funroll-loops` for loop unrolling
- **DO NOT use `-ffast-math`** (user requires strict numerical precision)
- Add OpenMP support (`-fopenmp`) which is currently commented out
- Add `-flto` for link-time optimization

**Impact:** 10-20% performance improvement in CPU-bound sections

---

## 2. CUDA Optimizations (HIGH PRIORITY)

### 2.1 Kernel Tuning for Ampere A100 (`src/cuda/cuda_kernels.cu`)

**File:** `src/cuda/cuda_kernels.cu:28-114`

**Current Configuration:**
- `BatchSize=64`, `GridDimY=32`, `NNB_per_block=128`

**Ampere-Specific Optimizations:**

1. **Increase `BatchSize` to 128 or 256** - A100 has 192KB shared memory per SM (vs 96KB on Volta)
   ```cpp
   // In def.h, change:
   #define BatchSize 128  // Was 64, A100 can handle larger batches
   #define GridDimY 64    // Was 32, better occupancy on A100
   ```

2. **Add `__launch_bounds__`** for better register allocation:
   ```cpp
   __global__ __launch_bounds__(128, 8)  // maxThreads=128, minBlocks=8
   void compute_forces(...)
   ```

3. **Enable asynchronous memory copies with `cudaMemcpyAsync`** using CUDA streams (partially done, enhance further)

4. **Use L2 cache persistence hints** (A100 feature):
   ```cpp
   cudaStreamAttrValue attr;
   attr.accessPolicyWindow.base_ptr = d_Jp;
   attr.accessPolicyWindow.num_bytes = J_count * sizeof(Jparticle);
   attr.accessPolicyWindow.hitRatio = 1.0f;
   attr.accessPolicyWindow.hitProp = cudaAccessPropertyPersisting;
   cudaStreamSetAttribute(stream, cudaStreamAttributeAccessPolicyWindow, &attr);
   ```

5. **Increase unroll factor**:
   ```cpp
   #pragma unroll 8  // Was 4, more aggressive for A100
   for (int jj=0; jj<current_batch_size; jj++){
   ```

**Expected Impact:** 30-50% improvement on GPU-bound operations

### 2.2 Memory Transfer Optimization (`src/cuda/cuda_my_acceleration.cu`)

**File:** `src/cuda/cuda_my_acceleration.cu:250-341`

**Issues:**
- Separate H2D transfers for J and I particles
- GPU data structures are reallocated when particle count changes

**Recommendations:**
1. Use CUDA streams more aggressively for overlapping H2D/D2H transfers with computation
2. Implement persistent GPU memory allocation to avoid reallocation overhead
3. Use pinned memory consistently (already partially done)

### 2.3 Reduction Kernel (`src/cuda/cuda_kernels.cu:336-415`)

**Current:** Good warp-level reduction using `__shfl_down_sync`

**Enhancement:**
- The `reduce_forces_kernel` at line 117-144 loops serially over rows - consider parallel reduction across the n dimension

### 2.4 NVCC Compilation Flags for A100

Update in `src/Makefile`:
```makefile
# Current NVCC flags
NVCCFLAGS = -O2 -std=c++11 ...

# Optimized for A100 (sm_80)
NVCCFLAGS = -O3 -std=c++11 \
    -arch=sm_80 \
    --use_fast_math=false \
    -Xcompiler -O3 \
    -maxrregcount=64 \
    --ptxas-options=-v \
    -Xptxas -dlcm=ca
```

---

## 3. MPI Communication Optimizations

### 3.1 Replace Blocking Sends with Non-blocking (`src/MPIRoutines.cpp`)

**File:** `src/MPIRoutines.cpp:73-119`

**Issue:** Multiple `InitialAssignmentOfTasks` functions use blocking `MPI_Send`

**Recommendation:** Convert to `MPI_Isend` + `MPI_Waitall` pattern:
```cpp
// Current (blocking)
for (int i=0; i<NumberOfWorker; i++) {
    MPI_Send(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD);
}

// Optimized (non-blocking)
std::vector<MPI_Request> requests(NumberOfWorker);
for (int i=0; i<NumberOfWorker; i++) {
    MPI_Isend(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[i]);
}
MPI_Waitall(NumberOfWorker, requests.data(), MPI_STATUSES_IGNORE);
```

**Impact:** Reduced synchronization overhead, better overlap of communication

### 3.2 Optimize Collective Operations (`src/cuda/CalculateRegularAcceleration.cpp`)

**File:** `src/cuda/CalculateRegularAcceleration.cpp:174-272`

**Current:** Uses `MPI_Gather` + `MPI_Gatherv` for particle data collection

**Recommendation:** This is already well-optimized. Consider using `MPI_Igatherv` for async gathering if computation can overlap.

---

## 4. OpenMP Parallelization (New)

### 4.1 Enable OpenMP in Makefile
Uncomment `-fopenmp` in `CXXFLAGS` and add OpenMP pragmas to key loops.

### 4.2 Parallelize Irregular Update Loop (`src/IrregularRoutines.cpp`)

**File:** `src/IrregularRoutines.cpp:264-271`

**Current (Sequential):**
```cpp
for (int ptcl_id : ThisLevelNode->ParticleList) {
    ptcl = &particles[ptcl_id];
    if (ptcl->NumberOfNeighbor != 0)
        ptcl->updateParticle();
    ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
    ptcl->CurrentTimeIrr = ptcl->CurrentBlockIrr*time_step;
}
```

**Optimized:**
```cpp
#pragma omp parallel for
for (size_t i = 0; i < ThisLevelNode->ParticleList.size(); i++) {
    int ptcl_id = ThisLevelNode->ParticleList[i];
    Particle* ptcl = &particles[ptcl_id];
    if (ptcl->NumberOfNeighbor != 0)
        ptcl->updateParticle();
    ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
    ptcl->CurrentTimeIrr = ptcl->CurrentBlockIrr*time_step;
}
```

### 4.3 Parallelize Regular Time Update (`src/RootRoutines.cpp`)

**File:** `src/RootRoutines.cpp:186-215` (`updateNextRegTime` function)

**Current:** Serial loop over all particles to find minimum time

**Optimized:** Use OpenMP reduction for parallel minimum finding:
```cpp
#pragma omp parallel
{
    ULL local_time = block_max;
    std::vector<int> local_list;

    #pragma omp for nowait
    for (int i = 0; i <= LastParticleIndex; i++) {
        Particle* ptcl = &particles[i];
        if (!ptcl->isActive) continue;
        ULL time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;
        if (time_tmp < local_time) {
            local_list.clear();
            local_time = time_tmp;
            local_list.push_back(ptcl->ParticleIndex);
        } else if (time_tmp == local_time) {
            local_list.push_back(ptcl->ParticleIndex);
        }
    }

    #pragma omp critical
    {
        if (local_time < time) {
            time = local_time;
            RegularList.clear();
            RegularList.insert(local_list.begin(), local_list.end());
        } else if (local_time == time) {
            RegularList.insert(local_list.begin(), local_list.end());
        }
    }
}
```

---

## 5. Memory Optimization

### 5.1 Reduce Dynamic Allocations in Acceleration Calculations

**File:** `src/Particle/ComputeAcceleration.cpp:43, 256, 430`

**Issue:** `unordered_set` is created inside hot functions:
```cpp
std::unordered_set<int> CMPtclsSet;           // Line 43
std::unordered_set<int> RealNeighbors;        // Line 256
std::unordered_set<int> hashTableOld;         // Line 430
std::unordered_set<int> hashTableNew;         // Line 432
```

**Recommendation:**
1. Reserve capacity upfront: `.reserve(NumberOfNeighbor)` is done but consider using thread-local storage
2. Consider using `std::vector<bool>` for neighbor lookup if indices are bounded

### 5.2 Optimize Particle Structure (`src/particle.h`)

**Current Size:** Very large struct with many fields

**Recommendation:** Separate hot and cold data:
- Create `ParticleHot` struct with frequently accessed fields (Position, Velocity, Mass, a_tot)
- Create `ParticleCold` struct with rarely accessed fields (FormationTime, WorldTime, etc.)

---

## 6. Algorithm Optimizations

### 6.1 Use Priority Queue for Regular Time Scheduling

**File:** `src/RootRoutines.cpp:186-215`

**Issue:** Linear scan O(N) to find particles with minimum next regular time

**Recommendation:** Use `std::priority_queue` or keep sorted structure:
```cpp
// Already have MULTIMAP option which uses std::multimap - enable it!
// Add -DMULTIMAP to CXXFLAGS
```

The `MULTIMAP` option is already implemented but disabled. Enable it for O(log N) lookups.

### 6.2 Skip List Optimization (`src/SkipList.h`)

Review the skip list implementation for memory locality optimizations.

---

## 7. Code Quality Improvements

### 7.1 Remove Dead Code

**Files with `#ifdef unuse` blocks:**
- `src/cuda/cuda_routines.cu:13-128`
- `src/QueueScheduler.h:334-552`

**Action:** Remove unused code blocks to improve maintainability.

### 7.2 Consolidate Duplicate Functions

**File:** `src/TimeStepRoutines.cpp`

**Issue:** `getNewTimeStepIrr` and `getNewTimeStep` are nearly identical

**Action:** Consolidate into single templated function.

---

## Files to Modify

1. `src/Makefile` - Compiler flags, NVCC flags
2. `src/def.h` - CUDA parameters (BatchSize, GridDimY)
3. `src/cuda/cuda_kernels.cu` - CUDA kernel tuning
4. `src/cuda/cuda_my_acceleration.cu` - Memory management
5. `src/IrregularRoutines.cpp` - OpenMP parallelization
6. `src/RootRoutines.cpp` - OpenMP parallelization, MULTIMAP
7. `src/MPIRoutines.cpp` - Non-blocking MPI
8. `src/Particle/ComputeAcceleration.cpp` - Memory optimization

**Files NOT to modify:** Anything under `SDAR/` or `SEVN/`

---

## Verification Plan

1. **Build verification:** `make clean && make USE_CUDA=1`
2. **Unit test:** Run `tests/test1` to verify correctness
3. **Performance test:** Compare wall-clock time on test_10 or test_20
4. **Numerical validation:** Compare output energy error with baseline
