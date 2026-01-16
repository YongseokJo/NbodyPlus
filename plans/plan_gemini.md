# Optimization Plan for ABYSS Codebase

## 1. Specific Optimization Recommendations

### A. GPU Memory Access optimization (Critical)
**Issue:** The current CUDA implementation uses an Array of Structures (AoS) layout for `Jparticle` (background particles).
- `Jparticle` size is 64 bytes (stride = 64).
- CUDA threads in a warp access memory in a strided pattern (Thread 0 reads address X, Thread 1 reads X+64, etc.), leading to **uncoalesced memory accesses**. This wastes ~75-87% of global memory bandwidth.
**Solution:** Switch to **Structure of Arrays (SoA)** for `Jparticle` in global memory.
- Store `posx`, `posy`, `posz`, `mass`, etc., in separate contiguous arrays.
- Threads loading data into shared memory will perform fully coalesced reads (32 x 4 bytes contiguous).
**Refactor:**
- Define `JparticleSoA` in `src/cuda/cuda_defs.h`.
- Update `compute_forces` and `CalculateAccelerationOnDevice` to specific SoA pointers.
- Perform AoS -> SoA conversion on the host before sending to GPU (or use a transpose kernel).

### B. MPI Scalability & Latency
**Issue:** `InitialAssignmentOfTasks` uses blocking `MPI_Send`.
- This serializes the dispatch process. Use `MPI_Isend` to dispatch tasks asynchronously.
**Solution:** Check `src/MPIRoutines.cpp` and ensure `MPI_Isend` is used properly with `MPI_Waitall` or similar management to avoid buffer reuse issues, or valid "fire and forget" if buffers are separate.

### C. Compiler Optimizations
**Issue:** `Makefile` uses `-O2` and `-std=c++11`.
**Solution:**
- Upgrade to `-O3` for aggressive loop optimizations.
- Add `-flto` (Link Time Optimization) to inline functions across translation units.
- Use `-march=native` (already present) but verify it's effective on the cluster.
- **CUDA:** Ensure `-use_fast_math` is enabled for `nvcc` to use hardware approximations for `rsqrt`, `sin`, `cos`.

### D. Algorithmic / Code Quality
- **Branch Divergence:** The `if (d2 > Ip.r2)` check in the kernel causes warp divergence.
    - *Mitigation:* Sort `J` particles spatially (e.g., simple coordinate sort on CPU before sending) to maximize the chance that a warp process particles either all "inside" or all "outside" the radius.
- **I/O:** `main.cpp` calls `fflush` frequently. This can be slow on network filesystems. Remove `fflush` from inner loops or tight paths.

---

## 2. Refactored Code (Critical Improvements)

### Step 1: Define SoA Structure (`src/cuda/cuda_defs.h`)

```cpp
// Add to src/cuda/cuda_defs.h

// Structure of Arrays for J particles (Optimized for Global Memory Coalescing)
struct JparticleSoA {
    CUDA_REAL* posx;
    CUDA_REAL* posy;
    CUDA_REAL* posz;
    CUDA_REAL* mass;
    CUDA_REAL* velx;
    CUDA_REAL* vely;
    CUDA_REAL* velz;
    int* index;
};
```

### Step 2: Optimized CUDA Kernel (`src/cuda/cuda_kernels.cu`)

Updated to accept individual arrays (SoA pointers) and load them into shared memory (AoS) for efficient processing.

```cuda-cpp
// src/cuda/cuda_kernels.cu

__global__ void compute_forces_soa(
    const Iparticle* __restrict__ d_Ip, 
    // Jparticle SoA pointers
    const CUDA_REAL* __restrict__ j_posx,
    const CUDA_REAL* __restrict__ j_posy,
    const CUDA_REAL* __restrict__ j_posz,
    const CUDA_REAL* __restrict__ j_mass,
    const CUDA_REAL* __restrict__ j_velx,
    const CUDA_REAL* __restrict__ j_vely,
    const CUDA_REAL* __restrict__ j_velz,
    const int*       __restrict__ j_index,
    
    CUDA_REAL* __restrict__ acc, 
    int* __restrict__ neighbor, 
    int* num_neighbor,
    int m, int n, int i_start
) {
    // ... [Setup code same as before] ...
    
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int tid = threadIdx.x;
    // ...

    for (int j = j_begin; j < j_end; j += BatchSize) {
        int current_batch_size = min(BatchSize, j_end - j);
        
        // Coalesced Load: Global (SoA) -> Shared (AoS)
        // Each thread reads one scalar from contiguous memory -> 100% Bandwidth Efficiency
        __shared__ Jparticle Jp_sh[BatchSize]; 

        if (tid < current_batch_size) {
            int global_idx = j + tid;
            Jp_sh[tid].posx = j_posx[global_idx];
            Jp_sh[tid].posy = j_posy[global_idx];
            Jp_sh[tid].posz = j_posz[global_idx];
            Jp_sh[tid].mass = j_mass[global_idx];
            Jp_sh[tid].velx = j_velx[global_idx];
            Jp_sh[tid].vely = j_vely[global_idx];
            Jp_sh[tid].velz = j_velz[global_idx];
            Jp_sh[tid].index = j_index[global_idx];
        }
        __syncthreads();

        // ... [Computation loop same as before, using Jp_sh] ...
    }
}
```

### Step 3: Implement Host-Side SoA Conversion (`src/cuda/CalculateRegularAcceleration.cpp`)

Modify `sendAllParticlesToGPU` (or the function calling `SendToDevice`) to repack data.

```cpp
// Inside src/cuda/CalculateRegularAcceleration.cpp or helper

void SendToDeviceSoA(const std::vector<Jparticle>& host_J_aos, const std::vector<Iparticle>& host_I) {
    size_t n = host_J_aos.size();
    
    // 1. Allocate SoA buffers on Host (temporary)
    std::vector<CUDA_REAL> h_jx(n), h_jy(n), h_jz(n), h_mass(n);
    std::vector<CUDA_REAL> h_jvx(n), h_jvy(n), h_jvz(n);
    std::vector<int>       h_jidx(n);

    // 2. Transpose AoS -> SoA (CPU Cache friendly if done linearly)
    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        h_jx[i]   = host_J_aos[i].posx;
        h_jy[i]   = host_J_aos[i].posy;
        h_jz[i]   = host_J_aos[i].posz;
        h_mass[i] = host_J_aos[i].mass;
        h_jvx[i]  = host_J_aos[i].velx;
        h_jvy[i]  = host_J_aos[i].vely;
        h_jvz[i]  = host_J_aos[i].velz;
        h_jidx[i] = host_J_aos[i].index;
    }

    // 3. Send to Device (Allocate device pointers once and reuse if possible)
    // Assume d_jx, d_jy... are global or static device pointers
    
    // Use cudaMemcpyAsync with streams
    cudaMemcpyAsync(d_jx, h_jx.data(), n * sizeof(CUDA_REAL), cudaMemcpyHostToDevice, stream);
    // ... repeat for all fields ...
    
    // 4. Launch Kernel with SoA pointers
    // compute_forces_soa<<<..., stream>>>(d_Ip, d_jx, d_jy, ..., n, ...);
}
```

### Step 4: MPI Optimization (`src/MPIRoutines.cpp`)

```cpp
// src/MPIRoutines.cpp - Optimized Task Assignment

void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG) {
    // Static request array to avoid reallocation
    static std::vector<MPI_Request> requests(NumberOfWorker * 2); 
    int req_count = 0;

    for (int i = 0; i < NumberOfWorker; i++) {
        if (i >= NumTask) break;
        
        // Use Non-Blocking Send
        MPI_Isend(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD, &requests[req_count++]);
        MPI_Isend(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD, &requests[req_count++]);
    }
    
    // Wait for all sends to complete (or manage buffer life-cycle properly)
    MPI_Waitall(req_count, requests.data(), MPI_STATUSES_IGNORE);
}
```

## 3. Performance Impact Estimates

| Optimization | Estimated Impact | Explanation |
| :--- | :--- | :--- |
| **SoA Memory Layout** | **2.0x - 4.0x Speedup** (Kernel) | Coalesced loads utilize ~100% memory bandwidth vs ~12% with stride-64 AoS. |
| **MPI`Isend`** | **10 - 30%** (Worker Dispatch) | Reduces idle time for workers waiting for tasks; critical for large core counts. |
| ** `-O3 -flto` ** | **5 - 15%** (Host Code) | Better vectorization and inlining for CPU-bound routines. |
| **Fast Math** | **1.2x - 1.5x** (Kernel) | Hardware intrinsics for `rsqrt` and trigonometry. |
