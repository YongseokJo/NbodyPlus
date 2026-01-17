# Phase 3 Research: GPU Integration

## Overview

This phase converts the GPU data structures and kernels from AoS to SoA layout. The GPU code handles regular force calculations in an N-body simulation.

## Current GPU Architecture

### Data Structures (`src/cuda/cuda_defs.h`)

**Target particle (i-particle)** — receives forces:
```cpp
struct i_particle_t {
    cuda_real_t pos_x, pos_y, pos_z, radius_sq;
    cuda_real_t vel_x, vel_y, vel_z, dt_reg;
};  // 8 fields × 4/8 bytes
```

**Source particle (j-particle)** — exerts forces:
```cpp
struct j_particle_t {
    cuda_real_t pos_x, pos_y, pos_z, mass;
    cuda_real_t vel_x, vel_y, vel_z;
    int index;
    int pad;  // alignment padding (double precision only)
};  // 9 fields
```

### Kernel Inventory (`src/cuda/cuda_kernels.cu`, `cuda_kernels.h`)

| Kernel | Signature | Purpose |
|--------|-----------|---------|
| `compute_forces_kernel` | `(i_particle_t*, j_particle_t*, result*, ...)` | Main force calculation |
| `reduce_forces_kernel` | `(diff*, result*, n, m)` | Warp-shuffle reduction |
| `gather_neighbor_kernel` | `(neighbor_block*, num_neighbor*, ...)` | Neighbor list gathering |
| `gather_numneighbor_kernel` | `(numneighbor_block*, ...)` | Neighbor count reduction |
| `initialize_arrays` | `(result*, diff*, n, m, subset*)` | Array zeroing |

### Host-Device Transfer (`src/cuda/cuda_acceleration.cu`)

**GPU struct** holds device pointers:
```cpp
struct GPU {
    int id;
    cudaStream_t stream;
    j_particle_t* dJ;
    i_particle_t* dI;
    cuda_real_t* d_result;
    int* d_neighbor;
    // ... pinned host buffers
};
```

**Transfer function**: `_ReceiveFromHost(vector<j_particle_t>& hJ, vector<i_particle_t>& hI)`
- Partitions j-particles across multiple GPUs (`block_range`)
- Replicates i-particles to all GPUs
- Uses async transfers with streams

**Result retrieval**: In `GetAcceleration()`:
- Results written to `ptcl->acc_irregular[dim][order]`
- Neighbor lists written to global `new_neighbors` array

### Multi-GPU Support

- `deviceCount` GPUs supported
- Each GPU gets a slice of j-particles
- i-particles replicated to all GPUs
- Results accumulated from all GPUs

### Shared Memory Usage

`compute_forces_kernel` uses shared memory for j-particles:
```cpp
__shared__ j_particle_t Jp_sh[BATCH_SIZE];
```

This tiling pattern loads `BATCH_SIZE` j-particles into shared memory for reuse across i-particle calculations.

## SoA Conversion Strategy

### Option A: SoA Device Pointers (Recommended)

Replace struct arrays with separate device arrays:

```cpp
struct GPUSoA {
    // j-particles (source)
    cuda_real_t* d_j_pos_x;
    cuda_real_t* d_j_pos_y;
    cuda_real_t* d_j_pos_z;
    cuda_real_t* d_j_mass;
    cuda_real_t* d_j_vel_x;
    cuda_real_t* d_j_vel_y;
    cuda_real_t* d_j_vel_z;
    int* d_j_index;

    // i-particles (target)
    cuda_real_t* d_i_pos_x;
    cuda_real_t* d_i_pos_y;
    cuda_real_t* d_i_pos_z;
    cuda_real_t* d_i_radius_sq;
    cuda_real_t* d_i_vel_x;
    cuda_real_t* d_i_vel_y;
    cuda_real_t* d_i_vel_z;
    cuda_real_t* d_i_dt_reg;
};
```

**Pros:**
- Coalesced memory access in kernels
- Direct mapping to CPU SoA (`ParticleData`)
- Each stream can transfer arrays independently

**Cons:**
- More `cudaMemcpy` calls (mitigate with streams)
- Shared memory tiling needs adjustment

### Option B: Keep AoS Structs on GPU

Keep `i_particle_t` and `j_particle_t` as-is, convert at host boundary.

**Pros:**
- Minimal kernel changes
- Shared memory tiling works as-is

**Cons:**
- Requires CPU-side conversion each call
- Conversion overhead (violates Pitfall #1)

**Decision**: Use Option A (full SoA on GPU) for consistency with Phase 1 design.

### Shared Memory Adaptation

Current AoS tiling:
```cpp
__shared__ j_particle_t Jp_sh[BATCH_SIZE];
Jp_sh[tid] = d_Jp[j + tid];
// Access: Jp_sh[jj].pos_x
```

SoA tiling:
```cpp
__shared__ cuda_real_t s_pos_x[BATCH_SIZE];
__shared__ cuda_real_t s_pos_y[BATCH_SIZE];
__shared__ cuda_real_t s_pos_z[BATCH_SIZE];
__shared__ cuda_real_t s_mass[BATCH_SIZE];
__shared__ cuda_real_t s_vel_x[BATCH_SIZE];
__shared__ cuda_real_t s_vel_y[BATCH_SIZE];
__shared__ cuda_real_t s_vel_z[BATCH_SIZE];
__shared__ int s_index[BATCH_SIZE];

s_pos_x[tid] = d_j_pos_x[j + tid];
s_pos_y[tid] = d_j_pos_y[j + tid];
// ...
// Access: s_pos_x[jj]
```

More shared memory registers, but global loads are now coalesced.

## Fields Required on GPU

From CPU `ParticleData` → GPU arrays:

**j-particles (source)**: 8 fields
- `pos_x, pos_y, pos_z` (position)
- `vel_x, vel_y, vel_z` (velocity)
- `mass`
- `particle_index` → `index`

**i-particles (target)**: 8 fields
- `new_pos_x, new_pos_y, new_pos_z` (predicted position)
- `new_vel_x, new_vel_y, new_vel_z` (predicted velocity)
- `neighbor_radius_sq`
- `time_step_reg` → `dt_reg`

**Output**: Written back to CPU
- `acc_irr[dim][order]` (acceleration and jerk, 6 values)
- `new_num_neighbors`
- `neighbors` (via global array)

## Pitfalls to Address

### Pitfall #4: GPU Transfer Fragmentation

**Problem**: 16+ separate `cudaMemcpy` calls vs. original 2.

**Mitigation**:
1. Use CUDA streams for overlapping transfers
2. Group arrays by access pattern
3. Consider pinned memory for host arrays

```cpp
// Use streams to overlap transfers
cudaMemcpyAsync(d_j_pos_x, h_pos_x, n*sizeof(cuda_real_t), H2D, stream);
cudaMemcpyAsync(d_j_pos_y, h_pos_y, n*sizeof(cuda_real_t), H2D, stream);
// Kernel can start as soon as all prereqs arrive
```

### Pitfall #8: Mixed Float/Double Precision

**Problem**: `cuda_real_t` may be float while CPU uses double.

**Mitigation**:
1. Define explicit conversion functions
2. Document precision at each interface
3. Consider `#ifdef CUDA_FLOAT` paths

## Requirement Mapping

| REQ-ID | Requirement | Implementation |
|--------|-------------|----------------|
| GPU-01 | `ParticleDataGPU` container | New struct with device pointers |
| GPU-02 | `copy_from_host()` | Transfer from `ParticleData` to device arrays |
| GPU-03 | `copy_to_host()` | Transfer results back |
| GPU-04 | Update `compute_forces_kernel` | Change to SoA parameters |
| GPU-05 | Update `predict_particle_kernel` | N/A — prediction is CPU-only in current code |

**Note**: GPU-05 may not apply. Current codebase does prediction on CPU, not GPU. The `predict_particle_kernel` mentioned in requirements does not exist in the current code. The GPU only computes forces.

## Implementation Order

1. **Plan 01**: Create `ParticleDataGPU` device container (GPU-01)
   - Struct with device pointers for all GPU-needed fields
   - `allocate()` / `deallocate()` for device memory

2. **Plan 02**: Implement host-device transfers (GPU-02, GPU-03)
   - `copy_from_host(ParticleData&)` for j and i particles
   - `copy_to_host(...)` for results
   - Use CUDA streams for overlap

3. **Plan 03**: Update `compute_forces_kernel` (GPU-04)
   - Change kernel signature to accept separate arrays
   - Update shared memory tiling
   - Update force calculation inner loop

4. **Plan 04**: Update host code in `cuda_acceleration.cu`
   - Replace `GPU` struct with `ParticleDataGPU`
   - Update `_ReceiveFromHost` to use SoA
   - Update `GetAcceleration` result handling

## File Inventory

| File | Changes Needed |
|------|----------------|
| `src/cuda/cuda_defs.h` | Add SoA device pointer struct |
| `src/cuda/cuda_kernels.h` | Update kernel declarations |
| `src/cuda/cuda_kernels.cu` | Update kernel implementations |
| `src/cuda/cuda_acceleration.cu` | Update transfers and host code |
| `src/Makefile` | No changes (files already compiled) |

## Validation

After Phase 3:
- Compile with `make` (no errors)
- Run `workflow/bin/submit.sh --tag soa_gpu_test`
- Verify `dE_over_E0_mean` ≤ 3.2e-5
- Compare `total_wall_s` vs baseline

---
*Generated: 2026-01-17*
