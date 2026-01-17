# Plan 03: Update compute_forces Kernel for SoA

## Goal

Update the `compute_forces_kernel` to accept SoA arrays instead of AoS struct arrays.

## Requirement

GPU-04: Update `compute_forces_kernel` for SoA

## Deliverable

Updated `src/cuda/cuda_kernels.cu` and `src/cuda/cuda_kernels.h`

## Tasks

### Task 1: Update kernel declaration in cuda_kernels.h

Change the kernel signature from AoS to SoA parameters.

**Current signature**:
```cpp
__global__ void compute_forces_kernel(
    const i_particle_t* __restrict__ d_i_particles,
    const j_particle_t* __restrict__ d_j_particles,
    cuda_real_t* __restrict__ acc,
    int* __restrict__ neighbor,
    int* num_neighbor,
    int m, int n, int i_start
);
```

**New signature**:
```cpp
__global__ void compute_forces_kernel(
    // I-particle arrays (target)
    const cuda_real_t* __restrict__ d_i_pos_x,
    const cuda_real_t* __restrict__ d_i_pos_y,
    const cuda_real_t* __restrict__ d_i_pos_z,
    const cuda_real_t* __restrict__ d_i_vel_x,
    const cuda_real_t* __restrict__ d_i_vel_y,
    const cuda_real_t* __restrict__ d_i_vel_z,
    const cuda_real_t* __restrict__ d_i_radius_sq,
    // J-particle arrays (source)
    const cuda_real_t* __restrict__ d_j_pos_x,
    const cuda_real_t* __restrict__ d_j_pos_y,
    const cuda_real_t* __restrict__ d_j_pos_z,
    const cuda_real_t* __restrict__ d_j_vel_x,
    const cuda_real_t* __restrict__ d_j_vel_y,
    const cuda_real_t* __restrict__ d_j_vel_z,
    const cuda_real_t* __restrict__ d_j_mass,
    const int* __restrict__ d_j_index,
    // Output arrays
    cuda_real_t* __restrict__ acc,
    int* __restrict__ neighbor,
    int* num_neighbor,
    int m, int n, int i_start
);
```

### Task 2: Update kernel implementation in cuda_kernels.cu

Update the `compute_forces` kernel to use SoA access patterns.

**Key changes**:

1. **Shared memory**: Change from struct array to separate arrays
```cpp
// Old
__shared__ j_particle_t Jp_sh[BATCH_SIZE];

// New
__shared__ cuda_real_t s_j_pos_x[BATCH_SIZE];
__shared__ cuda_real_t s_j_pos_y[BATCH_SIZE];
__shared__ cuda_real_t s_j_pos_z[BATCH_SIZE];
__shared__ cuda_real_t s_j_vel_x[BATCH_SIZE];
__shared__ cuda_real_t s_j_vel_y[BATCH_SIZE];
__shared__ cuda_real_t s_j_vel_z[BATCH_SIZE];
__shared__ cuda_real_t s_j_mass[BATCH_SIZE];
__shared__ int s_j_index[BATCH_SIZE];
```

2. **I-particle load**: Load from separate arrays
```cpp
// Old
i_particle_t Ip = d_Ip[i_ptcl];

// New
cuda_real_t i_pos_x = d_i_pos_x[i_ptcl];
cuda_real_t i_pos_y = d_i_pos_y[i_ptcl];
cuda_real_t i_pos_z = d_i_pos_z[i_ptcl];
cuda_real_t i_vel_x = d_i_vel_x[i_ptcl];
cuda_real_t i_vel_y = d_i_vel_y[i_ptcl];
cuda_real_t i_vel_z = d_i_vel_z[i_ptcl];
cuda_real_t i_radius_sq = d_i_radius_sq[i_ptcl];
```

3. **Shared memory load**: Load separate fields
```cpp
// Old
if (tid < current_batch_size) {
    Jp_sh[tid] = d_Jp[j + tid];
}

// New
if (tid < current_batch_size) {
    int jidx = j + tid;
    s_j_pos_x[tid] = d_j_pos_x[jidx];
    s_j_pos_y[tid] = d_j_pos_y[jidx];
    s_j_pos_z[tid] = d_j_pos_z[jidx];
    s_j_vel_x[tid] = d_j_vel_x[jidx];
    s_j_vel_y[tid] = d_j_vel_y[jidx];
    s_j_vel_z[tid] = d_j_vel_z[jidx];
    s_j_mass[tid] = d_j_mass[jidx];
    s_j_index[tid] = d_j_index[jidx];
}
```

4. **Inner loop access**: Use shared arrays directly
```cpp
// Old
cuda_real_t dx = Jp_sh[jj].pos_x - Ip.pos_x;

// New
cuda_real_t dx = s_j_pos_x[jj] - i_pos_x;
cuda_real_t dy = s_j_pos_y[jj] - i_pos_y;
cuda_real_t dz = s_j_pos_z[jj] - i_pos_z;
// ... similarly for velocity differences and mass
```

### Task 3: Update compute_forces macro

Update the `#define compute_forces` alias to match the new signature.

## Full Updated Kernel

```cpp
__global__ void compute_forces(
    // I-particle arrays (target)
    const cuda_real_t* __restrict__ d_i_pos_x,
    const cuda_real_t* __restrict__ d_i_pos_y,
    const cuda_real_t* __restrict__ d_i_pos_z,
    const cuda_real_t* __restrict__ d_i_vel_x,
    const cuda_real_t* __restrict__ d_i_vel_y,
    const cuda_real_t* __restrict__ d_i_vel_z,
    const cuda_real_t* __restrict__ d_i_radius_sq,
    // J-particle arrays (source)
    const cuda_real_t* __restrict__ d_j_pos_x,
    const cuda_real_t* __restrict__ d_j_pos_y,
    const cuda_real_t* __restrict__ d_j_pos_z,
    const cuda_real_t* __restrict__ d_j_vel_x,
    const cuda_real_t* __restrict__ d_j_vel_y,
    const cuda_real_t* __restrict__ d_j_vel_z,
    const cuda_real_t* __restrict__ d_j_mass,
    const int* __restrict__ d_j_index,
    // Output arrays
    cuda_real_t* __restrict__ acc,
    int* __restrict__ neighbor,
    int* num_neighbor,
    int m, int n, int i_start)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int tid = threadIdx.x;
    int idx_save_size = gridDim.y * m;

    int j_begin = blockIdx.y * n / gridDim.y;
    int j_end = (blockIdx.y + 1) * n / gridDim.y;
    if (blockIdx.y == gridDim.y - 1) j_end = n;

    // SoA shared memory
    __shared__ cuda_real_t s_j_pos_x[BATCH_SIZE];
    __shared__ cuda_real_t s_j_pos_y[BATCH_SIZE];
    __shared__ cuda_real_t s_j_pos_z[BATCH_SIZE];
    __shared__ cuda_real_t s_j_vel_x[BATCH_SIZE];
    __shared__ cuda_real_t s_j_vel_y[BATCH_SIZE];
    __shared__ cuda_real_t s_j_vel_z[BATCH_SIZE];
    __shared__ cuda_real_t s_j_mass[BATCH_SIZE];
    __shared__ int s_j_index[BATCH_SIZE];

    while (i < m + BATCH_SIZE) {
        int i_ptcl = (i < m) ? i + i_start : m - 1 + i_start;

        // Load i-particle data (SoA)
        cuda_real_t i_pos_x = d_i_pos_x[i_ptcl];
        cuda_real_t i_pos_y = d_i_pos_y[i_ptcl];
        cuda_real_t i_pos_z = d_i_pos_z[i_ptcl];
        cuda_real_t i_vel_x = d_i_vel_x[i_ptcl];
        cuda_real_t i_vel_y = d_i_vel_y[i_ptcl];
        cuda_real_t i_vel_z = d_i_vel_z[i_ptcl];
        cuda_real_t i_r2 = d_i_radius_sq[i_ptcl];

        int NumNeighbor = 0;
        int idx_save = i * gridDim.y + blockIdx.y;
        cuda_real_t ax = 0, ay = 0, az = 0;
        cuda_real_t jx = 0, jy = 0, jz = 0;
        int* BlockNeighbor = &neighbor[NNB_PER_BLOCK * idx_save];

        for (int j = j_begin; j < j_end; j += BATCH_SIZE) {
            int current_batch_size = min(BATCH_SIZE, j_end - j);

            __syncthreads();
            if (tid < current_batch_size) {
                int jidx = j + tid;
                s_j_pos_x[tid] = d_j_pos_x[jidx];
                s_j_pos_y[tid] = d_j_pos_y[jidx];
                s_j_pos_z[tid] = d_j_pos_z[jidx];
                s_j_vel_x[tid] = d_j_vel_x[jidx];
                s_j_vel_y[tid] = d_j_vel_y[jidx];
                s_j_vel_z[tid] = d_j_vel_z[jidx];
                s_j_mass[tid] = d_j_mass[jidx];
                s_j_index[tid] = d_j_index[jidx];
            }
            __syncthreads();

            #pragma unroll 4
            for (int jj = 0; jj < current_batch_size; jj++) {
                if (i < m) {
                    cuda_real_t dx = s_j_pos_x[jj] - i_pos_x;
                    cuda_real_t dy = s_j_pos_y[jj] - i_pos_y;
                    cuda_real_t dz = s_j_pos_z[jj] - i_pos_z;
                    cuda_real_t d2 = dx * dx + dy * dy + dz * dz;

                    if (d2 > i_r2) {
                        cuda_real_t dvx = s_j_vel_x[jj] - i_vel_x;
                        cuda_real_t dvy = s_j_vel_y[jj] - i_vel_y;
                        cuda_real_t dvz = s_j_vel_z[jj] - i_vel_z;
                        cuda_real_t inv_sqrt_d2 = rsqrt(d2);
                        cuda_real_t inv_d2 = inv_sqrt_d2 * inv_sqrt_d2;
                        cuda_real_t scale = s_j_mass[jj] * inv_sqrt_d2 * inv_d2;
                        cuda_real_t common_factor = 3.0 * (dx * dvx + dy * dvy + dz * dvz) * inv_d2;

                        ax += scale * dx;
                        ay += scale * dy;
                        az += scale * dz;
                        jx += scale * (dvx - common_factor * dx);
                        jy += scale * (dvy - common_factor * dy);
                        jz += scale * (dvz - common_factor * dz);
                    } else {
                        BlockNeighbor[NumNeighbor++] = s_j_index[jj];
                        assert(NumNeighbor < NNB_PER_BLOCK);
                    }
                }
            }
        }

        if (i < m) {
            acc[idx_save] = ax;
            acc[idx_save + idx_save_size] = ay;
            acc[idx_save + 2 * idx_save_size] = az;
            acc[idx_save + 3 * idx_save_size] = jx;
            acc[idx_save + 4 * idx_save_size] = jy;
            acc[idx_save + 5 * idx_save_size] = jz;
            num_neighbor[idx_save] = NumNeighbor;
        }
        i += gridDim.x * blockDim.x;
    }
}
```

## Verification

- [ ] Kernel compiles without errors
- [ ] Shared memory usage calculated: 8 arrays × BATCH_SIZE × sizeof(type)
- [ ] Global memory loads are coalesced (consecutive threads load consecutive addresses)
- [ ] Inner loop unchanged (same physics)

## Dependencies

- Plan 01, 02 (ParticleDataGPU container)

## Pitfalls Addressed

- Coalesced memory access pattern for better GPU bandwidth utilization

---
*Generated: 2026-01-17*
