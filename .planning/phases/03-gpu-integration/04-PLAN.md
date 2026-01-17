# Plan 04: Update cuda_acceleration.cu for SoA

## Goal

Update the host-side GPU dispatch code to use the new SoA `ParticleDataGPU` container and integrate with the CPU `ParticleData` class.

## Requirements

- GPU-02: Host-to-device transfer (integration with `ParticleData`)
- GPU-03: Device-to-host transfer
- GPU-04: Kernel invocation updates

## Deliverable

Updated `src/cuda/cuda_acceleration.cu`

## Tasks

### Task 1: Replace GPU struct with ParticleDataGPU

Replace the current `struct GPU` with a vector of `ParticleDataGPU`.

**Current code**:
```cpp
struct GPU {
    int id = -1;
    cudaStream_t stream = 0;
    j_particle_t* dJ = nullptr;
    i_particle_t* dI = nullptr;
    cuda_real_t* d_result = nullptr;
    // ...
};
static std::vector<GPU> gpu;
```

**New code**:
```cpp
#include "particle_data_gpu.h"

struct GPUContext {
    int id = -1;
    cudaStream_t stream = 0;
    ParticleDataGPU data;
};
static std::vector<GPUContext> gpu;
```

### Task 2: Update _ReceiveFromHost function

Change the function signature to accept `ParticleData` reference instead of AoS vectors.

**Current signature**:
```cpp
void _ReceiveFromHost(
    std::vector<j_particle_t>& hJ,
    std::vector<i_particle_t>& hI
);
```

**New signature**:
```cpp
void _ReceiveFromHost(
    ParticleData& particles,
    const std::vector<int>& j_indices,  // indices of source particles
    const std::vector<int>& i_indices   // indices of target particles
);
```

**Implementation**:
```cpp
void _ReceiveFromHost(
    ParticleData& particles,
    const std::vector<int>& j_indices,
    const std::vector<int>& i_indices)
{
    const size_t nI = i_indices.size();
    const size_t nJ = j_indices.size();
    NNB = nJ;
    isend++;
    assert(NNB <= NBODY_MAX);

    if ((first) || (new_size(NNB) > J_capacity)) {
        J_capacity = new_size(NNB);
        I_capacity = J_capacity;

        if (!first) {
            for (int i = 0; i < deviceCount; i++) {
                cudaSetDevice(i);
                gpu[i].data.deallocate();
            }
        } else {
            first = false;
        }

        for (int i = 0; i < deviceCount; i++) {
            cudaSetDevice(i);
            gpu[i].data.allocate(J_capacity, I_capacity);
        }
    }

    // Partition j-particles across GPUs
    for (int i = 0; i < deviceCount; ++i) {
        cudaSetDevice(i);
        size_t aJ, bJ;
        block_range(nJ, i, deviceCount, aJ, bJ);
        size_t numJ = (bJ > aJ) ? (bJ - aJ) : 0;

        // Copy j-particles using indices
        gpu[i].data.copy_j_particles_to_device(
            particles.new_pos_x(), particles.new_pos_y(), particles.new_pos_z(),
            particles.new_vel_x(), particles.new_vel_y(), particles.new_vel_z(),
            particles.mass(), particles.particle_index(),
            numJ, aJ, gpu[i].stream);

        // Copy i-particles (replicated to all GPUs)
        gpu[i].data.copy_i_particles_to_device(
            particles.new_pos_x(), particles.new_pos_y(), particles.new_pos_z(),
            particles.new_vel_x(), particles.new_vel_y(), particles.new_vel_z(),
            particles.neighbor_radius_sq(), particles.time_step_reg(),
            nI, gpu[i].stream);

        gpu[i].data.j_start = aJ;
        gpu[i].data.j_count = numJ;
    }
}
```

### Task 3: Update GetAcceleration kernel launch

Update the kernel invocation to pass SoA arrays.

**Current code**:
```cpp
compute_forces<<<gridDim, blockDim, 0, gpu[i].stream>>>(
    gpu[i].dI,
    gpu[i].dJ,
    gpu[i].d_result_block,
    gpu[i].d_neighbor_block,
    gpu[i].d_neighbor_count_block,
    NumTarget,
    gpu[i].J_count,
    TargetStart
);
```

**New code**:
```cpp
compute_forces<<<gridDim, blockDim, 0, gpu[i].stream>>>(
    // I-particle arrays
    gpu[i].data.d_i_pos_x,
    gpu[i].data.d_i_pos_y,
    gpu[i].data.d_i_pos_z,
    gpu[i].data.d_i_vel_x,
    gpu[i].data.d_i_vel_y,
    gpu[i].data.d_i_vel_z,
    gpu[i].data.d_i_radius_sq,
    // J-particle arrays
    gpu[i].data.d_j_pos_x,
    gpu[i].data.d_j_pos_y,
    gpu[i].data.d_j_pos_z,
    gpu[i].data.d_j_vel_x,
    gpu[i].data.d_j_vel_y,
    gpu[i].data.d_j_vel_z,
    gpu[i].data.d_j_mass,
    gpu[i].data.d_j_index,
    // Output
    gpu[i].data.d_result_block,
    gpu[i].data.d_neighbor_block,
    gpu[i].data.d_neighbor_count_block,
    NumTarget,
    gpu[i].data.j_count,
    TargetStart
);
```

### Task 4: Update result retrieval

Update `GetAcceleration` to copy results using SoA methods and write to `ParticleData`.

**Current code**:
```cpp
toHost(gpu[i].h_result, gpu[i].d_result, NUM_FORCE_COMPONENTS * NumTarget, gpu[i].stream);

// Later:
for (int k = 0; k < 3; k++) {
    ptcl->acc_irregular[k][0] += static_cast<double>(gpu[i].h_result[j * NUM_FORCE_COMPONENTS + k]);
    ptcl->acc_irregular[k][1] += static_cast<double>(gpu[i].h_result[j * NUM_FORCE_COMPONENTS + k + 3]);
}
```

**New code**:
```cpp
gpu[i].data.copy_results_to_host(
    gpu[i].data.h_result,
    gpu[i].data.h_neighbor,
    gpu[i].data.h_neighbor_count,
    NumTarget,
    gpu[i].stream);

// Later (assuming ParticleData reference passed to GetAcceleration):
for (int dim = 0; dim < 3; dim++) {
    double acc_val = particles.get_acc_irr(idx, dim, 0);
    acc_val += static_cast<double>(gpu[i].data.h_result[j * NUM_FORCE_COMPONENTS + dim]);
    particles.set_acc_irr(idx, dim, 0, acc_val);

    double jrk_val = particles.get_acc_irr(idx, dim, 1);
    jrk_val += static_cast<double>(gpu[i].data.h_result[j * NUM_FORCE_COMPONENTS + dim + 3]);
    particles.set_acc_irr(idx, dim, 1, jrk_val);
}
```

### Task 5: Update extern C interface

Update the external C interface functions.

```cpp
extern "C" {
    void SendToDevice(ParticleData& particles,
                      const std::vector<int>& j_indices,
                      const std::vector<int>& i_indices) {
        _ReceiveFromHost(particles, j_indices, i_indices);
    }

    void CalculateAccelerationOnDevice(int* NumTargetTotal,
                                       std::vector<int>& RegularList,
                                       ParticleData& particles) {
        GetAcceleration(*NumTargetTotal, RegularList, particles);
    }
}
```

### Task 6: Update include directives

Add necessary includes at top of file:

```cpp
#include "particle_data.h"
#include "particle_data_gpu.h"
```

## Verification

- [ ] File compiles with nvcc
- [ ] `_ReceiveFromHost` accepts `ParticleData` reference
- [ ] Kernel launch uses SoA arrays
- [ ] Results written back to `ParticleData` via accessors
- [ ] Multi-GPU partitioning works correctly

## Dependencies

- Plan 01, 02 (ParticleDataGPU container)
- Plan 03 (updated kernel)

## Pitfalls Addressed

- **#4 GPU Transfer Fragmentation**: Uses async transfers with streams
- **#8 Mixed Precision**: Conversion handled in `copy_*` functions

## Notes

The external interface change (`SendToDevice`, `CalculateAccelerationOnDevice`) will require updates in calling code. This may need coordination with Phase 4 (CPU routines) to avoid breaking the build during transition.

**Transition strategy**:
1. Keep old interface temporarily with deprecation warning
2. Add new interface alongside
3. Update callers in Phase 4
4. Remove old interface

---
*Generated: 2026-01-17*
