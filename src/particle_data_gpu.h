#ifndef PARTICLE_DATA_GPU_H
#define PARTICLE_DATA_GPU_H

#include "def.h"
#include "cuda/cuda_defs.h"
#include <cuda_runtime.h>

// ============================================================================
// ParticleDataGPU: SoA container for GPU particle data
// ============================================================================
// Holds device pointers for particle fields needed in GPU kernels.
// Use allocate()/deallocate() for device memory management.
// Use copy_from_host()/copy_to_host() for data transfer.
// ============================================================================

struct ParticleDataGPU {
    // ========================================================================
    // Device memory management
    // ========================================================================
    void allocate(size_t j_capacity, size_t i_capacity);
    void deallocate();

    // ========================================================================
    // Host-device transfer
    // ========================================================================
    void copy_j_particles_to_device(
        const double* h_pos_x, const double* h_pos_y, const double* h_pos_z,
        const double* h_vel_x, const double* h_vel_y, const double* h_vel_z,
        const double* h_mass, const int* h_index,
        size_t count, size_t offset, cudaStream_t stream);

    void copy_i_particles_to_device(
        const double* h_pos_x, const double* h_pos_y, const double* h_pos_z,
        const double* h_vel_x, const double* h_vel_y, const double* h_vel_z,
        const double* h_radius_sq, const double* h_dt_reg,
        size_t count, cudaStream_t stream);

    void copy_results_to_host(
        cuda_real_t* h_result, int* h_neighbor, int* h_neighbor_count,
        size_t num_targets, cudaStream_t stream);

    // ========================================================================
    // J-particle device arrays (source particles — exert forces)
    // ========================================================================
    cuda_real_t* d_j_pos_x;
    cuda_real_t* d_j_pos_y;
    cuda_real_t* d_j_pos_z;
    cuda_real_t* d_j_mass;
    cuda_real_t* d_j_vel_x;
    cuda_real_t* d_j_vel_y;
    cuda_real_t* d_j_vel_z;
    int* d_j_index;

    // ========================================================================
    // I-particle device arrays (target particles — receive forces)
    // ========================================================================
    cuda_real_t* d_i_pos_x;
    cuda_real_t* d_i_pos_y;
    cuda_real_t* d_i_pos_z;
    cuda_real_t* d_i_radius_sq;
    cuda_real_t* d_i_vel_x;
    cuda_real_t* d_i_vel_y;
    cuda_real_t* d_i_vel_z;
    cuda_real_t* d_i_dt_reg;

    // ========================================================================
    // Output arrays (on device)
    // ========================================================================
    cuda_real_t* d_result;        // [6 * i_capacity] acc + jerk
    cuda_real_t* d_result_block;  // [6 * GRID_DIM_Y * i_capacity] per-block results
    int* d_neighbor;              // [i_capacity * MAX_NUM_NEIGHBOR]
    int* d_neighbor_block;        // [GRID_DIM_Y * NNB_PER_BLOCK * i_capacity]
    int* d_neighbor_count;        // [i_capacity]
    int* d_neighbor_count_block;  // [GRID_DIM_Y * i_capacity]

    // ========================================================================
    // Pinned host buffers (for async transfer)
    // ========================================================================
    cuda_real_t* h_result;        // [6 * i_capacity]
    int* h_neighbor;              // [i_capacity * MAX_NUM_NEIGHBOR]
    int* h_neighbor_count;        // [i_capacity]

    // ========================================================================
    // Capacity tracking
    // ========================================================================
    size_t j_capacity_;
    size_t i_capacity_;

    // ========================================================================
    // J-particle range for multi-GPU partitioning
    // ========================================================================
    size_t j_start;
    size_t j_count;
};

#endif // PARTICLE_DATA_GPU_H
