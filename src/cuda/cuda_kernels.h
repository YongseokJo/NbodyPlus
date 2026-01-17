#ifndef CUDA_KERNELS_H
#define CUDA_KERNELS_H

#include "../def.h"
#include "cuda_defs.h"

// ============================================================================
// CUDA kernel declarations
// ============================================================================

// Initialization kernel
__global__ void initialize_arrays(
    cuda_real_t* result,
    cuda_real_t* diff,
    int n,
    int m,
    int* subset
);

// Legacy alias
#define initialize initialize_arrays

// Pairwise difference computation
__global__ void compute_pairwise_diff_subset(
    const cuda_real_t* ptcl,
    cuda_real_t* diff,
    int n,
    int m,
    const int* subset,
    int start
);

// Magnitude computation with neighbor detection
__global__ void compute_magnitudes_subset(
    const cuda_real_t* r2,
    const cuda_real_t* diff,
    cuda_real_t* magnitudes,
    int n,
    int m,
    int* subset,
    bool* neighbor2,
    int start
);

// Force computation (subset version)
__global__ void compute_forces_subset(
    const cuda_real_t* ptcl,
    cuda_real_t* diff,
    const cuda_real_t* magnitudes,
    int n,
    int m,
    const int* subset
);

// Neighbor assignment
__global__ void assign_neighbor(
    int* neighbor,
    int* num_neighbor,
    const cuda_real_t* r2,
    const cuda_real_t* magnitudes,
    int n,
    int m,
    const int* subset
);

// Force reduction
__global__ void reduce_forces(
    const cuda_real_t* diff,
    cuda_real_t* result,
    int n,
    int m
);

// Debug print
__global__ void print_forces_subset(
    cuda_real_t* result,
    int m,
    int n
);

// ============================================================================
// Main force computation kernel (SoA version)
// ============================================================================
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
    int m,
    int n,
    int i_start
);

// Legacy alias
#define compute_forces compute_forces_kernel

// ============================================================================
// Neighbor gathering kernels
// ============================================================================
__global__ void gather_neighbor_kernel(
    const int* neighbor_block,
    const int* num_neighbor,
    int* gathered_neighbor,
    int m
);

#define gather_neighbor gather_neighbor_kernel

__global__ void gather_num_neighbor_kernel(
    const int* numneighbor_block,
    int* gathered_numneighbor,
    int m
);

#define gather_numneighbor gather_num_neighbor_kernel

// ============================================================================
// Force reduction kernel (optimized with warp shuffle)
// ============================================================================
__global__ void reduce_forces_kernel(
    const cuda_real_t* diff,    // [6 * m * n] input
    cuda_real_t* result,        // [6 * m] output
    int n,                      // rows in each component
    int m                       // columns
);

#endif
