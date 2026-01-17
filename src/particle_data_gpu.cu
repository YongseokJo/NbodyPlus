#include "particle_data_gpu.h"
#include "cuda/cuda_defs.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cassert>

// ============================================================================
// Helper: CUDA error checking
// ============================================================================
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                __FILE__, __LINE__, cudaGetErrorString(err)); \
        assert(false); \
    } \
} while(0)

// ============================================================================
// allocate: Allocate device memory for all arrays
// ============================================================================
void ParticleDataGPU::allocate(size_t j_capacity, size_t i_capacity) {
    j_capacity_ = j_capacity;
    i_capacity_ = i_capacity;

    // J-particle arrays
    CUDA_CHECK(cudaMalloc(&d_j_pos_x, j_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_j_pos_y, j_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_j_pos_z, j_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_j_mass, j_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_j_vel_x, j_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_j_vel_y, j_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_j_vel_z, j_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_j_index, j_capacity * sizeof(int)));

    // I-particle arrays
    CUDA_CHECK(cudaMalloc(&d_i_pos_x, i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_i_pos_y, i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_i_pos_z, i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_i_radius_sq, i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_i_vel_x, i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_i_vel_y, i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_i_vel_z, i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_i_dt_reg, i_capacity * sizeof(cuda_real_t)));

    // Output arrays
    CUDA_CHECK(cudaMalloc(&d_result, NUM_FORCE_COMPONENTS * i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_result_block, NUM_FORCE_COMPONENTS * GRID_DIM_Y * i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMalloc(&d_neighbor, i_capacity * MAX_NUM_NEIGHBOR * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_neighbor_block, GRID_DIM_Y * NNB_PER_BLOCK * i_capacity * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_neighbor_count, i_capacity * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_neighbor_count_block, GRID_DIM_Y * i_capacity * sizeof(int)));

    // Pinned host buffers
    CUDA_CHECK(cudaMallocHost(&h_result, NUM_FORCE_COMPONENTS * i_capacity * sizeof(cuda_real_t)));
    CUDA_CHECK(cudaMallocHost(&h_neighbor, i_capacity * MAX_NUM_NEIGHBOR * sizeof(int)));
    CUDA_CHECK(cudaMallocHost(&h_neighbor_count, i_capacity * sizeof(int)));
}

// ============================================================================
// deallocate: Free all device and pinned memory
// ============================================================================
void ParticleDataGPU::deallocate() {
    // J-particle arrays
    cudaFree(d_j_pos_x);
    cudaFree(d_j_pos_y);
    cudaFree(d_j_pos_z);
    cudaFree(d_j_mass);
    cudaFree(d_j_vel_x);
    cudaFree(d_j_vel_y);
    cudaFree(d_j_vel_z);
    cudaFree(d_j_index);

    // I-particle arrays
    cudaFree(d_i_pos_x);
    cudaFree(d_i_pos_y);
    cudaFree(d_i_pos_z);
    cudaFree(d_i_radius_sq);
    cudaFree(d_i_vel_x);
    cudaFree(d_i_vel_y);
    cudaFree(d_i_vel_z);
    cudaFree(d_i_dt_reg);

    // Output arrays
    cudaFree(d_result);
    cudaFree(d_result_block);
    cudaFree(d_neighbor);
    cudaFree(d_neighbor_block);
    cudaFree(d_neighbor_count);
    cudaFree(d_neighbor_count_block);

    // Pinned host buffers
    cudaFreeHost(h_result);
    cudaFreeHost(h_neighbor);
    cudaFreeHost(h_neighbor_count);
}

// ============================================================================
// copy_j_particles_to_device: Transfer j-particles from host (with double→float)
// ============================================================================
void ParticleDataGPU::copy_j_particles_to_device(
    const double* h_pos_x, const double* h_pos_y, const double* h_pos_z,
    const double* h_vel_x, const double* h_vel_y, const double* h_vel_z,
    const double* h_mass, const int* h_index,
    size_t count, size_t offset, cudaStream_t stream)
{
    j_start = offset;
    j_count = count;

    // For double→float conversion, we need temporary pinned buffers
    // or use kernel-based conversion. For simplicity, use direct copy
    // if cuda_real_t is double, or convert on CPU first.

#ifdef CUDA_FLOAT
    // Need conversion — use temporary pinned buffer
    cuda_real_t* tmp;
    CUDA_CHECK(cudaMallocHost(&tmp, count * sizeof(cuda_real_t)));

    // Convert and copy each array
    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_pos_x[offset + i]);
    CUDA_CHECK(cudaMemcpyAsync(d_j_pos_x, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_pos_y[offset + i]);
    CUDA_CHECK(cudaMemcpyAsync(d_j_pos_y, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_pos_z[offset + i]);
    CUDA_CHECK(cudaMemcpyAsync(d_j_pos_z, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_vel_x[offset + i]);
    CUDA_CHECK(cudaMemcpyAsync(d_j_vel_x, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_vel_y[offset + i]);
    CUDA_CHECK(cudaMemcpyAsync(d_j_vel_y, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_vel_z[offset + i]);
    CUDA_CHECK(cudaMemcpyAsync(d_j_vel_z, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_mass[offset + i]);
    CUDA_CHECK(cudaMemcpyAsync(d_j_mass, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    CUDA_CHECK(cudaFreeHost(tmp));
#else
    // Direct copy (double precision GPU)
    CUDA_CHECK(cudaMemcpyAsync(d_j_pos_x, h_pos_x + offset, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_j_pos_y, h_pos_y + offset, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_j_pos_z, h_pos_z + offset, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_j_vel_x, h_vel_x + offset, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_j_vel_y, h_vel_y + offset, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_j_vel_z, h_vel_z + offset, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_j_mass, h_mass + offset, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
#endif

    CUDA_CHECK(cudaMemcpyAsync(d_j_index, h_index + offset, count * sizeof(int), cudaMemcpyHostToDevice, stream));
}

// ============================================================================
// copy_i_particles_to_device: Transfer i-particles from host
// ============================================================================
void ParticleDataGPU::copy_i_particles_to_device(
    const double* h_pos_x, const double* h_pos_y, const double* h_pos_z,
    const double* h_vel_x, const double* h_vel_y, const double* h_vel_z,
    const double* h_radius_sq, const double* h_dt_reg,
    size_t count, cudaStream_t stream)
{
#ifdef CUDA_FLOAT
    cuda_real_t* tmp;
    CUDA_CHECK(cudaMallocHost(&tmp, count * sizeof(cuda_real_t)));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_pos_x[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_pos_x, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_pos_y[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_pos_y, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_pos_z[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_pos_z, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_vel_x[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_vel_x, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_vel_y[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_vel_y, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_vel_z[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_vel_z, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_radius_sq[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_radius_sq, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    for (size_t i = 0; i < count; i++) tmp[i] = static_cast<cuda_real_t>(h_dt_reg[i]);
    CUDA_CHECK(cudaMemcpyAsync(d_i_dt_reg, tmp, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));

    CUDA_CHECK(cudaFreeHost(tmp));
#else
    CUDA_CHECK(cudaMemcpyAsync(d_i_pos_x, h_pos_x, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_i_pos_y, h_pos_y, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_i_pos_z, h_pos_z, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_i_vel_x, h_vel_x, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_i_vel_y, h_vel_y, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_i_vel_z, h_vel_z, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_i_radius_sq, h_radius_sq, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_i_dt_reg, h_dt_reg, count * sizeof(cuda_real_t), cudaMemcpyHostToDevice, stream));
#endif
}

// ============================================================================
// copy_results_to_host: Transfer results from device
// ============================================================================
void ParticleDataGPU::copy_results_to_host(
    cuda_real_t* h_result_out, int* h_neighbor_out, int* h_neighbor_count_out,
    size_t num_targets, cudaStream_t stream)
{
    CUDA_CHECK(cudaMemcpyAsync(h_result_out, d_result,
        NUM_FORCE_COMPONENTS * num_targets * sizeof(cuda_real_t),
        cudaMemcpyDeviceToHost, stream));

    CUDA_CHECK(cudaMemcpyAsync(h_neighbor_out, d_neighbor,
        num_targets * MAX_NUM_NEIGHBOR * sizeof(int),
        cudaMemcpyDeviceToHost, stream));

    CUDA_CHECK(cudaMemcpyAsync(h_neighbor_count_out, d_neighbor_count,
        num_targets * sizeof(int),
        cudaMemcpyDeviceToHost, stream));
}
