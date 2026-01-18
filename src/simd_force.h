#ifndef SIMD_FORCE_H
#define SIMD_FORCE_H

#include <cstddef>

// Forward declarations
struct Particle;
class ParticleData;

// ============================================================================
// SIMD Force Calculation Module
// ============================================================================
// Provides AVX-512 vectorized force calculation for irregular force loops.
// Processes 8 neighbor particles simultaneously using 512-bit vectors.
//
// Usage:
//   NeighborBatch batch;
//   gather_neighbor_data(particles, neighbors, offset, count, dt, batch);
//   double a[3], adot[3];
//   compute_force_vectorized(pos, vel, batch, a, adot);
// ============================================================================

// Maximum number of neighbors in a batch
static constexpr int NEIGHBOR_BATCH_SIZE = 512;

// ============================================================================
// NeighborBatch: Pre-gathered neighbor data in SoA layout
// ============================================================================
// Stores neighbor particle data in contiguous, aligned arrays for efficient
// SIMD processing. All arrays are 64-byte aligned for AVX-512.
struct NeighborBatch {
    // Predicted positions (64-byte aligned for AVX-512)
    alignas(64) double pos_x[NEIGHBOR_BATCH_SIZE];
    alignas(64) double pos_y[NEIGHBOR_BATCH_SIZE];
    alignas(64) double pos_z[NEIGHBOR_BATCH_SIZE];

    // Predicted velocities (64-byte aligned)
    alignas(64) double vel_x[NEIGHBOR_BATCH_SIZE];
    alignas(64) double vel_y[NEIGHBOR_BATCH_SIZE];
    alignas(64) double vel_z[NEIGHBOR_BATCH_SIZE];

    // Particle masses (64-byte aligned)
    alignas(64) double mass[NEIGHBOR_BATCH_SIZE];

    // Particle indices for bookkeeping (tracking which particle is in each slot)
    alignas(64) int indices[NEIGHBOR_BATCH_SIZE];

    // Number of active neighbors in this batch
    int count;

    NeighborBatch() : count(0) {}
};

// ============================================================================
// Pre-gather Functions
// ============================================================================

// Gather active neighbor data into aligned buffers for SIMD processing.
// Skips inactive particles and predicts positions/velocities to the target time.
//
// Parameters:
//   particles      - Global particle array
//   neighbor_indices - Neighbor index array
//   offset         - Offset into neighbor_indices for this particle
//   num_neighbors  - Number of neighbors to process
//   target_time    - Target prediction time
//   batch          - Output: pre-gathered neighbor data
//   cm_indices     - Output: set of CM particle indices encountered (pass empty set)
//   cm_indices_count - Output: number of CM indices added
//
// CM particle handling:
//   When an inactive particle is encountered that has a valid cm_particle_index,
//   that CM index is added to cm_indices for later processing.
void gather_neighbor_data(
    const Particle* particles,
    const int* neighbor_indices,
    int offset,
    int num_neighbors,
    double target_time,
    NeighborBatch& batch,
    int* cm_indices,
    int& cm_indices_count
);

// Gather CM particle data for the CM loop.
// Similar to gather_neighbor_data but works from a list of CM indices.
//
// Parameters:
//   particles      - Global particle array
//   cm_indices     - Array of CM particle indices
//   cm_count       - Number of CM particles
//   self_pid       - PID of the particle being calculated (skip if encountered)
//   target_time    - Target prediction time
//   batch          - Output: pre-gathered CM particle data
void gather_cm_particle_data(
    const Particle* particles,
    const int* cm_indices,
    int cm_count,
    int self_pid,
    double target_time,
    NeighborBatch& batch
);

// ============================================================================
// Force Calculation Functions
// ============================================================================

// Compute gravitational force and jerk from gathered neighbors.
// Uses AVX-512 if available, falls back to scalar otherwise.
//
// Parameters:
//   pos            - Target particle position [3]
//   vel            - Target particle velocity [3]
//   batch          - Pre-gathered neighbor data
//   a_out          - Output: acceleration [3]
//   adot_out       - Output: jerk [3]
void compute_force_vectorized(
    const double pos[3],
    const double vel[3],
    const NeighborBatch& batch,
    double a_out[3],
    double adot_out[3]
);

// ============================================================================
// Platform-specific implementations
// ============================================================================

#ifdef __AVX512F__
// AVX-512 vectorized implementation
// Processes 8 neighbors per iteration using 512-bit vectors.
void compute_force_avx512(
    const double pos[3],
    const double vel[3],
    const NeighborBatch& batch,
    double a_out[3],
    double adot_out[3]
);
#endif

// Scalar fallback implementation
// Used when AVX-512 is not available or for remainder elements.
void compute_force_scalar(
    const double pos[3],
    const double vel[3],
    const NeighborBatch& batch,
    double a_out[3],
    double adot_out[3]
);

#endif // SIMD_FORCE_H
