#ifndef CUDA_DEFS_H
#define CUDA_DEFS_H

#include "../def.h"

// ============================================================================
// CUDA constants
// ============================================================================
#define _PROFILE
#define NUM_FORCE_COMPONENTS 6      // acc(3) + jerk(3)
#define NUM_POS_COMPONENTS 2
#define NUM_RESULT_COMPONENTS 7     // 6 force + 1 extra

// Round up to next power of 2 (minimum 1024)
#define new_size(A) ((A > 1024) ? int(pow(2, ceil(log(A) / log(2.0)))) : 1024)

// ============================================================================
// GPU particle structures for force calculation
// ============================================================================

// Target particle (i-particle) - receives forces
struct i_particle_t {
    cuda_real_t pos_x;      // Position x
    cuda_real_t pos_y;      // Position y
    cuda_real_t pos_z;      // Position z
    cuda_real_t radius_sq;  // Neighbor search radius squared

    cuda_real_t vel_x;      // Velocity x
    cuda_real_t vel_y;      // Velocity y
    cuda_real_t vel_z;      // Velocity z
    cuda_real_t dt_reg;     // Regular time step
};

// Source particle (j-particle) - exerts forces
struct j_particle_t {
    cuda_real_t pos_x;      // Position x
    cuda_real_t pos_y;      // Position y
    cuda_real_t pos_z;      // Position z
    cuda_real_t mass;       // Particle mass

    cuda_real_t vel_x;      // Velocity x
    cuda_real_t vel_y;      // Velocity y
    cuda_real_t vel_z;      // Velocity z
    int index;              // Particle index
#ifndef CUDA_FLOAT
    int pad;                // Padding for alignment (double precision)
#endif
};

#endif
