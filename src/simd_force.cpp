#include "simd_force.h"
#include "particle.h"
#include "global.h"
#include <cmath>
#include <cstring>

#ifdef __AVX512F__
#include <immintrin.h>
#endif

// ============================================================================
// Pre-gather Functions
// ============================================================================

void gather_neighbor_data(
    const Particle* particles,
    const int* neighbor_indices,
    int offset,
    int num_neighbors,
    double target_time,
    NeighborBatch& batch,
    int* cm_indices,
    int& cm_indices_count
) {
    batch.count = 0;
    cm_indices_count = 0;

    for (int i = 0; i < num_neighbors && batch.count < NEIGHBOR_BATCH_SIZE; i++) {
        int neighbor_idx = neighbor_indices[offset + i];
        const Particle* ptcl = &particles[neighbor_idx];

        // Handle inactive particles - check for CM particle
        if (!ptcl->is_active) {
            if (ptcl->cm_particle_index != -1) {
                // Add CM particle index if not already present
                bool found = false;
                for (int j = 0; j < cm_indices_count; j++) {
                    if (cm_indices[j] == ptcl->cm_particle_index) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    cm_indices[cm_indices_count++] = ptcl->cm_particle_index;
                }
            }
            continue;
        }

        // Predict position and velocity to target time
        double dt = target_time - ptcl->current_time_irr;
        double pos[3], vel[3];

        // Inline prediction for efficiency (avoid virtual call overhead)
        double dt_scaled = dt * enzo_time_step;
        if (dt_scaled == 0.0) {
            pos[0] = ptcl->position[0];
            pos[1] = ptcl->position[1];
            pos[2] = ptcl->position[2];
            vel[0] = ptcl->velocity[0];
            vel[1] = ptcl->velocity[1];
            vel[2] = ptcl->velocity[2];
        } else {
            for (int dim = 0; dim < 3; dim++) {
                pos[dim] = ((ptcl->acc_total[dim][1] * dt_scaled / 3.0 + ptcl->acc_total[dim][0])
                           * dt_scaled / 2.0 + ptcl->velocity[dim]) * dt_scaled + ptcl->position[dim];
                vel[dim] = (ptcl->acc_total[dim][1] * dt_scaled / 2.0 + ptcl->acc_total[dim][0])
                          * dt_scaled + ptcl->velocity[dim];
            }
        }

        // Store in batch arrays
        int idx = batch.count;
        batch.pos_x[idx] = pos[0];
        batch.pos_y[idx] = pos[1];
        batch.pos_z[idx] = pos[2];
        batch.vel_x[idx] = vel[0];
        batch.vel_y[idx] = vel[1];
        batch.vel_z[idx] = vel[2];
        batch.mass[idx] = ptcl->mass;
        batch.indices[idx] = neighbor_idx;
        batch.count++;
    }
}

void gather_cm_particle_data(
    const Particle* particles,
    const int* cm_indices,
    int cm_count,
    int self_pid,
    double target_time,
    NeighborBatch& batch
) {
    batch.count = 0;

    for (int i = 0; i < cm_count && batch.count < NEIGHBOR_BATCH_SIZE; i++) {
        int cm_idx = cm_indices[i];
        const Particle* ptcl = &particles[cm_idx];

        // Skip self
        if (ptcl->pid == self_pid) {
            continue;
        }

        // CM particles should always be active
        if (!ptcl->is_active) {
            continue;
        }

        // Predict position and velocity to target time
        double dt = target_time - ptcl->current_time_irr;
        double pos[3], vel[3];

        double dt_scaled = dt * enzo_time_step;
        if (dt_scaled == 0.0) {
            pos[0] = ptcl->position[0];
            pos[1] = ptcl->position[1];
            pos[2] = ptcl->position[2];
            vel[0] = ptcl->velocity[0];
            vel[1] = ptcl->velocity[1];
            vel[2] = ptcl->velocity[2];
        } else {
            for (int dim = 0; dim < 3; dim++) {
                pos[dim] = ((ptcl->acc_total[dim][1] * dt_scaled / 3.0 + ptcl->acc_total[dim][0])
                           * dt_scaled / 2.0 + ptcl->velocity[dim]) * dt_scaled + ptcl->position[dim];
                vel[dim] = (ptcl->acc_total[dim][1] * dt_scaled / 2.0 + ptcl->acc_total[dim][0])
                          * dt_scaled + ptcl->velocity[dim];
            }
        }

        // Store in batch arrays
        int idx = batch.count;
        batch.pos_x[idx] = pos[0];
        batch.pos_y[idx] = pos[1];
        batch.pos_z[idx] = pos[2];
        batch.vel_x[idx] = vel[0];
        batch.vel_y[idx] = vel[1];
        batch.vel_z[idx] = vel[2];
        batch.mass[idx] = ptcl->mass;
        batch.indices[idx] = cm_idx;
        batch.count++;
    }
}

// ============================================================================
// Force Calculation: Vectorized dispatcher
// ============================================================================

void compute_force_vectorized(
    const double pos[3],
    const double vel[3],
    const NeighborBatch& batch,
    double a_out[3],
    double adot_out[3]
) {
#ifdef __AVX512F__
    compute_force_avx512(pos, vel, batch, a_out, adot_out);
#else
    compute_force_scalar(pos, vel, batch, a_out, adot_out);
#endif
}

// ============================================================================
// Force Calculation: Scalar implementation
// ============================================================================

void compute_force_scalar(
    const double pos[3],
    const double vel[3],
    const NeighborBatch& batch,
    double a_out[3],
    double adot_out[3]
) {
    // Initialize outputs
    a_out[0] = a_out[1] = a_out[2] = 0.0;
    adot_out[0] = adot_out[1] = adot_out[2] = 0.0;

    for (int i = 0; i < batch.count; i++) {
        // Compute displacement
        double dx = batch.pos_x[i] - pos[0];
        double dy = batch.pos_y[i] - pos[1];
        double dz = batch.pos_z[i] - pos[2];

        // Compute velocity difference
        double dvx = batch.vel_x[i] - vel[0];
        double dvy = batch.vel_y[i] - vel[1];
        double dvz = batch.vel_z[i] - vel[2];

        // r^2 and v dot r
        double r2 = dx*dx + dy*dy + dz*dz;
        double vdotr = dx*dvx + dy*dvy + dz*dvz;

        // m / r^3
        double r = std::sqrt(r2);
        double r3 = r2 * r;
        double m_r3 = batch.mass[i] / r3;

        // Accumulate acceleration
        a_out[0] += m_r3 * dx;
        a_out[1] += m_r3 * dy;
        a_out[2] += m_r3 * dz;

        // Jerk coefficient: -3 * v.r / r^2
        double coeff = -3.0 * vdotr / r2;

        // Accumulate jerk: m_r3 * (dv + coeff * dx)
        adot_out[0] += m_r3 * (dvx + coeff * dx);
        adot_out[1] += m_r3 * (dvy + coeff * dy);
        adot_out[2] += m_r3 * (dvz + coeff * dz);
    }
}

// ============================================================================
// Force Calculation: AVX-512 implementation
// ============================================================================

#ifdef __AVX512F__

void compute_force_avx512(
    const double pos[3],
    const double vel[3],
    const NeighborBatch& batch,
    double a_out[3],
    double adot_out[3]
) {
    // Initialize accumulators
    __m512d ax = _mm512_setzero_pd();
    __m512d ay = _mm512_setzero_pd();
    __m512d az = _mm512_setzero_pd();
    __m512d adx = _mm512_setzero_pd();
    __m512d ady = _mm512_setzero_pd();
    __m512d adz = _mm512_setzero_pd();

    // Broadcast target position and velocity
    __m512d px = _mm512_set1_pd(pos[0]);
    __m512d py = _mm512_set1_pd(pos[1]);
    __m512d pz = _mm512_set1_pd(pos[2]);
    __m512d vx = _mm512_set1_pd(vel[0]);
    __m512d vy = _mm512_set1_pd(vel[1]);
    __m512d vz = _mm512_set1_pd(vel[2]);

    // Constants
    __m512d three = _mm512_set1_pd(3.0);
    __m512d half = _mm512_set1_pd(0.5);
    __m512d three_halves = _mm512_set1_pd(1.5);

    // Process 8 neighbors at a time
    int i = 0;
    for (; i + 8 <= batch.count; i += 8) {
        // Load 8 neighbor positions
        __m512d nb_px = _mm512_load_pd(&batch.pos_x[i]);
        __m512d nb_py = _mm512_load_pd(&batch.pos_y[i]);
        __m512d nb_pz = _mm512_load_pd(&batch.pos_z[i]);

        // Load 8 neighbor velocities
        __m512d nb_vx = _mm512_load_pd(&batch.vel_x[i]);
        __m512d nb_vy = _mm512_load_pd(&batch.vel_y[i]);
        __m512d nb_vz = _mm512_load_pd(&batch.vel_z[i]);

        // Load 8 masses
        __m512d nb_m = _mm512_load_pd(&batch.mass[i]);

        // Compute displacement: dx = nb_pos - pos
        __m512d dx = _mm512_sub_pd(nb_px, px);
        __m512d dy = _mm512_sub_pd(nb_py, py);
        __m512d dz = _mm512_sub_pd(nb_pz, pz);

        // Compute velocity difference: dv = nb_vel - vel
        __m512d dvx = _mm512_sub_pd(nb_vx, vx);
        __m512d dvy = _mm512_sub_pd(nb_vy, vy);
        __m512d dvz = _mm512_sub_pd(nb_vz, vz);

        // r^2 = dx^2 + dy^2 + dz^2
        __m512d r2 = _mm512_fmadd_pd(dx, dx, _mm512_fmadd_pd(dy, dy, _mm512_mul_pd(dz, dz)));

        // v dot r = dx*dvx + dy*dvy + dz*dvz
        __m512d vdotr = _mm512_fmadd_pd(dx, dvx, _mm512_fmadd_pd(dy, dvy, _mm512_mul_pd(dz, dvz)));

        // Compute 1/sqrt(r2) using fast rsqrt with Newton-Raphson refinement
        // rsqrt14 gives ~14 bits of precision, NR gives ~28 bits (enough for double)
        __m512d r_inv = _mm512_rsqrt14_pd(r2);

        // Newton-Raphson: r_inv = r_inv * (1.5 - 0.5 * r2 * r_inv^2)
        __m512d r_inv_sq = _mm512_mul_pd(r_inv, r_inv);
        r_inv = _mm512_mul_pd(r_inv, _mm512_fnmadd_pd(half, _mm512_mul_pd(r2, r_inv_sq), three_halves));

        // Second Newton-Raphson iteration for better precision
        r_inv_sq = _mm512_mul_pd(r_inv, r_inv);
        r_inv = _mm512_mul_pd(r_inv, _mm512_fnmadd_pd(half, _mm512_mul_pd(r2, r_inv_sq), three_halves));

        // r^(-3) = r_inv^3
        __m512d r3_inv = _mm512_mul_pd(r_inv, _mm512_mul_pd(r_inv, r_inv));

        // m / r^3
        __m512d m_r3 = _mm512_mul_pd(nb_m, r3_inv);

        // Accumulate acceleration: a += m_r3 * dx
        ax = _mm512_fmadd_pd(m_r3, dx, ax);
        ay = _mm512_fmadd_pd(m_r3, dy, ay);
        az = _mm512_fmadd_pd(m_r3, dz, az);

        // Jerk coefficient: coeff = -3 * vdotr / r2
        // Note: using positive 3 and subtracting later for better FMA usage
        __m512d coeff = _mm512_mul_pd(three, _mm512_div_pd(vdotr, r2));

        // Jerk: adot += m_r3 * (dv - coeff * dx)
        // Using fnmadd: result = -(a*b) + c = c - a*b
        adx = _mm512_fmadd_pd(m_r3, _mm512_fnmadd_pd(coeff, dx, dvx), adx);
        ady = _mm512_fmadd_pd(m_r3, _mm512_fnmadd_pd(coeff, dy, dvy), ady);
        adz = _mm512_fmadd_pd(m_r3, _mm512_fnmadd_pd(coeff, dz, dvz), adz);
    }

    // Reduce vector accumulators to scalars
    a_out[0] = _mm512_reduce_add_pd(ax);
    a_out[1] = _mm512_reduce_add_pd(ay);
    a_out[2] = _mm512_reduce_add_pd(az);
    adot_out[0] = _mm512_reduce_add_pd(adx);
    adot_out[1] = _mm512_reduce_add_pd(ady);
    adot_out[2] = _mm512_reduce_add_pd(adz);

    // Handle remainder with scalar code
    for (; i < batch.count; i++) {
        double dx = batch.pos_x[i] - pos[0];
        double dy = batch.pos_y[i] - pos[1];
        double dz = batch.pos_z[i] - pos[2];

        double dvx = batch.vel_x[i] - vel[0];
        double dvy = batch.vel_y[i] - vel[1];
        double dvz = batch.vel_z[i] - vel[2];

        double r2 = dx*dx + dy*dy + dz*dz;
        double vdotr = dx*dvx + dy*dvy + dz*dvz;

        double r = std::sqrt(r2);
        double r3 = r2 * r;
        double m_r3 = batch.mass[i] / r3;

        a_out[0] += m_r3 * dx;
        a_out[1] += m_r3 * dy;
        a_out[2] += m_r3 * dz;

        double coeff = -3.0 * vdotr / r2;

        adot_out[0] += m_r3 * (dvx + coeff * dx);
        adot_out[1] += m_r3 * (dvy + coeff * dy);
        adot_out[2] += m_r3 * (dvz + coeff * dz);
    }
}

#endif // __AVX512F__
