#ifndef PARTICLE_DATA_MPI_H
#define PARTICLE_DATA_MPI_H

#include "particle_data.h"
#include <mpi.h>

// ============================================================================
// ParticleDataMPI: MPI shared memory extension of ParticleData
// ============================================================================
// This class extends ParticleData to support MPI shared memory allocation.
// Each SoA array gets its own MPI_Win for fine-grained memory sharing.
//
// Usage:
//   ParticleDataMPI data;
//   data.allocate_shared(1000, shared_comm);  // Allocate via MPI shared memory
//   // ... use data normally via inherited accessors ...
//   data.sync_all();  // Synchronize all windows
//   data.deallocate_shared();  // Free MPI windows
// ============================================================================

class ParticleDataMPI : public ParticleData {
public:
    // ========================================================================
    // Lifecycle
    // ========================================================================
    ParticleDataMPI();
    ~ParticleDataMPI();

    // ========================================================================
    // MPI shared memory allocation
    // ========================================================================
    void allocate_shared(size_t capacity, MPI_Comm comm);
    void deallocate_shared();
    void sync_all();
    bool is_shared() const { return is_shared_; }

private:
    // ========================================================================
    // Private helper methods for type-specific allocation
    // ========================================================================
    void allocate_double_array_shared(double** ptr, MPI_Win* win, size_t capacity);
    void allocate_ull_array_shared(ull_t** ptr, MPI_Win* win, size_t capacity);
    void allocate_int_array_shared(int** ptr, MPI_Win* win, size_t capacity);
    void allocate_bool_array_shared(bool** ptr, MPI_Win* win, size_t capacity);

    // ========================================================================
    // MPI state
    // ========================================================================
    MPI_Comm shared_comm_;
    bool is_shared_;

    // ========================================================================
    // MPI window handles for all SoA arrays (66 windows total)
    // ========================================================================

    // Position windows (3)
    MPI_Win win_pos_[3];

    // Velocity windows (3)
    MPI_Win win_vel_[3];

    // Mass window (1)
    MPI_Win win_mass_;

    // Acceleration windows (36 total: 3 types × 3 dims × 4 orders)
    MPI_Win win_acc_total_[3][4];
    MPI_Win win_acc_reg_[3][4];
    MPI_Win win_acc_irr_[3][4];

    // New position/velocity windows (6)
    MPI_Win win_new_pos_[3];
    MPI_Win win_new_vel_[3];

    // Neighbor radius squared window (1)
    MPI_Win win_neighbor_radius_sq_;

    // Timestep double windows (4)
    MPI_Win win_timestep_[4];  // current_time_irr, current_time_reg, time_step_irr, time_step_reg

    // Block time ull_t windows (6)
    MPI_Win win_block_[6];

    // Neighbor int windows (3)
    MPI_Win win_neighbor_int_[3];

    // ID/index int windows (3)
    MPI_Win win_id_[3];

    // Time level int windows (2)
    MPI_Win win_time_level_[2];

    // Bool flag windows (3)
    MPI_Win win_flag_[3];

    // Other double windows (2)
    MPI_Win win_other_double_[2];  // radius, delta_mass
};

#endif // PARTICLE_DATA_MPI_H
