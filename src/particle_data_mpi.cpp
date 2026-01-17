#include "particle_data_mpi.h"
#include <cstddef>

// ============================================================================
// Helper: Compute padded capacity to avoid cache associativity issues
// ============================================================================
namespace {
size_t padded_capacity(size_t n) {
    // Avoid exact powers of 2 (cache associativity issues)
    if (n > 64 && (n & (n - 1)) == 0) {
        return n + 64;
    }
    // Align to 8 elements (64 bytes for doubles)
    return ((n + 7) / 8) * 8;
}
}

// ============================================================================
// Constructor: Initialize all MPI state
// ============================================================================
ParticleDataMPI::ParticleDataMPI()
    : ParticleData(),
      shared_comm_(MPI_COMM_NULL),
      is_shared_(false),
      win_mass_(MPI_WIN_NULL),
      win_neighbor_radius_sq_(MPI_WIN_NULL)
{
    // Initialize position/velocity window arrays
    for (int i = 0; i < 3; i++) {
        win_pos_[i] = MPI_WIN_NULL;
        win_vel_[i] = MPI_WIN_NULL;
        win_new_pos_[i] = MPI_WIN_NULL;
        win_new_vel_[i] = MPI_WIN_NULL;
    }

    // Initialize acceleration window arrays
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            win_acc_total_[d][o] = MPI_WIN_NULL;
            win_acc_reg_[d][o] = MPI_WIN_NULL;
            win_acc_irr_[d][o] = MPI_WIN_NULL;
        }
    }

    // Initialize timestep windows
    for (int i = 0; i < 4; i++) {
        win_timestep_[i] = MPI_WIN_NULL;
    }

    // Initialize block time windows
    for (int i = 0; i < 6; i++) {
        win_block_[i] = MPI_WIN_NULL;
    }

    // Initialize neighbor int windows
    for (int i = 0; i < 3; i++) {
        win_neighbor_int_[i] = MPI_WIN_NULL;
    }

    // Initialize ID/index windows
    for (int i = 0; i < 3; i++) {
        win_id_[i] = MPI_WIN_NULL;
    }

    // Initialize time level windows
    for (int i = 0; i < 2; i++) {
        win_time_level_[i] = MPI_WIN_NULL;
    }

    // Initialize bool flag windows
    for (int i = 0; i < 3; i++) {
        win_flag_[i] = MPI_WIN_NULL;
    }

    // Initialize other double windows
    for (int i = 0; i < 2; i++) {
        win_other_double_[i] = MPI_WIN_NULL;
    }
}

// ============================================================================
// Destructor: Free MPI shared memory if allocated
// ============================================================================
ParticleDataMPI::~ParticleDataMPI() {
    if (is_shared_) {
        deallocate_shared();
    }
}

// ============================================================================
// Private helpers: Type-specific MPI shared memory allocation
// ============================================================================

void ParticleDataMPI::allocate_double_array_shared(double** ptr, MPI_Win* win, size_t capacity) {
    int shared_rank;
    MPI_Comm_rank(shared_comm_, &shared_rank);

    if (shared_rank == 0) {
        MPI_Win_allocate_shared(sizeof(double) * capacity, sizeof(double),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    } else {
        MPI_Win_allocate_shared(0, sizeof(double),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    }

    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(*win, 0, &size, &disp_unit, ptr);
}

void ParticleDataMPI::allocate_ull_array_shared(ull_t** ptr, MPI_Win* win, size_t capacity) {
    int shared_rank;
    MPI_Comm_rank(shared_comm_, &shared_rank);

    if (shared_rank == 0) {
        MPI_Win_allocate_shared(sizeof(ull_t) * capacity, sizeof(ull_t),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    } else {
        MPI_Win_allocate_shared(0, sizeof(ull_t),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    }

    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(*win, 0, &size, &disp_unit, ptr);
}

void ParticleDataMPI::allocate_int_array_shared(int** ptr, MPI_Win* win, size_t capacity) {
    int shared_rank;
    MPI_Comm_rank(shared_comm_, &shared_rank);

    if (shared_rank == 0) {
        MPI_Win_allocate_shared(sizeof(int) * capacity, sizeof(int),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    } else {
        MPI_Win_allocate_shared(0, sizeof(int),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    }

    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(*win, 0, &size, &disp_unit, ptr);
}

void ParticleDataMPI::allocate_bool_array_shared(bool** ptr, MPI_Win* win, size_t capacity) {
    int shared_rank;
    MPI_Comm_rank(shared_comm_, &shared_rank);

    if (shared_rank == 0) {
        MPI_Win_allocate_shared(sizeof(bool) * capacity, sizeof(bool),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    } else {
        MPI_Win_allocate_shared(0, sizeof(bool),
                                MPI_INFO_NULL, shared_comm_, ptr, win);
    }

    MPI_Aint size;
    int disp_unit;
    MPI_Win_shared_query(*win, 0, &size, &disp_unit, ptr);
}

// ============================================================================
// allocate_shared: Allocate all 66 arrays via MPI shared memory
// ============================================================================
void ParticleDataMPI::allocate_shared(size_t requested_capacity, MPI_Comm comm) {
    // Store communicator
    shared_comm_ = comm;

    // Apply padding to avoid cache associativity issues
    capacity_ = padded_capacity(requested_capacity);
    count_ = 0;

    // ========================================================================
    // 1. Position arrays (3 windows)
    // ========================================================================
    allocate_double_array_shared(&pos_x_, &win_pos_[0], capacity_);
    allocate_double_array_shared(&pos_y_, &win_pos_[1], capacity_);
    allocate_double_array_shared(&pos_z_, &win_pos_[2], capacity_);

    // ========================================================================
    // 2. Velocity arrays (3 windows)
    // ========================================================================
    allocate_double_array_shared(&vel_x_, &win_vel_[0], capacity_);
    allocate_double_array_shared(&vel_y_, &win_vel_[1], capacity_);
    allocate_double_array_shared(&vel_z_, &win_vel_[2], capacity_);

    // ========================================================================
    // 3. Mass (1 window)
    // ========================================================================
    allocate_double_array_shared(&mass_, &win_mass_, capacity_);

    // ========================================================================
    // 4. Acceleration arrays (36 windows: 3 types × 3 dims × 4 orders)
    // ========================================================================
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            allocate_double_array_shared(&acc_total_[d][o], &win_acc_total_[d][o], capacity_);
            allocate_double_array_shared(&acc_reg_[d][o], &win_acc_reg_[d][o], capacity_);
            allocate_double_array_shared(&acc_irr_[d][o], &win_acc_irr_[d][o], capacity_);
        }
    }

    // ========================================================================
    // 5. New position/velocity (6 windows)
    // ========================================================================
    allocate_double_array_shared(&new_pos_x_, &win_new_pos_[0], capacity_);
    allocate_double_array_shared(&new_pos_y_, &win_new_pos_[1], capacity_);
    allocate_double_array_shared(&new_pos_z_, &win_new_pos_[2], capacity_);
    allocate_double_array_shared(&new_vel_x_, &win_new_vel_[0], capacity_);
    allocate_double_array_shared(&new_vel_y_, &win_new_vel_[1], capacity_);
    allocate_double_array_shared(&new_vel_z_, &win_new_vel_[2], capacity_);

    // ========================================================================
    // 6. Neighbor radius (1 window)
    // ========================================================================
    allocate_double_array_shared(&neighbor_radius_sq_, &win_neighbor_radius_sq_, capacity_);

    // ========================================================================
    // 7. Timestep doubles (4 windows)
    // ========================================================================
    allocate_double_array_shared(&current_time_irr_, &win_timestep_[0], capacity_);
    allocate_double_array_shared(&current_time_reg_, &win_timestep_[1], capacity_);
    allocate_double_array_shared(&time_step_irr_, &win_timestep_[2], capacity_);
    allocate_double_array_shared(&time_step_reg_, &win_timestep_[3], capacity_);

    // ========================================================================
    // 8. Block times - ull_t (6 windows)
    // ========================================================================
    allocate_ull_array_shared(&current_block_irr_, &win_block_[0], capacity_);
    allocate_ull_array_shared(&new_current_block_irr_, &win_block_[1], capacity_);
    allocate_ull_array_shared(&current_block_reg_, &win_block_[2], capacity_);
    allocate_ull_array_shared(&time_block_irr_, &win_block_[3], capacity_);
    allocate_ull_array_shared(&time_block_reg_, &win_block_[4], capacity_);
    allocate_ull_array_shared(&next_block_irr_, &win_block_[5], capacity_);

    // ========================================================================
    // 9. Neighbor ints (3 windows)
    // ========================================================================
    allocate_int_array_shared(&num_neighbors_, &win_neighbor_int_[0], capacity_);
    allocate_int_array_shared(&new_num_neighbors_, &win_neighbor_int_[1], capacity_);
    allocate_int_array_shared(&neighbors_offset_, &win_neighbor_int_[2], capacity_);

    // ========================================================================
    // 10. ID/index ints (3 windows)
    // ========================================================================
    allocate_int_array_shared(&pid_, &win_id_[0], capacity_);
    allocate_int_array_shared(&particle_index_, &win_id_[1], capacity_);
    allocate_int_array_shared(&particle_type_, &win_id_[2], capacity_);

    // ========================================================================
    // 11. Time level ints (2 windows)
    // ========================================================================
    allocate_int_array_shared(&time_level_irr_, &win_time_level_[0], capacity_);
    allocate_int_array_shared(&time_level_reg_, &win_time_level_[1], capacity_);

    // ========================================================================
    // 12. Bool flags (3 windows)
    // ========================================================================
    allocate_bool_array_shared(&is_active_, &win_flag_[0], capacity_);
    allocate_bool_array_shared(&is_up_to_date_, &win_flag_[1], capacity_);
    allocate_bool_array_shared(&is_cm_particle_, &win_flag_[2], capacity_);

    // ========================================================================
    // 13. Other doubles (2 windows)
    // ========================================================================
    allocate_double_array_shared(&radius_, &win_other_double_[0], capacity_);
    allocate_double_array_shared(&delta_mass_, &win_other_double_[1], capacity_);

    is_shared_ = true;
}

// ============================================================================
// deallocate_shared: Free all MPI windows
// ============================================================================
void ParticleDataMPI::deallocate_shared() {
    if (!is_shared_) return;

    // ========================================================================
    // Position windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        if (win_pos_[i] != MPI_WIN_NULL) MPI_Win_free(&win_pos_[i]);
    }

    // ========================================================================
    // Velocity windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        if (win_vel_[i] != MPI_WIN_NULL) MPI_Win_free(&win_vel_[i]);
    }

    // ========================================================================
    // Mass window (1)
    // ========================================================================
    if (win_mass_ != MPI_WIN_NULL) MPI_Win_free(&win_mass_);

    // ========================================================================
    // Acceleration windows (36)
    // ========================================================================
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            if (win_acc_total_[d][o] != MPI_WIN_NULL) MPI_Win_free(&win_acc_total_[d][o]);
            if (win_acc_reg_[d][o] != MPI_WIN_NULL) MPI_Win_free(&win_acc_reg_[d][o]);
            if (win_acc_irr_[d][o] != MPI_WIN_NULL) MPI_Win_free(&win_acc_irr_[d][o]);
        }
    }

    // ========================================================================
    // New position/velocity windows (6)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        if (win_new_pos_[i] != MPI_WIN_NULL) MPI_Win_free(&win_new_pos_[i]);
        if (win_new_vel_[i] != MPI_WIN_NULL) MPI_Win_free(&win_new_vel_[i]);
    }

    // ========================================================================
    // Neighbor radius window (1)
    // ========================================================================
    if (win_neighbor_radius_sq_ != MPI_WIN_NULL) MPI_Win_free(&win_neighbor_radius_sq_);

    // ========================================================================
    // Timestep windows (4)
    // ========================================================================
    for (int i = 0; i < 4; i++) {
        if (win_timestep_[i] != MPI_WIN_NULL) MPI_Win_free(&win_timestep_[i]);
    }

    // ========================================================================
    // Block time windows (6)
    // ========================================================================
    for (int i = 0; i < 6; i++) {
        if (win_block_[i] != MPI_WIN_NULL) MPI_Win_free(&win_block_[i]);
    }

    // ========================================================================
    // Neighbor int windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        if (win_neighbor_int_[i] != MPI_WIN_NULL) MPI_Win_free(&win_neighbor_int_[i]);
    }

    // ========================================================================
    // ID/index windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        if (win_id_[i] != MPI_WIN_NULL) MPI_Win_free(&win_id_[i]);
    }

    // ========================================================================
    // Time level windows (2)
    // ========================================================================
    for (int i = 0; i < 2; i++) {
        if (win_time_level_[i] != MPI_WIN_NULL) MPI_Win_free(&win_time_level_[i]);
    }

    // ========================================================================
    // Bool flag windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        if (win_flag_[i] != MPI_WIN_NULL) MPI_Win_free(&win_flag_[i]);
    }

    // ========================================================================
    // Other double windows (2)
    // ========================================================================
    for (int i = 0; i < 2; i++) {
        if (win_other_double_[i] != MPI_WIN_NULL) MPI_Win_free(&win_other_double_[i]);
    }

    // Reset all array pointers to nullptr
    pos_x_ = pos_y_ = pos_z_ = nullptr;
    vel_x_ = vel_y_ = vel_z_ = nullptr;
    mass_ = nullptr;
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            acc_total_[d][o] = nullptr;
            acc_reg_[d][o] = nullptr;
            acc_irr_[d][o] = nullptr;
        }
    }
    new_pos_x_ = new_pos_y_ = new_pos_z_ = nullptr;
    new_vel_x_ = new_vel_y_ = new_vel_z_ = nullptr;
    neighbor_radius_sq_ = nullptr;
    current_time_irr_ = current_time_reg_ = nullptr;
    time_step_irr_ = time_step_reg_ = nullptr;
    current_block_irr_ = new_current_block_irr_ = nullptr;
    current_block_reg_ = time_block_irr_ = time_block_reg_ = next_block_irr_ = nullptr;
    num_neighbors_ = new_num_neighbors_ = neighbors_offset_ = nullptr;
    pid_ = particle_index_ = particle_type_ = nullptr;
    time_level_irr_ = time_level_reg_ = nullptr;
    is_active_ = is_up_to_date_ = is_cm_particle_ = nullptr;
    radius_ = delta_mass_ = nullptr;

    shared_comm_ = MPI_COMM_NULL;
    is_shared_ = false;
    capacity_ = 0;
    count_ = 0;
}

// ============================================================================
// sync_all: Synchronize all MPI windows
// ============================================================================
void ParticleDataMPI::sync_all() {
    if (!is_shared_) return;

    // ========================================================================
    // Position windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        MPI_Win_sync(win_pos_[i]);
    }

    // ========================================================================
    // Velocity windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        MPI_Win_sync(win_vel_[i]);
    }

    // ========================================================================
    // Mass window (1)
    // ========================================================================
    MPI_Win_sync(win_mass_);

    // ========================================================================
    // Acceleration windows (36)
    // ========================================================================
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            MPI_Win_sync(win_acc_total_[d][o]);
            MPI_Win_sync(win_acc_reg_[d][o]);
            MPI_Win_sync(win_acc_irr_[d][o]);
        }
    }

    // ========================================================================
    // New position/velocity windows (6)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        MPI_Win_sync(win_new_pos_[i]);
        MPI_Win_sync(win_new_vel_[i]);
    }

    // ========================================================================
    // Neighbor radius window (1)
    // ========================================================================
    MPI_Win_sync(win_neighbor_radius_sq_);

    // ========================================================================
    // Timestep windows (4)
    // ========================================================================
    for (int i = 0; i < 4; i++) {
        MPI_Win_sync(win_timestep_[i]);
    }

    // ========================================================================
    // Block time windows (6)
    // ========================================================================
    for (int i = 0; i < 6; i++) {
        MPI_Win_sync(win_block_[i]);
    }

    // ========================================================================
    // Neighbor int windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        MPI_Win_sync(win_neighbor_int_[i]);
    }

    // ========================================================================
    // ID/index windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        MPI_Win_sync(win_id_[i]);
    }

    // ========================================================================
    // Time level windows (2)
    // ========================================================================
    for (int i = 0; i < 2; i++) {
        MPI_Win_sync(win_time_level_[i]);
    }

    // ========================================================================
    // Bool flag windows (3)
    // ========================================================================
    for (int i = 0; i < 3; i++) {
        MPI_Win_sync(win_flag_[i]);
    }

    // ========================================================================
    // Other double windows (2)
    // ========================================================================
    for (int i = 0; i < 2; i++) {
        MPI_Win_sync(win_other_double_[i]);
    }

    // Memory barrier
    MPI_Barrier(shared_comm_);
}
