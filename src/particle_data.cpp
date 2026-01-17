#include "particle_data.h"

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
// Constructor: Initialize all pointers to nullptr
// ============================================================================
ParticleData::ParticleData()
    : capacity_(0), count_(0),
      // Critical fields
      pos_x_(nullptr), pos_y_(nullptr), pos_z_(nullptr),
      vel_x_(nullptr), vel_y_(nullptr), vel_z_(nullptr),
      mass_(nullptr),
      new_pos_x_(nullptr), new_pos_y_(nullptr), new_pos_z_(nullptr),
      new_vel_x_(nullptr), new_vel_y_(nullptr), new_vel_z_(nullptr),
      neighbor_radius_sq_(nullptr),
      // Important fields - timestep doubles
      current_time_irr_(nullptr), current_time_reg_(nullptr),
      time_step_irr_(nullptr), time_step_reg_(nullptr),
      // Important fields - block times
      current_block_irr_(nullptr), new_current_block_irr_(nullptr),
      current_block_reg_(nullptr), time_block_irr_(nullptr),
      time_block_reg_(nullptr), next_block_irr_(nullptr),
      // Important fields - neighbor ints
      num_neighbors_(nullptr), new_num_neighbors_(nullptr), neighbors_offset_(nullptr),
      // Low priority fields - IDs/indices
      pid_(nullptr), particle_index_(nullptr), particle_type_(nullptr),
      // Low priority fields - time levels
      time_level_irr_(nullptr), time_level_reg_(nullptr),
      // Low priority fields - bools
      is_active_(nullptr), is_up_to_date_(nullptr), is_cm_particle_(nullptr),
      // Low priority fields - other doubles
      radius_(nullptr), delta_mass_(nullptr)
{
    // Initialize acceleration pointer arrays to nullptr
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            acc_total_[d][o] = nullptr;
            acc_reg_[d][o] = nullptr;
            acc_irr_[d][o] = nullptr;
        }
    }
}

// ============================================================================
// Destructor: Free all allocated memory
// ============================================================================
ParticleData::~ParticleData() {
    deallocate();
}

// ============================================================================
// Allocate: Allocate all 66 arrays with padded capacity
// ============================================================================
void ParticleData::allocate(size_t requested_capacity) {
    // Clean up existing allocation if any
    if (capacity_ > 0) {
        deallocate();
    }

    // Apply padding to avoid cache associativity issues
    capacity_ = padded_capacity(requested_capacity);
    count_ = 0;

    // ========================================================================
    // Critical double arrays (43 total)
    // ========================================================================
    // Position (3)
    pos_x_ = new double[capacity_]();
    pos_y_ = new double[capacity_]();
    pos_z_ = new double[capacity_]();

    // Velocity (3)
    vel_x_ = new double[capacity_]();
    vel_y_ = new double[capacity_]();
    vel_z_ = new double[capacity_]();

    // Mass (1)
    mass_ = new double[capacity_]();

    // Accelerations (36 total: 3 types × 3 dims × 4 orders)
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            acc_total_[d][o] = new double[capacity_]();
            acc_reg_[d][o] = new double[capacity_]();
            acc_irr_[d][o] = new double[capacity_]();
        }
    }

    // New position/velocity (6)
    new_pos_x_ = new double[capacity_]();
    new_pos_y_ = new double[capacity_]();
    new_pos_z_ = new double[capacity_]();
    new_vel_x_ = new double[capacity_]();
    new_vel_y_ = new double[capacity_]();
    new_vel_z_ = new double[capacity_]();

    // Neighbor radius (1)
    neighbor_radius_sq_ = new double[capacity_]();

    // ========================================================================
    // Important fields (13 total)
    // ========================================================================
    // Timestep doubles (4)
    current_time_irr_ = new double[capacity_]();
    current_time_reg_ = new double[capacity_]();
    time_step_irr_ = new double[capacity_]();
    time_step_reg_ = new double[capacity_]();

    // Block times - ull_t (6)
    current_block_irr_ = new ull_t[capacity_]();
    new_current_block_irr_ = new ull_t[capacity_]();
    current_block_reg_ = new ull_t[capacity_]();
    time_block_irr_ = new ull_t[capacity_]();
    time_block_reg_ = new ull_t[capacity_]();
    next_block_irr_ = new ull_t[capacity_]();

    // Neighbor ints (3)
    num_neighbors_ = new int[capacity_]();
    new_num_neighbors_ = new int[capacity_]();
    neighbors_offset_ = new int[capacity_]();

    // ========================================================================
    // Low priority fields (10 total)
    // ========================================================================
    // IDs and indices (3)
    pid_ = new int[capacity_]();
    particle_index_ = new int[capacity_]();
    particle_type_ = new int[capacity_]();

    // Time levels (2)
    time_level_irr_ = new int[capacity_]();
    time_level_reg_ = new int[capacity_]();

    // Bool flags (3)
    is_active_ = new bool[capacity_]();
    is_up_to_date_ = new bool[capacity_]();
    is_cm_particle_ = new bool[capacity_]();

    // Other doubles (2)
    radius_ = new double[capacity_]();
    delta_mass_ = new double[capacity_]();
}

// ============================================================================
// Deallocate: Free all arrays and reset to nullptr
// ============================================================================
void ParticleData::deallocate() {
    if (capacity_ == 0) return;

    // ========================================================================
    // Critical double arrays
    // ========================================================================
    delete[] pos_x_; pos_x_ = nullptr;
    delete[] pos_y_; pos_y_ = nullptr;
    delete[] pos_z_; pos_z_ = nullptr;

    delete[] vel_x_; vel_x_ = nullptr;
    delete[] vel_y_; vel_y_ = nullptr;
    delete[] vel_z_; vel_z_ = nullptr;

    delete[] mass_; mass_ = nullptr;

    // Accelerations
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            delete[] acc_total_[d][o]; acc_total_[d][o] = nullptr;
            delete[] acc_reg_[d][o]; acc_reg_[d][o] = nullptr;
            delete[] acc_irr_[d][o]; acc_irr_[d][o] = nullptr;
        }
    }

    delete[] new_pos_x_; new_pos_x_ = nullptr;
    delete[] new_pos_y_; new_pos_y_ = nullptr;
    delete[] new_pos_z_; new_pos_z_ = nullptr;
    delete[] new_vel_x_; new_vel_x_ = nullptr;
    delete[] new_vel_y_; new_vel_y_ = nullptr;
    delete[] new_vel_z_; new_vel_z_ = nullptr;

    delete[] neighbor_radius_sq_; neighbor_radius_sq_ = nullptr;

    // ========================================================================
    // Important fields
    // ========================================================================
    delete[] current_time_irr_; current_time_irr_ = nullptr;
    delete[] current_time_reg_; current_time_reg_ = nullptr;
    delete[] time_step_irr_; time_step_irr_ = nullptr;
    delete[] time_step_reg_; time_step_reg_ = nullptr;

    delete[] current_block_irr_; current_block_irr_ = nullptr;
    delete[] new_current_block_irr_; new_current_block_irr_ = nullptr;
    delete[] current_block_reg_; current_block_reg_ = nullptr;
    delete[] time_block_irr_; time_block_irr_ = nullptr;
    delete[] time_block_reg_; time_block_reg_ = nullptr;
    delete[] next_block_irr_; next_block_irr_ = nullptr;

    delete[] num_neighbors_; num_neighbors_ = nullptr;
    delete[] new_num_neighbors_; new_num_neighbors_ = nullptr;
    delete[] neighbors_offset_; neighbors_offset_ = nullptr;

    // ========================================================================
    // Low priority fields
    // ========================================================================
    delete[] pid_; pid_ = nullptr;
    delete[] particle_index_; particle_index_ = nullptr;
    delete[] particle_type_; particle_type_ = nullptr;

    delete[] time_level_irr_; time_level_irr_ = nullptr;
    delete[] time_level_reg_; time_level_reg_ = nullptr;

    delete[] is_active_; is_active_ = nullptr;
    delete[] is_up_to_date_; is_up_to_date_ = nullptr;
    delete[] is_cm_particle_; is_cm_particle_ = nullptr;

    delete[] radius_; radius_ = nullptr;
    delete[] delta_mass_; delta_mass_ = nullptr;

    capacity_ = 0;
    count_ = 0;
}
