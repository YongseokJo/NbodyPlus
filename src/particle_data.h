#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <cstddef>
#include "def.h"

// Forward declaration for sync methods
struct Particle;

// ============================================================================
// ParticleData: Structure of Arrays (SoA) container for particle data
// ============================================================================
// This class stores particle data in SoA layout for improved cache efficiency
// and GPU memory coalescing. All fields from the original Particle struct are
// converted to separate arrays.
//
// Usage:
//   ParticleData data;
//   data.allocate(1000);  // Allocate for 1000 particles
//   data.set_position(0, x, y, z);
//   double px = data.get_pos_x(0);
//   data.deallocate();
// ============================================================================

class ParticleData {
public:
    // ========================================================================
    // Lifecycle
    // ========================================================================
    ParticleData();
    ~ParticleData();

    // Disable copy semantics (memory is manually managed)
    ParticleData(const ParticleData&) = delete;
    ParticleData& operator=(const ParticleData&) = delete;

    // Memory management
    void allocate(size_t capacity);
    void deallocate();

    // Size accessors
    size_t size() const { return count_; }
    size_t capacity() const { return capacity_; }
    void set_size(size_t n) { count_ = n; }

    // ========================================================================
    // Position accessors (3 arrays)
    // ========================================================================
    double* pos_x() { return pos_x_; }
    double* pos_y() { return pos_y_; }
    double* pos_z() { return pos_z_; }
    const double* pos_x() const { return pos_x_; }
    const double* pos_y() const { return pos_y_; }
    const double* pos_z() const { return pos_z_; }

    double get_pos_x(size_t i) const { return pos_x_[i]; }
    double get_pos_y(size_t i) const { return pos_y_[i]; }
    double get_pos_z(size_t i) const { return pos_z_[i]; }
    void set_pos_x(size_t i, double val) { pos_x_[i] = val; }
    void set_pos_y(size_t i, double val) { pos_y_[i] = val; }
    void set_pos_z(size_t i, double val) { pos_z_[i] = val; }

    // Convenience 3D accessors
    void get_position(size_t i, double& x, double& y, double& z) const {
        x = pos_x_[i]; y = pos_y_[i]; z = pos_z_[i];
    }
    void set_position(size_t i, double x, double y, double z) {
        pos_x_[i] = x; pos_y_[i] = y; pos_z_[i] = z;
    }

    // ========================================================================
    // Velocity accessors (3 arrays)
    // ========================================================================
    double* vel_x() { return vel_x_; }
    double* vel_y() { return vel_y_; }
    double* vel_z() { return vel_z_; }
    const double* vel_x() const { return vel_x_; }
    const double* vel_y() const { return vel_y_; }
    const double* vel_z() const { return vel_z_; }

    double get_vel_x(size_t i) const { return vel_x_[i]; }
    double get_vel_y(size_t i) const { return vel_y_[i]; }
    double get_vel_z(size_t i) const { return vel_z_[i]; }
    void set_vel_x(size_t i, double val) { vel_x_[i] = val; }
    void set_vel_y(size_t i, double val) { vel_y_[i] = val; }
    void set_vel_z(size_t i, double val) { vel_z_[i] = val; }

    // Convenience 3D accessors
    void get_velocity(size_t i, double& x, double& y, double& z) const {
        x = vel_x_[i]; y = vel_y_[i]; z = vel_z_[i];
    }
    void set_velocity(size_t i, double x, double y, double z) {
        vel_x_[i] = x; vel_y_[i] = y; vel_z_[i] = z;
    }

    // ========================================================================
    // Mass accessor (1 array)
    // ========================================================================
    double* mass() { return mass_; }
    const double* mass() const { return mass_; }
    double get_mass(size_t i) const { return mass_[i]; }
    void set_mass(size_t i, double val) { mass_[i] = val; }

    // ========================================================================
    // Acceleration accessors (36 arrays total: 3 types × 3 dims × 4 orders)
    // Access pattern: acc_total_[dim][order] where dim=0,1,2 and order=0,1,2,3
    // ========================================================================
    double* acc_total(int dim, int order) { return acc_total_[dim][order]; }
    const double* acc_total(int dim, int order) const { return acc_total_[dim][order]; }
    double get_acc_total(size_t i, int dim, int order) const { return acc_total_[dim][order][i]; }
    void set_acc_total(size_t i, int dim, int order, double val) { acc_total_[dim][order][i] = val; }

    double* acc_reg(int dim, int order) { return acc_reg_[dim][order]; }
    const double* acc_reg(int dim, int order) const { return acc_reg_[dim][order]; }
    double get_acc_reg(size_t i, int dim, int order) const { return acc_reg_[dim][order][i]; }
    void set_acc_reg(size_t i, int dim, int order, double val) { acc_reg_[dim][order][i] = val; }

    double* acc_irr(int dim, int order) { return acc_irr_[dim][order]; }
    const double* acc_irr(int dim, int order) const { return acc_irr_[dim][order]; }
    double get_acc_irr(size_t i, int dim, int order) const { return acc_irr_[dim][order][i]; }
    void set_acc_irr(size_t i, int dim, int order, double val) { acc_irr_[dim][order][i] = val; }

    // ========================================================================
    // New position/velocity accessors (6 arrays) - prediction output
    // ========================================================================
    double* new_pos_x() { return new_pos_x_; }
    double* new_pos_y() { return new_pos_y_; }
    double* new_pos_z() { return new_pos_z_; }
    const double* new_pos_x() const { return new_pos_x_; }
    const double* new_pos_y() const { return new_pos_y_; }
    const double* new_pos_z() const { return new_pos_z_; }

    double get_new_pos_x(size_t i) const { return new_pos_x_[i]; }
    double get_new_pos_y(size_t i) const { return new_pos_y_[i]; }
    double get_new_pos_z(size_t i) const { return new_pos_z_[i]; }
    void set_new_pos_x(size_t i, double val) { new_pos_x_[i] = val; }
    void set_new_pos_y(size_t i, double val) { new_pos_y_[i] = val; }
    void set_new_pos_z(size_t i, double val) { new_pos_z_[i] = val; }

    void get_new_position(size_t i, double& x, double& y, double& z) const {
        x = new_pos_x_[i]; y = new_pos_y_[i]; z = new_pos_z_[i];
    }
    void set_new_position(size_t i, double x, double y, double z) {
        new_pos_x_[i] = x; new_pos_y_[i] = y; new_pos_z_[i] = z;
    }

    double* new_vel_x() { return new_vel_x_; }
    double* new_vel_y() { return new_vel_y_; }
    double* new_vel_z() { return new_vel_z_; }
    const double* new_vel_x() const { return new_vel_x_; }
    const double* new_vel_y() const { return new_vel_y_; }
    const double* new_vel_z() const { return new_vel_z_; }

    double get_new_vel_x(size_t i) const { return new_vel_x_[i]; }
    double get_new_vel_y(size_t i) const { return new_vel_y_[i]; }
    double get_new_vel_z(size_t i) const { return new_vel_z_[i]; }
    void set_new_vel_x(size_t i, double val) { new_vel_x_[i] = val; }
    void set_new_vel_y(size_t i, double val) { new_vel_y_[i] = val; }
    void set_new_vel_z(size_t i, double val) { new_vel_z_[i] = val; }

    void get_new_velocity(size_t i, double& x, double& y, double& z) const {
        x = new_vel_x_[i]; y = new_vel_y_[i]; z = new_vel_z_[i];
    }
    void set_new_velocity(size_t i, double x, double y, double z) {
        new_vel_x_[i] = x; new_vel_y_[i] = y; new_vel_z_[i] = z;
    }

    // ========================================================================
    // Neighbor radius squared accessor (1 array)
    // ========================================================================
    double* neighbor_radius_sq() { return neighbor_radius_sq_; }
    const double* neighbor_radius_sq() const { return neighbor_radius_sq_; }
    double get_neighbor_radius_sq(size_t i) const { return neighbor_radius_sq_[i]; }
    void set_neighbor_radius_sq(size_t i, double val) { neighbor_radius_sq_[i] = val; }

    // ========================================================================
    // Timestep double accessors (4 arrays)
    // ========================================================================
    double* current_time_irr() { return current_time_irr_; }
    const double* current_time_irr() const { return current_time_irr_; }
    double get_current_time_irr(size_t i) const { return current_time_irr_[i]; }
    void set_current_time_irr(size_t i, double val) { current_time_irr_[i] = val; }

    double* current_time_reg() { return current_time_reg_; }
    const double* current_time_reg() const { return current_time_reg_; }
    double get_current_time_reg(size_t i) const { return current_time_reg_[i]; }
    void set_current_time_reg(size_t i, double val) { current_time_reg_[i] = val; }

    double* time_step_irr() { return time_step_irr_; }
    const double* time_step_irr() const { return time_step_irr_; }
    double get_time_step_irr(size_t i) const { return time_step_irr_[i]; }
    void set_time_step_irr(size_t i, double val) { time_step_irr_[i] = val; }

    double* time_step_reg() { return time_step_reg_; }
    const double* time_step_reg() const { return time_step_reg_; }
    double get_time_step_reg(size_t i) const { return time_step_reg_[i]; }
    void set_time_step_reg(size_t i, double val) { time_step_reg_[i] = val; }

    // ========================================================================
    // Block time accessors (ull_t - 6 arrays)
    // ========================================================================
    ull_t* current_block_irr() { return current_block_irr_; }
    const ull_t* current_block_irr() const { return current_block_irr_; }
    ull_t get_current_block_irr(size_t i) const { return current_block_irr_[i]; }
    void set_current_block_irr(size_t i, ull_t val) { current_block_irr_[i] = val; }

    ull_t* new_current_block_irr() { return new_current_block_irr_; }
    const ull_t* new_current_block_irr() const { return new_current_block_irr_; }
    ull_t get_new_current_block_irr(size_t i) const { return new_current_block_irr_[i]; }
    void set_new_current_block_irr(size_t i, ull_t val) { new_current_block_irr_[i] = val; }

    ull_t* current_block_reg() { return current_block_reg_; }
    const ull_t* current_block_reg() const { return current_block_reg_; }
    ull_t get_current_block_reg(size_t i) const { return current_block_reg_[i]; }
    void set_current_block_reg(size_t i, ull_t val) { current_block_reg_[i] = val; }

    ull_t* time_block_irr() { return time_block_irr_; }
    const ull_t* time_block_irr() const { return time_block_irr_; }
    ull_t get_time_block_irr(size_t i) const { return time_block_irr_[i]; }
    void set_time_block_irr(size_t i, ull_t val) { time_block_irr_[i] = val; }

    ull_t* time_block_reg() { return time_block_reg_; }
    const ull_t* time_block_reg() const { return time_block_reg_; }
    ull_t get_time_block_reg(size_t i) const { return time_block_reg_[i]; }
    void set_time_block_reg(size_t i, ull_t val) { time_block_reg_[i] = val; }

    ull_t* next_block_irr() { return next_block_irr_; }
    const ull_t* next_block_irr() const { return next_block_irr_; }
    ull_t get_next_block_irr(size_t i) const { return next_block_irr_[i]; }
    void set_next_block_irr(size_t i, ull_t val) { next_block_irr_[i] = val; }

    // ========================================================================
    // Neighbor int accessors (3 arrays)
    // ========================================================================
    int* num_neighbors() { return num_neighbors_; }
    const int* num_neighbors() const { return num_neighbors_; }
    int get_num_neighbors(size_t i) const { return num_neighbors_[i]; }
    void set_num_neighbors(size_t i, int val) { num_neighbors_[i] = val; }

    int* new_num_neighbors() { return new_num_neighbors_; }
    const int* new_num_neighbors() const { return new_num_neighbors_; }
    int get_new_num_neighbors(size_t i) const { return new_num_neighbors_[i]; }
    void set_new_num_neighbors(size_t i, int val) { new_num_neighbors_[i] = val; }

    int* neighbors_offset() { return neighbors_offset_; }
    const int* neighbors_offset() const { return neighbors_offset_; }
    int get_neighbors_offset(size_t i) const { return neighbors_offset_[i]; }
    void set_neighbors_offset(size_t i, int val) { neighbors_offset_[i] = val; }

    // ========================================================================
    // ID/index int accessors (3 arrays)
    // ========================================================================
    int* pid() { return pid_; }
    const int* pid() const { return pid_; }
    int get_pid(size_t i) const { return pid_[i]; }
    void set_pid(size_t i, int val) { pid_[i] = val; }

    int* particle_index() { return particle_index_; }
    const int* particle_index() const { return particle_index_; }
    int get_particle_index(size_t i) const { return particle_index_[i]; }
    void set_particle_index(size_t i, int val) { particle_index_[i] = val; }

    int* particle_type() { return particle_type_; }
    const int* particle_type() const { return particle_type_; }
    int get_particle_type(size_t i) const { return particle_type_[i]; }
    void set_particle_type(size_t i, int val) { particle_type_[i] = val; }

    // ========================================================================
    // Time level int accessors (2 arrays)
    // ========================================================================
    int* time_level_irr() { return time_level_irr_; }
    const int* time_level_irr() const { return time_level_irr_; }
    int get_time_level_irr(size_t i) const { return time_level_irr_[i]; }
    void set_time_level_irr(size_t i, int val) { time_level_irr_[i] = val; }

    int* time_level_reg() { return time_level_reg_; }
    const int* time_level_reg() const { return time_level_reg_; }
    int get_time_level_reg(size_t i) const { return time_level_reg_[i]; }
    void set_time_level_reg(size_t i, int val) { time_level_reg_[i] = val; }

    // ========================================================================
    // Bool flag accessors (3 arrays)
    // ========================================================================
    bool* is_active() { return is_active_; }
    const bool* is_active() const { return is_active_; }
    bool get_is_active(size_t i) const { return is_active_[i]; }
    void set_is_active(size_t i, bool val) { is_active_[i] = val; }

    bool* is_up_to_date() { return is_up_to_date_; }
    const bool* is_up_to_date() const { return is_up_to_date_; }
    bool get_is_up_to_date(size_t i) const { return is_up_to_date_[i]; }
    void set_is_up_to_date(size_t i, bool val) { is_up_to_date_[i] = val; }

    bool* is_cm_particle() { return is_cm_particle_; }
    const bool* is_cm_particle() const { return is_cm_particle_; }
    bool get_is_cm_particle(size_t i) const { return is_cm_particle_[i]; }
    void set_is_cm_particle(size_t i, bool val) { is_cm_particle_[i] = val; }

    // ========================================================================
    // Other double accessors (2 arrays)
    // ========================================================================
    double* radius() { return radius_; }
    const double* radius() const { return radius_; }
    double get_radius(size_t i) const { return radius_[i]; }
    void set_radius(size_t i, double val) { radius_[i] = val; }

    double* delta_mass() { return delta_mass_; }
    const double* delta_mass() const { return delta_mass_; }
    double get_delta_mass(size_t i) const { return delta_mass_[i]; }
    void set_delta_mass(size_t i, double val) { delta_mass_[i] = val; }

    // ========================================================================
    // Particle sync methods - copy between AoS Particle and SoA ParticleData
    // ========================================================================

    // Copy all fields from a Particle struct to ParticleData at index i
    void sync_from_particle(const Particle& p, size_t i);

    // Copy all fields from ParticleData at index i back to Particle struct
    void sync_to_particle(Particle& p, size_t i) const;

    // ========================================================================
    // Bulk 3-vector accessors - return pointers for contiguous access
    // ========================================================================

    // Get position as array: out[0]=x, out[1]=y, out[2]=z
    void get_position_vec(size_t i, double out[3]) const {
        out[0] = pos_x_[i]; out[1] = pos_y_[i]; out[2] = pos_z_[i];
    }

    // Get velocity as array: out[0]=vx, out[1]=vy, out[2]=vz
    void get_velocity_vec(size_t i, double out[3]) const {
        out[0] = vel_x_[i]; out[1] = vel_y_[i]; out[2] = vel_z_[i];
    }

    // Get new_position as array
    void get_new_position_vec(size_t i, double out[3]) const {
        out[0] = new_pos_x_[i]; out[1] = new_pos_y_[i]; out[2] = new_pos_z_[i];
    }

    // Get new_velocity as array
    void get_new_velocity_vec(size_t i, double out[3]) const {
        out[0] = new_vel_x_[i]; out[1] = new_vel_y_[i]; out[2] = new_vel_z_[i];
    }

    // Set position from array
    void set_position_vec(size_t i, const double in[3]) {
        pos_x_[i] = in[0]; pos_y_[i] = in[1]; pos_z_[i] = in[2];
    }

    // Set velocity from array
    void set_velocity_vec(size_t i, const double in[3]) {
        vel_x_[i] = in[0]; vel_y_[i] = in[1]; vel_z_[i] = in[2];
    }

    // Set new_position from array
    void set_new_position_vec(size_t i, const double in[3]) {
        new_pos_x_[i] = in[0]; new_pos_y_[i] = in[1]; new_pos_z_[i] = in[2];
    }

    // Set new_velocity from array
    void set_new_velocity_vec(size_t i, const double in[3]) {
        new_vel_x_[i] = in[0]; new_vel_y_[i] = in[1]; new_vel_z_[i] = in[2];
    }

    // ========================================================================
    // Acceleration array accessors - get/set [3][4] acceleration arrays
    // ========================================================================

    // Get acceleration total as [3][4] array
    void get_acc_total_array(size_t i, double out[3][4]) const {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                out[d][o] = acc_total_[d][o][i];
            }
        }
    }

    // Set acceleration total from [3][4] array
    void set_acc_total_array(size_t i, const double in[3][4]) {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                acc_total_[d][o][i] = in[d][o];
            }
        }
    }

    // Get acceleration regular as [3][4] array
    void get_acc_reg_array(size_t i, double out[3][4]) const {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                out[d][o] = acc_reg_[d][o][i];
            }
        }
    }

    // Set acceleration regular from [3][4] array
    void set_acc_reg_array(size_t i, const double in[3][4]) {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                acc_reg_[d][o][i] = in[d][o];
            }
        }
    }

    // Get acceleration irregular as [3][4] array
    void get_acc_irr_array(size_t i, double out[3][4]) const {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                out[d][o] = acc_irr_[d][o][i];
            }
        }
    }

    // Set acceleration irregular from [3][4] array
    void set_acc_irr_array(size_t i, const double in[3][4]) {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                acc_irr_[d][o][i] = in[d][o];
            }
        }
    }

    // ========================================================================
    // Acceleration accumulation helpers
    // ========================================================================

    // Zero all acceleration components for particle i (total, reg, irr)
    void zero_all_acc(size_t i) {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                acc_total_[d][o][i] = 0.0;
                acc_reg_[d][o][i] = 0.0;
                acc_irr_[d][o][i] = 0.0;
            }
        }
    }

    // Zero total acceleration for particle i
    void zero_acc_total(size_t i) {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                acc_total_[d][o][i] = 0.0;
            }
        }
    }

    // Zero regular acceleration for particle i
    void zero_acc_reg(size_t i) {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                acc_reg_[d][o][i] = 0.0;
            }
        }
    }

    // Zero irregular acceleration for particle i
    void zero_acc_irr(size_t i) {
        for (int d = 0; d < 3; d++) {
            for (int o = 0; o < 4; o++) {
                acc_irr_[d][o][i] = 0.0;
            }
        }
    }

    // Add to total acceleration [dim][order] for particle i
    void add_to_acc_total(size_t i, int dim, int order, double val) {
        acc_total_[dim][order][i] += val;
    }

    // Add to regular acceleration [dim][order] for particle i
    void add_to_acc_reg(size_t i, int dim, int order, double val) {
        acc_reg_[dim][order][i] += val;
    }

    // Add to irregular acceleration [dim][order] for particle i
    void add_to_acc_irr(size_t i, int dim, int order, double val) {
        acc_irr_[dim][order][i] += val;
    }

    // Add 3-vector to total acceleration at order 0 (position-dependent force)
    void add_to_acc_total_vec(size_t i, const double a[3]) {
        acc_total_[0][0][i] += a[0];
        acc_total_[1][0][i] += a[1];
        acc_total_[2][0][i] += a[2];
    }

    // Add 3-vector to regular acceleration at order 0
    void add_to_acc_reg_vec(size_t i, const double a[3]) {
        acc_reg_[0][0][i] += a[0];
        acc_reg_[1][0][i] += a[1];
        acc_reg_[2][0][i] += a[2];
    }

    // Add 3-vector to irregular acceleration at order 0
    void add_to_acc_irr_vec(size_t i, const double a[3]) {
        acc_irr_[0][0][i] += a[0];
        acc_irr_[1][0][i] += a[1];
        acc_irr_[2][0][i] += a[2];
    }

protected:
    // ========================================================================
    // Capacity and count (protected for MPI subclass access)
    // ========================================================================
    size_t capacity_;
    size_t count_;

    // ========================================================================
    // Array pointers (protected for MPI subclass access)
    // ========================================================================
    // Critical fields (43 double arrays)
    // Position (3)
    double* pos_x_;
    double* pos_y_;
    double* pos_z_;

    // Velocity (3)
    double* vel_x_;
    double* vel_y_;
    double* vel_z_;

    // Mass (1)
    double* mass_;

    // Accelerations (36 total: 3 types × 3 dims × 4 orders)
    // [dim][order] where dim=0,1,2 (x,y,z) and order=0,1,2,3 (acc,jerk,snap,crackle)
    double* acc_total_[3][4];
    double* acc_reg_[3][4];
    double* acc_irr_[3][4];

    // New position/velocity (6)
    double* new_pos_x_;
    double* new_pos_y_;
    double* new_pos_z_;
    double* new_vel_x_;
    double* new_vel_y_;
    double* new_vel_z_;

    // Neighbor radius (1)
    double* neighbor_radius_sq_;

    // Important fields (13 arrays)
    // Timestep doubles (4)
    double* current_time_irr_;
    double* current_time_reg_;
    double* time_step_irr_;
    double* time_step_reg_;

    // Block times (ull_t) (6)
    ull_t* current_block_irr_;
    ull_t* new_current_block_irr_;
    ull_t* current_block_reg_;
    ull_t* time_block_irr_;
    ull_t* time_block_reg_;
    ull_t* next_block_irr_;

    // Neighbor ints (3)
    int* num_neighbors_;
    int* new_num_neighbors_;
    int* neighbors_offset_;

    // Low priority fields (10 arrays)
    // IDs and indices (3)
    int* pid_;
    int* particle_index_;
    int* particle_type_;

    // Time levels (2)
    int* time_level_irr_;
    int* time_level_reg_;

    // Bool flags (3)
    bool* is_active_;
    bool* is_up_to_date_;
    bool* is_cm_particle_;

    // Other doubles (2)
    double* radius_;
    double* delta_mass_;
};

#endif // PARTICLE_DATA_H
