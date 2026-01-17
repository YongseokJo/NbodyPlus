#ifndef PARTICLE_H
#define PARTICLE_H

#include "def.h"
#include <cmath>
#include "cstring"
#include <vector>
#include "cuda/cuda_defs.h"
#include "global_state.h"

// SDAR
#include "Common/Float.h"
#include <iostream>
#include <iomanip>
#ifdef SEVN
#include "star.h"
#ifdef SEVN_BINARY
#include "binstar.h"
#endif
#endif

// External declarations
extern double initial_neighbor_radius;
extern double enzo_time_step;
extern GlobalVariable *g_state;
struct Particle;
extern Particle *particles;
extern double time_step;

// ============================================================================
// Binary interrupt state enumeration
// ============================================================================
enum class binary_interrupt_state_t : int {
    none              = 0,
    form              = 1,
    exchange          = 2,
    collision_candidate = 3,
    collision         = 4,
    threebody         = 5,  // particles from 3-body interaction
    manybody          = 6,  // many-body (>3) group terminated
    kicked            = 7,  // kicked or exploded as PISN
    terminated        = 8,  // CM particle only; terminated
    merger            = 9   // CM particle only; merger inside
};

// Scoped alias for backward compatibility with SDAR code
using BinaryInterruptState = binary_interrupt_state_t;

#define BINARY_STATE_ID_SHIFT 4
#define BINARY_INTERRUPT_STATE_MASK 0xF

// ============================================================================
// Particle structure
// ============================================================================
struct Group;

struct Particle {
    // ========================================================================
    // Core particle properties
    // ========================================================================
    int pid;                    // Particle ID
    int particle_index;         // Index in particle array
    int particle_type;          // Star type (see STAR_TYPE_* constants)
    double position[DIM];       // Position in 3D space (x, y, z)
    double velocity[DIM];       // Velocity in 3D space (vx, vy, vz)
    double mass;                // Particle mass

    // ========================================================================
    // Acceleration arrays [dimension][derivative_order]
    // Order: 0=acc, 1=jerk, 2=snap, 3=crackle
    // ========================================================================
    double acc_total[DIM][HERMITE_ORDER];     // Total acceleration
    double acc_regular[DIM][HERMITE_ORDER];   // Regular (far-field) acceleration
    double acc_irregular[DIM][HERMITE_ORDER]; // Irregular (near-field) acceleration

    // ========================================================================
    // Neighbor information
    // ========================================================================
    int num_neighbors;          // Current number of neighbors
    int new_num_neighbors;      // New neighbor count after update
    int neighbors_offset;       // Offset into global neighbor array

    // ========================================================================
    // Time stepping - irregular
    // ========================================================================
    double current_time_irr;    // Current irregular time
    ull_t current_block_irr;    // Current irregular block number
    ull_t new_current_block_irr;// New irregular block after update
    ull_t next_block_irr;       // Next irregular block
    double time_step_irr;       // Irregular time step
    ull_t time_block_irr;       // Irregular time block
    int time_level_irr;         // Irregular time level

    // ========================================================================
    // Time stepping - regular
    // ========================================================================
    double current_time_reg;    // Current regular time
    ull_t current_block_reg;    // Current regular block number
    double time_step_reg;       // Regular time step
    ull_t time_block_reg;       // Regular time block
    int time_level_reg;         // Regular time level

    // ========================================================================
    // Predicted quantities
    // ========================================================================
    double new_position[DIM];   // Predicted position
    double new_velocity[DIM];   // Predicted velocity
    double neighbor_radius_sq;  // Squared neighbor search radius
    double background_acc[DIM]; // Background acceleration

    // ========================================================================
    // SDAR/few-body properties
    // ========================================================================
    bool is_active;             // Whether particle is active
    bool is_up_to_date;         // Whether particle state is current
    double radius;              // Physical radius (also used by SEVN)
    double delta_mass;          // Mass to distribute to gas cells
    double time_check;          // Time to check next interrupt
    long long int binary_state; // Binary state (low bits: interrupt state, high bits: pair ID)
    double spin_param[DIM];     // Dimensionless spin parameter
    Group* group_info;          // Pointer to group information
    bool is_cm_particle;        // Is this a center-of-mass particle?
    int cm_particle_index;      // Index of parent CM particle
    int num_members;            // Number of group members (CM only)
    int members[10];            // Member particle indices (CM only)
    int new_num_members;        // New member count
    int new_members[10];        // New member indices

#ifdef SEVN
    StarSEVN* stellar_evolution;    // SEVN stellar evolution object
#ifdef SEVN_BINARY
    Binstar* binary_evolution;      // SEVN binary evolution object
#endif
    double formation_time;          // Formation time in Myr
    double world_time;              // World time in Myr
#endif

    // ========================================================================
    // Constructor
    // ========================================================================
    Particle() {
        position[0] = position[1] = position[2] = 0.0;
        velocity[0] = velocity[1] = velocity[2] = 0.0;
        pid = -1;
        mass = 0;
        neighbor_radius_sq = -1;
        num_neighbors = 0;
        new_num_neighbors = 0;
        neighbors_offset = -1;
        particle_type = STAR_TYPE_NO_FEEDBACK;
        current_time_irr = 0.0;
        current_time_reg = 0.0;
        current_block_irr = 0;
        current_block_reg = 0;
        next_block_irr = 0;
        time_step_irr = 0;
        time_step_reg = 0;
        time_level_irr = 0;
        time_level_reg = 0;
        time_block_irr = 0;
        time_block_reg = 0;

        for (int i = 0; i < DIM; i++) {
            velocity[i] = 0.0;
            position[i] = 0.0;
            new_position[i] = 0.0;
            new_velocity[i] = 0.0;
            background_acc[i] = 0.0;
            for (int j = 0; j < HERMITE_ORDER; j++) {
                acc_total[i][j] = 0.0;
                acc_regular[i][j] = 0.0;
                acc_irregular[i][j] = 0.0;
            }
            spin_param[i] = 0.0;
        }

        is_active = true;
        particle_index = -1;
        radius = 0.0;
        delta_mass = 0.0;
        time_check = NUMERIC_FLOAT_MAX;
        set_binary_interrupt_state(binary_interrupt_state_t::none);
        group_info = nullptr;
        is_cm_particle = false;
        is_up_to_date = true;
        cm_particle_index = -1;
        num_members = 0;
        new_num_members = 0;

#ifdef SEVN
        stellar_evolution = nullptr;
#ifdef SEVN_BINARY
        binary_evolution = nullptr;
#endif
        formation_time = 0.0;
        world_time = 0.0;
#endif
    }

    // ========================================================================
    // Initialization
    // ========================================================================
    void initialize(double *data, int p_id) {
        this->pid = p_id;
        this->position[0] = data[0];
        this->position[1] = data[1];
        this->position[2] = data[2];
        this->velocity[0] = data[3];
        this->velocity[1] = data[4];
        this->velocity[2] = data[5];
        this->mass = data[6];
        this->particle_type = STAR_TYPE_NO_FEEDBACK;
        this->current_time_reg = 0;
        this->current_time_irr = 0;
        this->neighbor_radius_sq = initial_neighbor_radius * initial_neighbor_radius;
        this->num_neighbors = 0;
        this->new_num_neighbors = 0;
        this->is_active = true;
        this->particle_index = p_id;
        this->neighbors_offset = this->particle_index * MAX_NUM_NEIGHBOR;
        this->delta_mass = 0.0;
        this->time_check = NUMERIC_FLOAT_MAX;
        this->set_binary_interrupt_state(binary_interrupt_state_t::none);
        this->group_info = nullptr;
        this->is_cm_particle = false;
        this->spin_param[0] = 0.0;
        this->spin_param[1] = 0.0;
        this->spin_param[2] = 0.0;
        this->cm_particle_index = -1;
        this->is_up_to_date = true;
        this->num_members = 0;
        this->new_num_members = 0;

#ifndef SEVN
        this->particle_type = STAR_TYPE_NO_FEEDBACK;
        this->radius = 2.25461e-8 / POSITION_UNIT * pow(this->mass * 1e9, 1.0/3.0);
#endif
    }

    // ========================================================================
    // Clear/reset particle
    // ========================================================================
    void clear() {
        pid = -1;
        particle_index = -1;
        num_neighbors = 0;
        new_num_neighbors = 0;
        neighbors_offset = -1;
        is_up_to_date = true;
        is_active = false;
        group_info = nullptr;
        is_cm_particle = false;
        cm_particle_index = -1;
        num_members = 0;
        new_num_members = 0;
        set_binary_interrupt_state(binary_interrupt_state_t::none);
        particle_type = STAR_TYPE_NO_FEEDBACK;
        spin_param[0] = 0.0;
        spin_param[1] = 0.0;
        spin_param[2] = 0.0;
    }

    // ========================================================================
    // Unit normalization
    // ========================================================================
    void normalize_particle() {
        this->mass *= 1e9;
        this->mass /= MASS_UNIT;
        for (int dim = 0; dim < DIM; dim++) {
            this->position[dim] *= 1000;  // kpc to pc
            this->position[dim] /= POSITION_UNIT;
            this->velocity[dim] *= 1e5 * YR_TO_SEC / PC_TO_CM;  // km/s to pc/yr
            this->velocity[dim] /= VELOCITY_UNIT;
        }
    }

    // ========================================================================
    // Update particle position/velocity
    // ========================================================================
    void update_particle() {
        if (this->is_cm_particle) {
            Particle* member_ptcl;
            for (int i = 0; i < this->num_members; i++) {
                member_ptcl = &particles[this->members[i]];
                for (int dim = 0; dim < DIM; dim++) {
                    member_ptcl->position[dim] += this->new_position[dim] - this->position[dim];
                    member_ptcl->velocity[dim] += this->new_velocity[dim] - this->velocity[dim];
                }
            }
        }

        for (int dim = 0; dim < DIM; dim++) {
            this->position[dim] = this->new_position[dim];
            this->velocity[dim] = this->new_velocity[dim];
        }
    }

    // ========================================================================
    // Second-order prediction
    // ========================================================================
    template <typename T>
    void predict_particle_second_order(double dt, T pos[], T vel[]) {
        dt = dt * enzo_time_step;

        if (dt == 0) {
            for (int dim = 0; dim < DIM; dim++) {
                pos[dim] = static_cast<T>(position[dim]);
                vel[dim] = static_cast<T>(velocity[dim]);
            }
        } else {
            for (int dim = 0; dim < DIM; dim++) {
                pos[dim] = static_cast<T>(((acc_total[dim][1] * dt / 3 + acc_total[dim][0]) * dt / 2 + velocity[dim]) * dt + position[dim]);
                vel[dim] = static_cast<T>((acc_total[dim][1] * dt / 2 + acc_total[dim][0]) * dt + velocity[dim]);
            }
        }
    }

    void predict_particle_second_order(double dt, std::vector<j_particle_t>& jparticles,
                                        std::vector<i_particle_t>& iparticles,
                                        std::vector<int>& local_regular_list) {
        dt = dt * enzo_time_step;

        i_particle_t iptcl;
        j_particle_t jptcl;

        if (dt == 0) {
            jptcl.pos_x = static_cast<cuda_real_t>(position[0]);
            jptcl.pos_y = static_cast<cuda_real_t>(position[1]);
            jptcl.pos_z = static_cast<cuda_real_t>(position[2]);
            jptcl.vel_x = static_cast<cuda_real_t>(velocity[0]);
            jptcl.vel_y = static_cast<cuda_real_t>(velocity[1]);
            jptcl.vel_z = static_cast<cuda_real_t>(velocity[2]);
            jptcl.mass = static_cast<cuda_real_t>(mass);
            jptcl.index = particle_index;
            jparticles.push_back(jptcl);
        } else {
            jptcl.pos_x = static_cast<cuda_real_t>(((acc_total[0][1] * dt / 3 + acc_total[0][0]) * dt / 2 + velocity[0]) * dt + position[0]);
            jptcl.pos_y = static_cast<cuda_real_t>(((acc_total[1][1] * dt / 3 + acc_total[1][0]) * dt / 2 + velocity[1]) * dt + position[1]);
            jptcl.pos_z = static_cast<cuda_real_t>(((acc_total[2][1] * dt / 3 + acc_total[2][0]) * dt / 2 + velocity[2]) * dt + position[2]);
            jptcl.vel_x = static_cast<cuda_real_t>((acc_total[0][1] * dt / 2 + acc_total[0][0]) * dt + velocity[0]);
            jptcl.vel_y = static_cast<cuda_real_t>((acc_total[1][1] * dt / 2 + acc_total[1][0]) * dt + velocity[1]);
            jptcl.vel_z = static_cast<cuda_real_t>((acc_total[2][1] * dt / 2 + acc_total[2][0]) * dt + velocity[2]);
            jptcl.mass = static_cast<cuda_real_t>(mass);
            jptcl.index = particle_index;
            jparticles.push_back(jptcl);
        }

        if (current_block_reg + time_block_reg == g_state->next_reg_time_block) {
            iptcl.pos_x = jptcl.pos_x;
            iptcl.pos_y = jptcl.pos_y;
            iptcl.pos_z = jptcl.pos_z;
            iptcl.vel_x = jptcl.vel_x;
            iptcl.vel_y = jptcl.vel_y;
            iptcl.vel_z = jptcl.vel_z;
            iptcl.radius_sq = static_cast<cuda_real_t>(neighbor_radius_sq);
            iptcl.dt_reg = static_cast<cuda_real_t>(time_block_reg * time_step * enzo_time_step);
            iparticles.push_back(iptcl);
            local_regular_list.push_back(particle_index);
        }
    }

    // ========================================================================
    // Method declarations (implemented elsewhere)
    // ========================================================================
    void correct_particle_fourth_order(double dt, double pos[], double vel[], double a[3][4]);
    void update_radius();
    void initialize_time_step();
    void compute_acceleration_irr();
    void compute_acceleration_reg();
    void calculate_time_step_irr();
    void calculate_time_step_irr_v2();
    void calculate_time_step_reg();
    void update_regular_particle_cuda();

    // SDAR group checking
    void check_new_group();
    void check_new_group_v2();
    void check_new_group_v3();
    void check_new_group_v4();

    // ========================================================================
    // Binary state manipulation
    // ========================================================================
    void set_binary_pair_id(const int id) {
        binary_state = (binary_state & BINARY_INTERRUPT_STATE_MASK) | (id << BINARY_STATE_ID_SHIFT);
    }

    void set_binary_interrupt_state(const binary_interrupt_state_t state) {
        binary_state = ((binary_state >> BINARY_STATE_ID_SHIFT) << BINARY_STATE_ID_SHIFT) | static_cast<int>(state);
    }

    binary_interrupt_state_t get_binary_interrupt_state() const {
        return static_cast<binary_interrupt_state_t>(binary_state & BINARY_INTERRUPT_STATE_MASK);
    }

    int get_binary_pair_id() const {
        return (binary_state >> BINARY_STATE_ID_SHIFT);
    }

    // ========================================================================
    // Accessors
    // ========================================================================
    double* get_position() { return position; }
    double* get_velocity() { return velocity; }

    // SDAR compatibility methods
    double* getPos() { return position; }
    double* getVel() { return velocity; }

    // SDAR compatibility member accessors (using methods instead of reference members to allow assignment operator)
    int& NumberOfMember() { return num_members; }
    int& NumberOfNeighbor() { return num_neighbors; }
    const int& NumberOfMember() const { return num_members; }
    const int& NumberOfNeighbor() const { return num_neighbors; }

    // ========================================================================
    // Output functions
    // ========================================================================
    static void print_column_title(std::ostream& fout, const int width = 20) {
        fout << std::setw(width) << "mass"
             << std::setw(width) << "pos.x"
             << std::setw(width) << "pos.y"
             << std::setw(width) << "pos.z"
             << std::setw(width) << "vel.x"
             << std::setw(width) << "vel.y"
             << std::setw(width) << "vel.z"
             << std::setw(width) << "radius"
             << std::setw(width) << "id";
    }
    // SDAR compatibility alias
    static void printColumnTitle(std::ostream& fout, const int width = 20) { print_column_title(fout, width); }

    void print_column(std::ostream& fout, const int width = 20) {
        fout << std::setw(width) << mass
             << std::setw(width) << position[0]
             << std::setw(width) << position[1]
             << std::setw(width) << position[2]
             << std::setw(width) << velocity[0]
             << std::setw(width) << velocity[1]
             << std::setw(width) << velocity[2]
             << std::setw(width) << radius
             << std::setw(width) << pid;
    }
    // SDAR compatibility alias
    void printColumn(std::ostream& fout, const int width = 20) { print_column(fout, width); }

    void copy_new_members(Particle* ptcl) {
        this->new_num_members = ptcl->new_num_members;
        std::memcpy(this->new_members, ptcl->new_members, sizeof(int) * ptcl->new_num_members);
    }

    void print_particle_info(FILE* file) {
        fprintf(file, "pid: %d, particle_index: %d, particle_type: %d\n", pid, particle_index, particle_type);
        fprintf(file, "position (pc): %e, %e, %e\n", position[0] * POSITION_UNIT, position[1] * POSITION_UNIT, position[2] * POSITION_UNIT);
        fprintf(file, "velocity (km/s): %e, %e, %e\n",
                velocity[0] * VELOCITY_UNIT / YR_TO_SEC * PC_TO_CM / 1e5,
                velocity[1] * VELOCITY_UNIT / YR_TO_SEC * PC_TO_CM / 1e5,
                velocity[2] * VELOCITY_UNIT / YR_TO_SEC * PC_TO_CM / 1e5);
        fprintf(file, "mass (Msun): %e\n", mass * MASS_UNIT);
        fprintf(file, "num_neighbors: %d, neighbor_radius (pc): %e\n", num_neighbors, sqrt(neighbor_radius_sq) * POSITION_UNIT);
        fprintf(file, "acc_total (0): %e, %e, %e\n", acc_total[0][0], acc_total[1][0], acc_total[2][0]);
        fprintf(file, "acc_regular (0): %e, %e, %e\n", acc_regular[0][0], acc_regular[1][0], acc_regular[2][0]);
        fprintf(file, "acc_irregular (0): %e, %e, %e\n", acc_irregular[0][0], acc_irregular[1][0], acc_irregular[2][0]);
        fprintf(file, "time_step_irr (Myr): %e, time_step_reg (Myr): %e\n", time_step_irr * enzo_time_step * 1e4, time_step_reg * enzo_time_step * 1e4);
    }
};

#endif
