#include <iostream>
#include <stdio.h>
#include "global.h"
#ifdef SEVN
#include <map>
#endif
#include <unordered_map>

// ============================================================================
// MPI communicator variables
// ============================================================================
int my_rank;
int num_processors;
int num_workers;
MPI_Comm shared_comm;

// ============================================================================
// MPI windows for shared memory
// ============================================================================
MPI_Win win;
Particle *particles;

MPI_Win win2;
GlobalVariable *g_state;

MPI_Win win3;
int* neighbors;

MPI_Win win4;
int* new_neighbors;

// ============================================================================
// MPI custom datatypes
// ============================================================================
MPI_Datatype queue_type_mpi;
MPI_Datatype iparticle_type_mpi;
MPI_Datatype jparticle_type_mpi;

// ============================================================================
// Particle array management
// ============================================================================
int last_particle_index;    // Index of last particle in array
int num_particles;          // Number of active particles (single + cm)
int new_cm_pid;             // Next PID for new CM particle

// ============================================================================
// Simulation parameters (from config file)
// ============================================================================
double eta;
int fixed_num_neighbors;
double initial_neighbor_radius;
double r_search;            // Few-body search radius (code units)
double t_search;            // Few-body search time (code units)

// ============================================================================
// Output settings
// ============================================================================
bool use_compression = true;    // Enable HDF5 compression by default
int compression_level = 6;      // GZIP compression level (1-9)

// ============================================================================
// Restart settings
// ============================================================================
bool restart_enabled = false;
std::string checkpoint_file;

// ============================================================================
// Few-body tracking
// ============================================================================
std::unordered_map<int, int> cm_particle_worker_map;
std::unordered_map<int, int> prev_cm_particle_worker_map;

// ============================================================================
// Time variables
// ============================================================================
double global_time;
double global_time_irr;
ull_t next_reg_time_block;
int time_block;
double time_step;
ull_t block_max;
double end_time;
double enzo_time_step;

// ============================================================================
// I/O variables
// ============================================================================
char* fname;
bool restart;
char* foutput;
double output_time;
int output_num;
double output_time_step;
char* config_file;

// ============================================================================
// Energy tracking
// ============================================================================
double energy_binary;       // Total binary energy in processor
double energy_binary_sd;    // Total slowdown binary energy
double energy_merger;       // Total merger energy
double energy_pn;           // Total post-Newtonian energy

// ============================================================================
// Output file handles
// ============================================================================
FILE* bin_output_file;
FILE* merger_output_file;
FILE* worker_output_file;

#ifdef SEVN
FILE* sevn_output_file;
IO* sevnio = nullptr;
std::multimap<double, int> sevn_list;
#endif

#ifdef PERFORMANCETRACE
Performance performance;
#endif

// ============================================================================
// Initialize default global values
// ============================================================================
void DefaultGlobal() {
    // Timesteps
    end_time = 1;
    enzo_time_step = end_time / 1e10;  // end_time should be Myr
    output_time_step = output_time_step / end_time;

    time_block = -30;
    block_max = static_cast<ull_t>(pow(2, -time_block));
    time_step = std::pow(2, time_block);

    end_time = 0.0;
    output_time_step = 0.0;

    global_time = 0.0;
    output_time = 0.0;

    eta = 0.01;
    fixed_num_neighbors = 100;
    initial_neighbor_radius = 0.011;

    // Few-body search parameters (converted to code units)
    // Default: r_search = 2.5e-4 pc, t_search = 1e-6 Myr
    r_search = 2.5e-4 / POSITION_UNIT;
    t_search = 1e-6 / (enzo_time_step * 1e4);

    energy_binary = 0.0;
    energy_binary_sd = 0.0;
    energy_merger = 0.0;
    energy_pn = 0.0;

#ifdef SEVN
    std::vector<std::string> args = {
        "empty",  // Not used
        "-tables", "/home/vinicius/install/sevn_custom/tables/SEVNtracks_parsec_ov04_AGB",
        "-xspinmode", "geneva",
        "-hardmode", "disabled",
        "-tmode", "disabled"
    };
    std::vector<char*> c_args;
    for (auto& arg : args) {
        c_args.push_back(&arg[0]);
    }

    sevnio = new IO;
    sevnio->load(c_args.size(), c_args.data());
#endif
}
