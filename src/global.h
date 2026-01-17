#ifndef GLOBAL_H
#define GLOBAL_H

#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "./FewBody/Group.h"
#include "performance.h"
#include "profiler.h"
#include <mpi.h>
#include <unordered_map>
#include <string>

#ifdef SEVN
#include "IO.h"
#include <map>
#endif

// ============================================================================
// Constants
// ============================================================================
const int ROOT = 0;

// MPI Tags
const int TASK_TAG      = 1;
const int PTCL_TAG      = 2;
const int TIME_TAG      = 3;
const int QUEUE_TAG     = 4;
const int ANY_TAG       = 100;
const int TERMINATE_TAG = 666;

// ============================================================================
// MPI communicators and process info
// ============================================================================
extern int my_rank;
extern int num_processors;
extern int num_workers;
extern MPI_Comm shared_comm;

// ============================================================================
// MPI windows for shared memory
// ============================================================================
extern MPI_Win win;
extern Particle *particles;

extern MPI_Win win2;
extern GlobalVariable *g_state;

extern MPI_Win win3;
extern int* neighbors;
extern int* neighbors_original;

extern MPI_Win win4;
extern int* new_neighbors;
extern int* new_neighbors_original;

// ============================================================================
// MPI custom datatypes
// ============================================================================
extern MPI_Datatype queue_type_mpi;
extern MPI_Datatype iparticle_type_mpi;
extern MPI_Datatype jparticle_type_mpi;

// ============================================================================
// Particle array management
// ============================================================================
extern int last_particle_index;
extern int num_particles;
extern int new_cm_pid;

// ============================================================================
// Simulation parameters (from config file)
// ============================================================================
extern double eta;                      // Time step accuracy parameter
extern int fixed_num_neighbors;         // Fixed neighbor count (if used)
extern double initial_neighbor_radius;  // Initial neighbor search radius
extern double r_search;                 // Few-body search radius (code units)
extern double t_search;                 // Few-body search time (code units)

// ============================================================================
// Output settings
// ============================================================================
extern bool use_compression;            // Enable HDF5 compression
extern int compression_level;           // GZIP compression level (1-9)

// ============================================================================
// Restart settings
// ============================================================================
extern bool restart_enabled;
extern std::string checkpoint_file;

// ============================================================================
// Few-body tracking
// ============================================================================
extern std::unordered_map<int, int> cm_particle_worker_map;
extern std::unordered_map<int, int> prev_cm_particle_worker_map;

// ============================================================================
// Time variables
// ============================================================================
extern double global_time;
extern double global_time_irr;
extern ull_t next_reg_time_block;
extern int time_block;
extern double time_step;
extern ull_t block_max;
extern double end_time;
extern double enzo_time_step;

// ============================================================================
// I/O variables
// ============================================================================
extern char* fname;
extern bool restart;
extern char* foutput;
extern double output_time;
extern int output_num;
extern double output_time_step;
extern char* config_file;

// ============================================================================
// Energy tracking
// ============================================================================
extern double energy_binary;
extern double energy_binary_sd;
extern double energy_merger;
extern double energy_pn;

// ============================================================================
// Output file handles
// ============================================================================
extern FILE* bin_output_file;
extern FILE* merger_output_file;
extern FILE* worker_output_file;

#ifdef SEVN
extern FILE* sevn_output_file;
extern IO* sevnio;
extern std::multimap<double, int> sevn_list;
#endif

// ============================================================================
// Performance tracking
// ============================================================================
#ifdef PERFORMANCETRACE
extern Performance performance;
#endif

#endif
