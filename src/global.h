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

const int ROOT = 0;

// Task
const int TASK_TAG = 1;
const int PTCL_TAG = 2;
const int TIME_TAG = 3;
const int QUEUE_TAG = 4;
const int ANY_TAG = 100;
const int TERMINATE_TAG = 666;

/* Communicators */
extern int MyRank;
extern int NumberOfProcessor;
extern int NumberOfWorker;

extern MPI_Comm shared_comm;

extern MPI_Win win;
extern Particle *particles;

extern MPI_Win win2;
extern GlobalVariable *global_variable;

extern MPI_Win win3;
extern int* Neighbors;
extern int* Neighbors_original;

extern MPI_Win win4;
extern int* NewNeighbors;
extern int* NewNeighbors_original;

extern MPI_Datatype QueueType;
extern MPI_Datatype IparticleType;
extern MPI_Datatype JparticleType;

// Particle array
extern int LastParticleIndex;
extern int NumberOfParticle;
extern int NewCMPID;

// Parameters determined in config file
extern double eta;
extern int FixNumNeighbor;
extern double InitialNeighborRadius;
extern double RSearch;              // Few-body search radius (code units)
extern double TSearch;              // Few-body search time (code units)

// Output settings
extern bool UseCompression;         // Enable HDF5 compression
extern int CompressionLevel;        // GZIP compression level (1-9)

// Restart settings
extern bool RestartEnabled;
extern std::string CheckpointFile;

// Few-body
extern std::unordered_map<int, int> CMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11
extern std::unordered_map<int, int> PrevCMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11

// Time
extern double global_time;
extern double global_time_irr;
extern ULL NextRegTimeBlock;
extern int time_block;
extern double time_step;
extern ULL block_max;
extern double endTime;
extern double EnzoTimeStep;

// i/o
extern char* fname;
extern bool restart;
extern char* foutput;
extern double outputTime;
extern int outNum;
extern double outputTimeStep;
extern char* config_file;

// Energy tracking
extern double E_binary;
extern double E_binary_SD;
extern double E_merger;
extern double E_PN;

extern FILE* binout;
extern FILE* mergerout;
#ifdef SEVN
extern FILE* SEVNout;
extern IO* sevnio;
extern std::multimap<double, int> SEVNList;
#endif
extern FILE* workerout;

#ifdef PERFORMANCETRACE
extern Performance performance;
#endif

#endif
