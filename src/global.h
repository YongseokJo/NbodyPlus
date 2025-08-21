#ifndef GLOBAL_H
#define GLOBAL_H
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "./FewBody/Group.h"
#include "performance.h"
#include <mpi.h>

#ifdef SEVN
#include "IO.h"
#include <map>
#endif


extern Particle *particles;
extern Particle *particles_original;

/* Communicators */
extern MPI_Win win;
extern MPI_Win win2;
extern MPI_Win win3;
extern MPI_Comm shared_comm;
extern int MyRank;
extern int NumberOfProcessor;
extern int NumberOfWorker;
const int ROOT = 0;
extern int NumberOfCommunication;
extern GlobalVariable *global_variable;
extern GlobalVariable *global_variable_original;
extern MPI_Datatype QueueType;
extern MPI_Datatype IparticleType;
extern MPI_Datatype JparticleType;

extern int *ActiveIndexToOriginalIndex_orginal;
extern int *ActiveIndexToOriginalIndex;

extern int LastParticleIndex;
extern int NumberOfParticle;
extern int NewPID;

extern int FixNumNeighbor;
extern double InitialNeighborRadius;


// Task
const int TASK_TAG = 1;
const int PTCL_TAG = 2;
const int TIME_TAG = 3;
const int QUEUE_TAG = 4;
const int ANY_TAG = 100;
const int TERMINATE_TAG = 666;

// Time
extern double global_time;
extern double global_time_irr;
extern ULL NextRegTimeBlock;
extern int time_block;
extern double time_step;
extern ULL block_max;
extern double eta;

// Enzo to Nbody
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTimeStep;


extern double E_binary;
extern double E_binary_SD;
extern double E_merger;
extern double E_PN;


// i/o
extern char* fname;
extern double endTime;
extern bool restart;
extern char* foutput;
extern double outputTime;
extern int outNum;
extern double outputTimeStep;
extern char* config_file;

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
