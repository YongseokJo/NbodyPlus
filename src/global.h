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
extern int shared_rank;
extern int MyRank;
extern int NumberOfProcessor;
extern int NumberOfWorker;
const int ROOT = 0;
extern int NumberOfCommunication;
extern GlobalVariable *global_variable;
extern GlobalVariable *global_variable_original;
extern MPI_Datatype QueueType;

#ifdef MultiNode
extern int update_rank;
extern int NumberOfNode;
extern int* ranks_update_comm;
extern MPI_Comm update_comm;
extern MPI_Datatype UpdateInitAcc1Type;
extern MPI_Datatype UpdateInitAcc2Type;
extern MPI_Datatype UpdateTimeType;
extern MPI_Datatype UpdateTimeCorrType;
extern MPI_Datatype UpdateIrrForceType;
extern MPI_Datatype UpdateBinaryType;
extern MPI_Datatype UpdateFBTermType;
extern MPI_Datatype UpdateNewCMType;
extern MPI_Datatype UpdateRegCudaType;
extern MPI_Datatype UpdateRegCudaUpdateType;
#ifdef SEVN
extern MPI_Datatype UpdateSEVN0Type;
extern MPI_Datatype UpdateSEVN1Type;
#endif
#endif

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
extern int Task[NumberOfTask];

// Time
extern double global_time;
extern double global_time_irr;
extern ULL NextRegTimeBlock;
extern int time_block;
extern double time_step;
extern ULL block_max;
extern double eta;


extern double binary_time;
extern double binary_time_prev;
extern ULL binary_block;

// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTimeStep;



// i/o
extern char* fname;
extern double inputTime;
extern double endTime;
extern bool restart;
extern char* foutput;
extern bool IsOutput;
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
