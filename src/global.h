#include <vector>
// #include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"


#ifdef time_trace
#include "TimeTrace.h"
extern TimeTracer _time;
#endif



extern std::vector<int> LevelList;
extern int NNB;
extern int newNNB;
extern int NumNeighborMax;

// Time
extern double global_time;
extern ULL NextRegTimeBlock;
extern int time_block;
extern double time_step;
extern ULL block_max;

// ComputationChain
extern std::vector<Particle*> ComputationChain;
extern Particle* FirstComputation;
extern std::vector<Particle*> ComputationList;
extern int ComputationTimeMarker;
extern std::vector<int> RegIndexList; // how about changing this to particle list


//extern bool debug;
extern char* fname;
extern double inputTime;
extern double endTime;
extern bool restart;


// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoTime;
extern double EnzoTimeStep;


// i/o
extern char* foutput;
extern bool IsOutput;
extern double outputTime;
extern double outputTimeStep;
extern int outNum;
//
//


// MPI variables
//extern MPI_Comm inter_comm;
//extern MPI_Comm nbody_comm;


