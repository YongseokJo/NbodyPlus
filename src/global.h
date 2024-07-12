#include <vector>
// #include <mpi.h>
#include <iostream>
#include "Particle/Particle.h"
#include "Binary/Binary.h"
#include <stdio.h>
#include <stdexcept>

#ifdef time_trace
#include "TimeTrace.h"
extern TimeTracer _time;
#endif



extern std::vector<int> LevelList;
extern int NNB;
extern int newNNB;
//extern int NumNeighborMax;

// Time
extern double global_time;
extern double global_time_irr;
extern ULL NextRegTimeBlock;
extern int time_block;
extern double time_step;
extern ULL block_max;


extern double binary_time;
extern double binary_time_prev;
extern ULL binary_block;

// ComputationChain
extern std::vector<Particle*> ComputationChain;
extern Particle* FirstComputation;
extern std::vector<Particle*> ComputationList;
extern int ComputationTimeMarker;
extern std::vector<Particle*> RegularList;
extern std::vector<Particle*> BinaryCandidateList;
extern std::vector<Binary*> BinaryList; // List of binaries to calculate

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

extern FILE* binout;
