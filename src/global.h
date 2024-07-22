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
extern REAL global_time;
extern REAL global_time_irr;
extern ULL NextRegTimeBlock;
extern int time_block;
extern REAL time_step;
extern ULL block_max;


extern REAL binary_time;
extern REAL binary_time_prev;
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
extern REAL inputTime;
extern REAL endTime;
extern bool restart;


// Enzo to Nbody
extern Particle* FirstEnzoParticle;
extern REAL EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern REAL EnzoTime;
extern REAL EnzoTimeStep;


// i/o
extern char* foutput;
extern bool IsOutput;
extern REAL outputTime;
extern REAL outputTimeStep;
extern int outNum;
//
//

extern FILE* binout;
