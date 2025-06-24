#ifndef QUEUE_H
#define QUEUE_H

#include <cstdint>

enum TaskName : int8_t {
    IrrForce = 0,
    RegForce = 1,
    IrrUpdate = 2, 
    RegUpdate = 3,
    RegCuda = 4,
    RegCudaUpdate = 5,
    InitAcc1 = 7,
    InitAcc2 = 8,
    InitTime = 9,
    TimeSync = 10, 
    SearchPrimordialGroup = 20,
    SearchGroup = 22,
    MakePrimordialGroup = 23,
    MakeGroup = 24,
    DeleteGroup = 25,
    ARIntegration = 26,
    MergeManyBody = 27,
#ifdef MultiNode
    UpdateInitAcc01 = 30,
    UpdateInitAcc23 = 31,
    SendNewNeighbors = 32,
    UpdateLastParticleIndex = 33,
    UpdateTimeVariables = 34,
    UpdateTimeCorrection = 35,
    UpdateIrregularForce = 36,
    UpdateBinaryMerger = 37,
    UpdateFBTermination = 38,
    UpdateNewGroup = 39,
    UpdateManybodyMerger = 40,
    UpdateActiveIndexToOriginalIndex = 41,
    UpdateBeforeRegCuda = 42,
    UpdateAfterRegCuda = 43,
    UpdateAfterRegCudaUpdate = 44,
#ifdef SEVN
    UpdateStellarEvolution0 = 45,
    UpdateStellarEvolution1 = 46,
#endif
#endif
    Synchronize = 100,
    Ends = -100,
    Error = -1
};

/* This should contain all the information that should be transmitted to workers */
struct Queue {
    TaskName task;
    int pid;
    double next_time;
    
    void print() {
        std::cout << "Task: " << task << ", PID: " << pid << ", Next Time: " << next_time << std::endl;
    }
};

#ifdef MultiNode
struct UpdateInitAcc1 {
    int pid;
    int numberofneighbors;
    double airr0[3];
    double airr1[3];
    double areg0[3];
    double areg1[3];
};

struct UpdateInitAcc2 {
    int pid;
    double airr2[3];
    double airr3[3];
    double areg2[3];
    double areg3[3];
};

struct UpdateTime {
    int pid;
    double timestep_reg;
    ULL timeblock_reg;
    int timelevel_reg;
    double timestep_irr;
    ULL timeblock_irr;
    int timelevel_irr;
    double currenttime_irr;
    double currenttime_reg;
    ULL currentblock_irr;
    ULL currentblock_reg;
};

struct UpdateTimeCorr {
    int pid;
    double timestep_irr;
    ULL timeblock_irr;
    int timelevel_irr;
};

struct UpdateIrrForce {
    int pid;
    int newnumberofneighbors;
    int newneighbors[5];
    double newposition[3];
    double newvelocity[3];
    double airr[3][4];
    double atot[3][4];
    ULL newcurrentblock_irr;
    int timelevel_irr;
    double timestep_irr;
    ULL timeblock_irr;
    ULL nextblock_irr;
};

struct UpdateBinary {
    int pid;
    long long int binary_state;
    double position[3];
    double velocity[3];
    double mass;
    double currenttime_irr;
};

struct UpdateFBTerm {
    int pid;
    long long int binary_state;
    ULL currentblock_irr;
    double currenttime_irr;
    ULL currentblock_reg;
    double currenttime_reg;
    ULL newcurrentblock_irr;
    ULL nextblock_irr;
    int timelevel_irr;
    double timestep_irr;
    ULL timeblock_irr;
    int timelevel_reg;
    double timestep_reg;
    ULL timeblock_reg;
    double radiusofneighbor;
    double airr[3][4];
    double areg[3][4];
    int numberofneighbors;
    int neighbors[MaxNumNeighbor];
    double position[3];
    double velocity[3];
};

struct UpdateNewCM {
    int PID; // not ParticleIndex, but Particle ID
    int numberofmembers;
    int members[5];
    double position[3];
    double velocity[3];
    double mass;
    double radiusofneighbor;
    double currenttime_irr;
    double currenttime_reg;
    ULL currentblock_irr;
    ULL currentblock_reg;
    ULL newcurrentblock_irr;
    ULL nextblock_irr;
    double timestep_irr;
    ULL timeblock_irr;
    int timelevel_irr;
    double timestep_reg;
    ULL timeblock_reg;
    int timelevel_reg;
    int numberofneighbors;
    int neighbors[MaxNumNeighbor];
    double airr[3][4];
    double areg[3][4];
};

struct UpdateRegCuda {
    int pid;
    double newposition[3];
    double newvelocity[3];
    double airr[3][2];
    double areg[3][4];
    double atot[3][4];
    int newnumberofneighbors;
    int newneighbors[MaxNumNeighbor];
};

struct UpdateRegCudaUpdate {
    int pid;
    ULL currentblock_reg;
    double currenttime_reg;
    int timelevel_reg;
    double timestep_reg;
    ULL timeblock_reg;
    int timelevel_irr;
    double timestep_irr;
    ULL timeblock_irr;
    double radiusofneighbor;
    ULL nextblock_irr;
};
#ifdef SEVN
struct UpdateSEVN0 {
    int pid;
    int particletype;
    double mass;
    double radius;
};

struct UpdateSEVN1 {
    int pid;
    double dm;
    double mass;
    double radius;
    double velocity[3];
    long long int binary_state;
    double a_spin[3];
    char isactive;
};
#endif
#endif
#endif