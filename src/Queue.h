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
    CalcAcc01MPI = 28,
    CalcAcc23MPI = 29,
    GetTotalEnergy = 30,
    PrepareGPUCalc = 31,
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
#endif