#ifndef QUEUE_H
#define QUEUE_H

#include <cstdint>
#include <iostream>

// ============================================================================
// Task type enumeration
// ============================================================================
enum task_name_t : int8_t {
    // Force calculations
    TASK_IRR_FORCE        = 0,
    TASK_REG_FORCE        = 1,
    TASK_IRR_UPDATE       = 2,
    TASK_REG_UPDATE       = 3,
    TASK_REG_CUDA         = 4,

    // Initialization
    TASK_INIT_ACC_1       = 7,
    TASK_INIT_ACC_2       = 8,
    TASK_INIT_TIME        = 9,
    TASK_TIME_SYNC        = 10,

    // Few-body / group operations
    TASK_SEARCH_PRIMORDIAL_GROUP = 20,
    TASK_SEARCH_GROUP     = 22,
    TASK_MAKE_PRIMORDIAL_GROUP = 23,
    TASK_MAKE_GROUP       = 24,
    TASK_DELETE_GROUP     = 25,
    TASK_AR_INTEGRATION   = 26,
    TASK_MERGE_MANYBODY   = 27,

    // MPI operations
    TASK_CALC_ACC_01_MPI  = 28,
    TASK_CALC_ACC_23_MPI  = 29,
    TASK_GET_TOTAL_ENERGY = 30,
    TASK_PREPARE_GPU_CALC = 31,

    // Profiling (Phase 21.5)
    TASK_SEND_PROFILING   = 40,   // Request workers to send profiling data to root

    // Control
    TASK_SYNCHRONIZE      = 100,
    TASK_END              = -100,
    TASK_ERROR            = -1
};

// ============================================================================
// Queue structure for task communication
// ============================================================================
struct Queue {
    task_name_t task;       // Task type to execute
    int pid;                // Particle ID (target of task)
    double next_time;       // Next time for time-based tasks

    void print() {
        std::cout << "Task: " << static_cast<int>(task)
                  << ", PID: " << pid
                  << ", Next Time: " << next_time << std::endl;
    }
};

#endif
