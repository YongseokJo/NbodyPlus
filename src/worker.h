#ifndef WORKER_H
#define WORKER_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "mpi.h"
#include "global.h"
#include "queue.h"

// ============================================================================
// Worker structure for MPI task distribution
// ============================================================================

#define MAX_QUEUE 1000

struct Worker {
    // ========================================================================
    // Member variables
    // ========================================================================
    int rank;                           // MPI rank of this worker
    bool on_duty;                       // Whether worker is currently processing
    std::unordered_set<int> cm_particle_ids;  // CM particles assigned to this worker
    bool is_cm_worker;                  // Whether this worker handles CM particles
    Queue queues[MAX_QUEUE];            // Queue of tasks
    short num_queues;                   // Number of pending queues
    short current_queue;                // Index of current queue

    // Async MPI infrastructure (Phase 11)
    MPI_Request send_request;           // Request handle for async send
    MPI_Request recv_request;           // Request handle for async receive
    int result_buffer;                  // Dedicated buffer for async receive

    // Legacy aliases
    #define MyRank rank
    #define onDuty on_duty
    #define CMPtclIDs cm_particle_ids
    #define isCMWorker is_cm_worker
    #define NumberOfQueues num_queues
    #define CurrentQueue current_queue

    // ========================================================================
    // Constructor
    // ========================================================================
    Worker() {
        _initialize();
    }

    // ========================================================================
    // Initialization methods
    // ========================================================================
    void initialize() {
        on_duty = false;
        if (cm_particle_ids.size() > 0) {
            is_cm_worker = true;
        }
        num_queues = 0;
        current_queue = 0;
    }

    void initialize(int worker_rank) {
        rank = worker_rank;
        on_duty = false;
        if (cm_particle_ids.size() > 0) {
            is_cm_worker = true;
        }
        num_queues = 0;
    }

    // ========================================================================
    // Queue management
    // ========================================================================
    void add_queue(Queue queue) {
        int index = current_queue + num_queues++;
        if (num_queues == MAX_QUEUE) {
            fprintf(stderr, "num_queues exceeds MAX_QUEUE!\n");
            exit(EXIT_FAILURE);
        }
        index %= MAX_QUEUE;
        queues[index] = queue;
    }

    // Legacy alias
    void addQueue(Queue queue) { add_queue(queue); }

    void run_queue() {
        if (on_duty) {
            std::cout << "Worker " << rank << " is already on duty" << std::endl;
            std::cerr << "Worker " << rank << " is already on duty" << std::endl;
            exit(1);
        }
        if (num_queues == 0) {
            fprintf(stderr, "There is no queue in this worker!\n");
            exit(EXIT_FAILURE);
        }
        send_task(queues[current_queue]);
    }

    // Legacy alias
    void runQueue() { run_queue(); }

    void remove_queue() {
        // Currently empty - placeholder for future implementation
    }

    // Legacy alias
    void removeQueue() { remove_queue(); }

    // ========================================================================
    // Callback handling
    // ========================================================================
    void callback() {
        int return_value;
        PROFILE_START(TimerID::MPIRecv);
        MPI_Recv(&return_value, 1, MPI_INT, this->rank, TERMINATE_TAG, MPI_COMM_WORLD, &_status);
        PROFILE_STOP(TimerID::MPIRecv);
        if (!on_duty) {
            fprintf(stderr, "Error: worker %d was not on duty.\n", this->rank);
            fprintf(stdout, "Error: worker %d was not on duty.\n", this->rank);
            exit(1);
        }
        on_duty = false;
        current_queue++;
        current_queue %= MAX_QUEUE;
        num_queues--;
    }

    // ========================================================================
    // Task sending
    // ========================================================================
    void send_task(Queue& queue) {
        PROFILE_START(TimerID::MPISend);
        MPI_Send(&queue, 1, queue_type_mpi, this->rank, QUEUE_TAG, MPI_COMM_WORLD);
        PROFILE_STOP(TimerID::MPISend);
        on_duty = true;
    }

    // Legacy alias
    void sendTask(Queue& queue) { send_task(queue); }

    // ========================================================================
    // Async task sending (Phase 11)
    // ========================================================================
    void send_task_async() {
        // Safety check: ensure previous send completed
        if (send_request != MPI_REQUEST_NULL) {
            PROFILE_START(TimerID::MPIWait);
            MPI_Wait(&send_request, MPI_STATUS_IGNORE);
            PROFILE_STOP(TimerID::MPIWait);
        }

        PROFILE_START(TimerID::MPIIsend);
        MPI_Isend(&queues[current_queue], 1, queue_type_mpi,
                  this->rank, QUEUE_TAG, MPI_COMM_WORLD, &send_request);
        PROFILE_STOP(TimerID::MPIIsend);
        on_duty = true;
    }

    // ========================================================================
    // Async receive posting (Phase 11)
    // ========================================================================
    void post_receive() {
        // Safety check: ensure previous receive completed
        if (recv_request != MPI_REQUEST_NULL) {
            PROFILE_START(TimerID::MPIWait);
            MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            PROFILE_STOP(TimerID::MPIWait);
        }

        PROFILE_START(TimerID::MPIIrecv);
        MPI_Irecv(&result_buffer, 1, MPI_INT, this->rank,
                  TERMINATE_TAG, MPI_COMM_WORLD, &recv_request);
        PROFILE_STOP(TimerID::MPIIrecv);
    }

    // Wait for pending async send to complete (Phase 11)
    void wait_send_complete() {
        if (send_request != MPI_REQUEST_NULL) {
            PROFILE_START(TimerID::MPIWait);
            MPI_Wait(&send_request, MPI_STATUS_IGNORE);
            PROFILE_STOP(TimerID::MPIWait);
            send_request = MPI_REQUEST_NULL;
        }
    }

    // Handle completion from async receive (Phase 11)
    // Called after QueueScheduler::waitQueueAsync() identifies this worker
    void callback_async() {
        // recv_request already completed via MPI_Waitany
        recv_request = MPI_REQUEST_NULL;

        if (!on_duty) {
            fprintf(stderr, "Error: worker %d was not on duty (async).\n", this->rank);
            fprintf(stdout, "Error: worker %d was not on duty (async).\n", this->rank);
            exit(1);
        }
        on_duty = false;
        current_queue++;
        current_queue %= MAX_QUEUE;
        num_queues--;
    }

    // ========================================================================
    // Accessors
    // ========================================================================
    Queue* get_current_queue() { return &queues[current_queue]; }
    Queue* getCurrentQueue() { return get_current_queue(); }

private:
    MPI_Request _request;
    MPI_Status _status;

    void _initialize() {
        rank = -1;
        on_duty = false;
        cm_particle_ids.clear();
        is_cm_worker = false;
        num_queues = 0;
        current_queue = 0;

        // Initialize async request handles (Phase 11)
        send_request = MPI_REQUEST_NULL;
        recv_request = MPI_REQUEST_NULL;
        result_buffer = 0;
    }
};

extern Worker* workers;

#endif
