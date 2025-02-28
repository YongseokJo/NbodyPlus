#ifndef QUEUE_MANAGER_H
#define QUEUE_MANAGER_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include "global.h"
#include "Worker.h"

class QueueManager
{
public:
    MPI_Request requests[MAX_COMMUNICATION];
    MPI_Status statuses[MAX_COMMUNICATION];
    int NumberOfCommunication;
    int completed_queues;

    QueueManager() {}

    void initialize(std::vector<int> list, TaskName task) {
        //MPI_Win_lock_all(0, win4);
        //MPI_Win_lock_all(0, win5);
        global_variable->QueueSize = list.size();
        for (int i=0; i<list.size(); i++) {
            queues[i] = list[i];
        }
        tasks[0] = task;
        completed_queues = 0;
        NumberOfCommunication = 0;
        MPI_Win_sync(win4);
        MPI_Win_sync(win5);
		MPI_Win_flush_all(win4);
		MPI_Win_flush_all(win5);
        //MPI_Win_unlock_all(win4);
        //MPI_Win_unlock_all(win5);
    }

    void initialize(std::unordered_set<int> list, TaskName task) {
        //MPI_Win_lock_all(0, win4);
        //MPI_Win_lock_all(0, win5);
        int i=0;
        global_variable->QueueSize = list.size();
        for (int elem:list) {
            queues[i++] = elem;
        }
        tasks[0] = task;
        completed_queues = 0;
        NumberOfCommunication = 0;
        MPI_Win_sync(win4);
        MPI_Win_sync(win5);
		MPI_Win_flush_all(win4);
		MPI_Win_flush_all(win5);
        //MPI_Win_unlock_all(win4);
        //MPI_Win_unlock_all(win5);
    }


    void signalWorkers() {
        for (int i=1; i<=NumberOfWorker; i++) {
            //fprintf(stdout, "signal workers %d!\n", i);
            MPI_Isend(NULL, 0, MPI_BYTE, i, TASK_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
            //MPI_Send(NULL, 0, MPI_BYTE, i, TASK_TAG, MPI_COMM_WORLD);
            //MPI_Send(&tasks[0], 1, MPI_INT, i, TASK_TAG, MPI_COMM_WORLD);
        }
        MPI_Waitall(NumberOfCommunication, requests, statuses);
        NumberOfCommunication = 0;
        MPI_Win_fence(0, win);
		MPI_Win_fence(0, win4);
		MPI_Win_fence(0, win5);
    }

    void signalWorker(int rank) {
        //fprintf(stdout, "signal worker %d!\n", rank);
        //fflush(stdout);
        MPI_Win_flush(rank, win4);
        MPI_Win_flush(rank, win5);
        MPI_Send(NULL, 0, MPI_BYTE, rank, TASK_TAG, MPI_COMM_WORLD);
		//MPI_Win_fence(0, win4);
		//MPI_Win_fence(0, win5);
    }

    void reportProgress () {
        while (completed_queues < global_variable->QueueSize)
        {
            MPI_Probe(MPI_ANY_SOURCE, FINISH_TAG, MPI_COMM_WORLD, &statuses[0]);
            MPI_Recv(NULL, 0, MPI_BYTE, statuses[0].MPI_SOURCE, FINISH_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            completed_queues++;
            //fprintf(stderr, "completed queues = (%d/%d)\n", completed_queues, global_variable->QueueSize);
        }
        MPI_Win_fence(0, win);
        //MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Win_fence(0, win);
        //MPI_Win_sync(win);
    }


};

#endif