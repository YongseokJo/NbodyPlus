#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <unordered_set>

#include "def.h"
#include "global.h"
#include "QueueScheduler.h"
#include "cuda/cuda_defs.h"

void RegularWorker(int NumTargetTotal, int Jstart, int Jend, int gpu_id);
void AllocateDeviceMemory(int N_i, int N_j, int gpu_id);
void SendToDeviceMPI(int N_i, int N_j, int gpu_id);
void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList);
void RegularRoot(int NumTargetTotal, CUDA_REAL Acceleration[], int NumNeighbor[], int *NeighborList);

#ifdef Unuse

void RegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueScheduler &queue_scheduler){
	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];
    
	int *ACListReceive, *NumNeighborReceive;
	CUDA_REAL *Acceleration;
	int J_start, J_end;

	Particle *ptcl;
	double new_time = NextRegTimeBlock*time_step;  // next regular time

	Acceleration		= new CUDA_REAL[6*ListSize];
	NumNeighborReceive  = new int[ListSize];
	ACListReceive		= new int[ListSize * MaxNumNeighbor];

	std::memset(Acceleration, 0, ListSize * 6 * sizeof(CUDA_REAL));

	// Run this in worker's routine in queue_scheduler
	// two workers launch sendAllParitclesToGPU using MPI
	// replace this with the queue_scheduler
	/* // This is the original code that uses MPI to send particles to GPU
	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	assert(mpi_size > deviceCount);  // at least ranks 0 (root), 1, 2 (workers)
	if (mpi_rank > 0 && mpi_rank <= deviceCount) {
		int gpu_id = mpi_rank - 1; // temporarily
		int N_start = (gpu_id * ListSize) / deviceCount;
		int N_end = ((gpu_id + 1) * ListSize) / deviceCount;
		sendAllParticlesToGPU(new_time, RegularList, IndexList, N_start, N_end, gpu_id);
		// But this code, only processeses bound to GPUs participate in to the preprocessing
	}
	// synchronize all ranks before proceeding
	MPI_Barrier(MPI_COMM_WORLD);
	*/

	Queue queue;

	// queue_scheduler.initialize(RegSend);
	/*
	queue.task = RegSend;
	queue.next_time = new_time;
	for (int p = 1; p <= deviceCount; p++) {
		queue.pid = p - 1;
		MPI_Send(&queue, 1, QueueType, p, QUEUE_TAG, MPI_COMM_WORLD);
		MPI_Send()
	}
	*/
	sendAllParticlesToGPU(new_time, RegularList, IndexList);


	queue_scheduler.initialize(RegCal);
	queue.task = RegCal;
	queue.next_time = -1;
	int* int_lists = new int[3];
	for (int p = 1; p <= deviceCount; p++) {
		int gpu_id = p - 1;
		J_start = (gpu_id * ListSize) / deviceCount;
		J_end = ((gpu_id + 1) * ListSize) / deviceCount;

		queue.pid = gpu_id;
		MPI_Send(&queue, 1, QueueType, p, QUEUE_TAG, MPI_COMM_WORLD);

		int_lists[0] = ListSize;
		int_lists[1] = J_start;
		int_lists[2] = J_end;
		MPI_Send(int_lists, 3, MPI_INT, p, 1, MPI_COMM_WORLD);
	}
	delete int_lists;
	// Start RegularWorker in the queue_scheduler
	// RegularWorker(ListSize, J_start, J_end, gpu_id);

	RegularRoot(ListSize, Acceleration, NumNeighborReceive, ACListReceive);

	int completed = 0;
	MPI_Status status;
	while (completed < deviceCount) {
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int completed_rank = status.MPI_SOURCE;
		int return_value;
		MPI_Recv(&return_value, 1, MPI_INT, completed_rank, TERMINATE_TAG, MPI_COMM_WORLD, &status);
		completed++;
	}

	for (int i=0; i<ListSize; i++) {
		ptcl = &particles[ActiveIndexToOriginalIndex[IndexList[i]]];
		ptcl->NewNumberOfNeighbor = NumNeighborReceive[i];
		std::memcpy(ptcl->NewNeighbors, &ACListReceive[i * MaxNumNeighbor], NumNeighborReceive[i] * sizeof(int));

		for (int dim=0; dim<Dim; dim++) {
			ptcl->a_irr[dim][0]  = static_cast<double>(Acceleration[_six*i+dim]);
			ptcl->a_irr[dim][1] = static_cast<double>(Acceleration[_six*i+dim+3]);
		}
	}


	queue_scheduler.initialize(RegCuda);
	queue_scheduler.takeQueueRegularList(RegularList);
	do
	{
		queue_scheduler.assignQueueAutoRegularList();
		queue_scheduler.runQueueAuto();
		queue_scheduler.waitQueue(0); // blocking wait
	} while (queue_scheduler.isComplete());

	delete[] IndexList;
	delete[] Acceleration;
	delete[] NumNeighborReceive;
	delete[] ACListReceive;
}


void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList) {

	// Create a vector of indices from 0 to LastParticleIndex
	std::vector<int> indices(global_variable->LastParticleIndex + 1);
	std::iota(indices.begin(), indices.end(), 0);

	// Shuffle the indices randomly
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);

	// variables for saving variables to send to GPU
	
	CUDA_REAL * Radius2;
	//int size = NumberOfParticle;
	int size=0, j=0, i=0;
	int N_target = RegularList.size();

	// allocate memory to the temporary variables
	int N_j = N_end - N_start;
	h_ptcl_j     = new CUDA_REAL[N_j*7];
	h_ptcl_i     = new CUDA_REAL[N_target*7];
	Radius2  = new CUDA_REAL[N_target];		
	
	Particle *ptcl;
	double Position[3];
	double Velocity[3];

	int J_start, J_end;
	// copy the data of particles to the arrays to be sent
	
	Queue queue;

	queue.task = RegSend;
	queue.next_time = -1;
	// queue.pid = N_j;


	for (int idx: indices) {
		ptcl       = &particles[idx];
		size++;
		if (!ptcl->isActive) {
			// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
			continue;
		}

		if (RegularList.find(idx) != RegularList.end()) {
			IndexList[j] = size;
			Radius2[j] = (CUDA_REAL)ptcl->RadiusOfNeighbor; // mass weight?

			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position, Velocity);
			h_ptcl_i[j] = (CUDA_REAL) Position[0];
			h_ptcl_i[j + N_target] = (CUDA_REAL) Position[1];
			h_ptcl_i[j + 2 * N_target] = (CUDA_REAL) Position[2];
			h_ptcl_i[j + 3 * N_target] = (CUDA_REAL) Velocity[0];
			h_ptcl_i[j + 4 * N_target] = (CUDA_REAL) Velocity[1];
			h_ptcl_i[j + 5 * N_target] = (CUDA_REAL) Velocity[2];
			j++;
		}

		ActiveIndexToOriginalIndex[size] = idx;
		size++;
	}

	assert(NumberOfParticle == size); // for debugging by EW 2025.1.25

	for (int p = 1; p <= deviceCount; p++) {
		int gpu_id = p - 1;
		J_start = (gpu_id * ListSize) / deviceCount;
		J_end = ((gpu_id + 1) * ListSize) / deviceCount;

		for (int i = J_start; i < J_end; ++i) {
			int idx = indices[i];
			ptcl = &particles[idx];

			if (!ptcl->isActive) {
				// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
				continue;
			}

			// h_ptcl_j[size + k * NumberOfParticle] = PREDICTED POSITION AND VELCOITY;
			h_ptcl_j[i + 6 * N_j] = (CUDA_REAL)ptcl->Mass;
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position, Velocity);
			h_ptcl_j[i] 			= (CUDA_REAL) Position[0];
			h_ptcl_j[i + 	 N_j]	= (CUDA_REAL) Position[1];
			h_ptcl_j[i + 2 * N_j]	= (CUDA_REAL) Position[2];
			h_ptcl_j[i + 3 * N_j]	= (CUDA_REAL) Velocity[0];
			h_ptcl_j[i + 4 * N_j]	= (CUDA_REAL) Velocity[1];
			h_ptcl_j[i + 5 * N_j]	= (CUDA_REAL) Velocity[2];
		}
		// write this for me
		queue.pid = gpu_id;
		MPI_Send(&queue, 1, QueueType, p, QUEUE_TAG, MPI_COMM_WORLD);
		MPI_Send(N_target, 1, MPI_INT, p, 1010, MPI_COMM_WORLD);
		MPI_Send(N_j, 1, MPI_INT, p, 1011, MPI_COMM_WORLD);
		MPI_Send(h_ptcl_i, 7 * N_target, MPI_CUDA, p, 1001, MPI_COMM_WORLD);

		// send the “j”‐particles block
		MPI_Send(h_ptcl_j, 7 * N_j, MPI_CUDA, p, 1002, MPI_COMM_WORLD);
		MPI_Send(IndexList, N_target, MPI_INT, p, 1003, MPI_COMM_WORLD);
		MPI_Send(Radius2, N_target, MPI_CUDA, p, 1004, MPI_COMM_WORLD);
	}

	int completed = 0;
	MPI_Status status;
	while (completed < deviceCount) {
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int completed_rank = status.MPI_SOURCE;
		int return_value;
		MPI_Recv(&return_value, 1, MPI_INT, completed_rank, TERMINATE_TAG, MPI_COMM_WORLD, &status);
		completed++;
	}

	delete[] h_ptcl_j;
	delete[] h_ptcl_i;
	delete[] Radius2;
}

#endif