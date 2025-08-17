#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "../global.h"
#include "../QueueScheduler.h"
#include "cuda_functions.h"
#include <cstring>
#include <omp.h>

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

void sendAllParticlesToGPU(double new_time, std::unordered_set<int>& RegularList, int *IndexList);

/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */
void calculateRegAccelerationOnGPU(std::unordered_set<int>& RegularList, QueueScheduler &queue_scheduler){

#ifdef PERFORMANCETRACE
	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
#endif

	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here
	CUDA_REAL (*AccRegReceive)[Dim]		= new CUDA_REAL[ListSize][Dim];
	CUDA_REAL (*AccRegDotReceive)[Dim]	= new CUDA_REAL[ListSize][Dim];
	double (*AccIrr)[Dim]				= new double[ListSize][Dim];
	double (*AccIrrDot)[Dim]			= new double[ListSize][Dim];

	int *NumNeighborReceive				= new int[ListSize];
	int *ACListReceive					= new int[ListSize * MaxNumNeighbor];

	Particle *ptcl;

	double new_time = NextRegTimeBlock*time_step;  // next regular time

	// (Query to MY) Do we have to initialize them to 0?
	std::memset(AccRegReceive,		0, ListSize * Dim * sizeof(CUDA_REAL));
	std::memset(AccRegDotReceive,	0, ListSize * Dim * sizeof(CUDA_REAL));
	std::memset(AccIrr,				0, ListSize * Dim * sizeof(double));
	std::memset(AccIrrDot,			0, ListSize * Dim * sizeof(double));

#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG
	std::cout << "sendAllParticlesToGPU starts" << std::endl;
#endif

#ifdef NSIGHT
	nvtxRangePushA("sendAllParticlesToGPU");
#endif
	sendAllParticlesToGPU(new_time, RegularList, IndexList);  // needs to be updated
#ifdef NSIGHT
	nvtxRangePop();
#endif

#ifdef DEBUG
	std::cout << "sendAllParticlesToGPU ended" << std::endl;
#endif

#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	performance.RegularSendAllParticlesToGPU +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG
	std::cout << "CalculateAccelerationOnDevice starts" << std::endl;
#endif
  
#ifdef NSIGHT
	nvtxRangePushA("CalculateAccelerationOnDevice");
#endif

	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive, AccRegDotReceive, NumNeighborReceive, ACListReceive);
  
#ifdef NSIGHT
	nvtxRangePop();
#endif
  
#ifdef DEBUG
	std::cout << "CalculateAccelerationOnDevice ended" << std::endl;
#endif

#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	performance.RegularGPU +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG
	std::cout << "Adjust Regular Gravity starts" << std::endl;
#endif

#ifdef NSIGHT
	nvtxRangePushA("RegCuda");
#endif

	for (int i=0; i<ListSize; i++) {
		ptcl = &particles[ActiveIndexToOriginalIndex[IndexList[i]]];

		ptcl->NewNumberOfNeighbor = NumNeighborReceive[i];
		std::memcpy(ptcl->NewNeighbors, &ACListReceive[i * MaxNumNeighbor], NumNeighborReceive[i] * sizeof(int));

		for (int dim=0; dim<Dim; dim++) {
#ifdef CUDA_FLOAT
			ptcl->a_irr[dim][0] = static_cast<double>(AccRegReceive[i][dim]);		// Just temporarilly save new reg acc here!
			ptcl->a_irr[dim][1] = static_cast<double>(AccRegDotReceive[i][dim]);	// Just temporarilly save new reg acc here!
#else
			ptcl->a_irr[dim][0] = AccRegReceive[i][dim];		// Just temporarilly save new reg acc here!
			ptcl->a_irr[dim][1] = AccRegDotReceive[i][dim];		// Just temporarilly save new reg acc here!
#endif
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

/*
	// Adjust Regular Gravity
	int i=0;
	TaskName task=RegCuda;
	Queue queue = {task, -1, -1.0};
	queue_scheduler.initialize(RegCuda);
	queue_scheduler.takeQueueRegularList(RegularList);
	do
	{
		queue_scheduler.assignQueueRegularList();

        for (auto worker = queue_scheduler.WorkersToGo.begin(); worker != queue_scheduler.WorkersToGo.end();)
        {
            if ((*worker)->NumberOfQueues > 0) // original
            {
				//std::cout << "(REG_CUDA) My Rank =" << (*worker)->MyRank << std::endl;
				// queue_scheduler.sendQueueforRegCuda(*worker);
				// MPI_Send(&task, 1, MPI_INT, (*worker)->MyRank, TASK_TAG, MPI_COMM_WORLD);
				// MPI_Send(&ActiveIndexToOriginalIndex[IndexList[i]], 1, MPI_INT, (*worker)->MyRank, PTCL_TAG, MPI_COMM_WORLD);
				queue.pid = ActiveIndexToOriginalIndex[IndexList[i]];
				MPI_Send(&queue, 1, QueueType, (*worker)->MyRank, QUEUE_TAG, MPI_COMM_WORLD);
				MPI_Send(&NumNeighborReceive[i], 1, MPI_INT, (*worker)->MyRank, 10, MPI_COMM_WORLD);
				MPI_Send(&ACListReceive[i * MaxNumNeighbor], NumNeighborReceive[i], MPI_INT, (*worker)->MyRank, 11, MPI_COMM_WORLD);
				MPI_Send(&AccRegReceive[i][0], 3, MPI_DOUBLE, (*worker)->MyRank, 12, MPI_COMM_WORLD);
				MPI_Send(&AccRegDotReceive[i][0], 3, MPI_DOUBLE, (*worker)->MyRank, 13, MPI_COMM_WORLD);
				((*worker))->onDuty = true;
                worker = queue_scheduler.WorkersToGo.erase(worker);
				i++;
            }
			else
			{
                ++worker;
			}
        }
		queue_scheduler.waitQueue(0); // blocking wait
	} while (queue_scheduler.isComplete());
*/
#ifdef NSIGHT
	nvtxRangePop();
#endif

#ifdef DEBUG
	std::cout << "Adjust Regular Gravity ended" << std::endl;
#endif

#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	performance.RegularAdjust +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

	delete[] IndexList;

	delete[] AccRegReceive;
	delete[] AccRegDotReceive;
	delete[] AccIrr;
	delete[] AccIrrDot;

	delete[] NumNeighborReceive;
	delete[] ACListReceive;

	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends




// (Query MY) Let's optimize this function later. Copying data to h_ptcl in _ReceiveFromHost of cuda_my_acceleation.cpp seems super inefficient. 2025.5.24
void sendAllParticlesToGPU(double new_time, std::unordered_set<int>& RegularList, int *IndexList) {

	// variables for saving variables to send to GPU
	CUDA_REAL *Mass				= new CUDA_REAL[NumberOfParticle];
	CUDA_REAL *Mdot				= new CUDA_REAL[NumberOfParticle];	
	CUDA_REAL *Radius2			= new CUDA_REAL[NumberOfParticle];
	CUDA_REAL(*Position)[Dim]	= new CUDA_REAL[NumberOfParticle][Dim];
	CUDA_REAL(*Velocity)[Dim]	= new CUDA_REAL[NumberOfParticle][Dim];

	int size=0, j=0;

	/* new code using OpenMP by EW 2025.8.6
	for (int i = 0; i <= LastParticleIndex; i++) {
		Particle* ptcl = &particles[i];

		if (!ptcl->isActive)
			continue;

		if (RegularList.find(i) != RegularList.end())
			IndexList[j++] = size;
		
		ActiveIndexToOriginalIndex[size] = i;
		size++;
	}

	#pragma omp parallel for
	for (int i = 0; i < size; i++) {
		Particle* ptcl = &particles[ActiveIndexToOriginalIndex[i]];
		Mass[i] = (CUDA_REAL)ptcl->Mass;
		Mdot[i] = 0; // particle[i]->Mass;
		Radius2[i] = (CUDA_REAL)ptcl->RadiusOfNeighbor; // mass weight?
		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeReg, Position[i], Velocity[i]);
		else
			ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, Position[i], Velocity[i]);
	}
	*/

	// /* // original code not using OpenMP by EW 2025.8.6
	Particle *ptcl;

	Queue queue = {PrepareGPUCalc, -1, new_time};
	MPI_Request requests[NumberOfWorker];
	for (int i = 0; i < NumberOfWorker; i++) {
		MPI_Isend(&queue, 1, QueueType, i+1, QUEUE_TAG, MPI_COMM_WORLD, &requests[i]);
	}

	std::vector<int> num_elements(NumberOfWorker, 0);
	int num = 0;
	while (num < NumberOfWorker) {
		int J_start = num * (LastParticleIndex + 1) / NumberOfWorker;
		int J_end = (num+1) * (LastParticleIndex + 1) / NumberOfWorker;

		int nn = 0;
		for (int i = J_start; i < J_end; i++) {
			ptcl = &particles[i];

			if (!ptcl->isActive)
				continue;

			nn++;

			if (RegularList.find(i) != RegularList.end()) {
				IndexList[j] = size;
				j++;
			}

			Mass[size]    = (CUDA_REAL)ptcl->Mass;
			Mdot[size]    = 0; //particle[i]->Mass;
			Radius2[size] = (CUDA_REAL)ptcl->RadiusOfNeighbor; // mass weight?

			ActiveIndexToOriginalIndex[size] = i;
			size++;
		}
		num_elements[num] = nn;
		num++;
	}
	assert(NumberOfParticle == size);
	std::vector<int> displs(NumberOfWorker, 0);
	for (int i = 1; i < NumberOfWorker; i++)
		displs[i] = displs[i-1] + num_elements[i-1];

	MPI_Waitall(NumberOfWorker, requests, MPI_STATUSES_IGNORE);
	int completed = 0;
	MPI_Status status;
	while (completed < NumberOfWorker) {
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int completed_rank = status.MPI_SOURCE;

		MPI_Recv(Position[displs[completed_rank - 1]], num_elements[completed_rank - 1] * Dim, MPI_FLOAT, completed_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Velocity[displs[completed_rank - 1]], num_elements[completed_rank - 1] * Dim, MPI_FLOAT, completed_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		completed++;
	}


	/*
	// copy the data of particles to the arrays to be sent
	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl       = &particles[i];

		if (!ptcl->isActive) {
			// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
			continue;
		}

		if (RegularList.find(i) != RegularList.end()) {
			IndexList[j] = size;
			j++;
		}

		Mass[size]    = (CUDA_REAL)ptcl->Mass;
		Mdot[size]    = 0; //particle[i]->Mass;
		Radius2[size] = (CUDA_REAL)ptcl->RadiusOfNeighbor; // mass weight?

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position[size], Velocity[size]);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, Position[size], Velocity[size]);

		ActiveIndexToOriginalIndex[size] = i;
		// std::cout << "(size , i) = "  << size << " " << i << std::endl;
		size++;
	}
	*/

	assert(NumberOfParticle == size); // for debugging by EW 2025.1.25

	// fprintf(stdout, "in sendAllParticlesToGPU, NumberOfParticle = %d, size=%d, TotalNumberOfParticle=%d\n", NumberOfParticle, size, LastParticleIndex+1);


	//fprintf(stdout, "Sending particles to GPU...\n");
	//fflush(stdout);
	// send the arrays to GPU
	SendToDevice(&size, Mass, Position, Velocity, Radius2, Mdot);

	//fprintf(stdout, "Done.\n");
	//fflush(stdout);
	// free the temporary variables
	delete[] Mass;
	delete[] Mdot;
	delete[] Radius2;
	delete[] Position;
	delete[] Velocity;
}