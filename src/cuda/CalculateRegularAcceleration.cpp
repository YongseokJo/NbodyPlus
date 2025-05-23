#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "../global.h"
#include "../QueueScheduler.h"
#include "cuda_functions.h"
#include <cstring>

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int* data, int NumTask, int TAG);
void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList);
void CalculateAccelerationOnDevice(int *NumTargetTotal, int *h_target_list, double acc[][3], double adot[][3], int NumNeighbor[], int *NeighborList);

/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */
void calculateRegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueScheduler &queue_scheduler){

#ifdef PERFORMANCETRACE
	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
#endif

	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];

	// variables for saving variables to send to GPU
	// only regular particle informations are stored here
	double (*AccRegReceive)[Dim];
	double (*AccRegDotReceive)[Dim];
	double (*AccIrr)[Dim];
	double (*AccIrrDot)[Dim];
#ifdef CUDA_FLOAT
	CUDA_REAL (*AccRegReceive_f)[Dim];
	CUDA_REAL (*AccRegDotReceive_f)[Dim];
#endif 


	int *ACListReceive;
	int *NumNeighborReceive;

	Particle *ptcl;

	double new_time = NextRegTimeBlock*time_step;  // next regular time


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	AccRegReceive    = new double[ListSize][Dim];
	AccRegDotReceive = new double[ListSize][Dim];
	AccIrr           = new double[ListSize][Dim];
	AccIrrDot        = new double[ListSize][Dim];

#ifdef CUDA_FLOAT
	AccRegReceive_f		= new CUDA_REAL[ListSize][Dim];
	AccRegDotReceive_f	= new CUDA_REAL[ListSize][Dim];
#endif 
	NumNeighborReceive  = new int[ListSize];

	ACListReceive = new int[ListSize * MaxNumNeighbor];

	/* // original code; (Query to MY) Do we have to initialize them to 0?
	for (int i=0; i<ListSize; i++) {
		for (int dim=0; dim<Dim; dim++) {
			AccRegReceive[i][dim]    = 0;
			AccRegDotReceive[i][dim] = 0;
			AccIrr[i][dim]           = 0;
			AccIrrDot[i][dim]        = 0;
#ifdef CUDA_FLOAT
			AccRegReceive_f[i][dim]		= 0;
			AccRegDotReceive_f[i][dim]	= 0;
#endif 
		}
	}
	*/
	// /* // faster than the above code
	std::memset(AccRegReceive, 0, ListSize * Dim * sizeof(double));
	std::memset(AccRegDotReceive, 0, ListSize * Dim * sizeof(double));
	std::memset(AccIrr, 0, ListSize * Dim * sizeof(double));
	std::memset(AccIrrDot, 0, ListSize * Dim * sizeof(double));
#ifdef CUDA_FLOAT
	std::memset(AccRegReceive_f, 0, ListSize * Dim * sizeof(CUDA_REAL));
	std::memset(AccRegDotReceive_f, 0, ListSize * Dim * sizeof(CUDA_REAL));
#endif
	// */

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
  
#ifdef CUDA_FLOAT
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive_f,	AccRegDotReceive_f, NumNeighborReceive, ACListReceive);
#else
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive,		AccRegDotReceive,	NumNeighborReceive, ACListReceive);
#endif
  
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

	/*
	for (int i=0; i<ListSize; i++) {
		for (int dim=0; dim<Dim; dim++) {
			AccRegReceive[i][dim]    = (CUDA_REAL) AccRegReceive_f[i][dim];
			AccRegDotReceive[i][dim] = (CUDA_REAL) AccRegDotReceive_f[i][dim];
		}
	}
	*/

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
			ptcl->a_irr[dim][0] = static_cast<double>(AccRegReceive_f[i][dim]);		// Just temporarilly save new reg acc here!
			ptcl->a_irr[dim][1] = static_cast<double>(AccRegDotReceive_f[i][dim]);	// Just temporarilly save new reg acc here!
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

#ifdef CUDA_FLOAT
	delete[] AccRegReceive_f;
	delete[] AccRegDotReceive_f;
#endif 

	delete[] NumNeighborReceive;
	delete[] ACListReceive;

	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends





void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList) {

	

#ifdef CUDA_FLOAT
	// variables for saving variables to send to GPU
	CUDA_REAL * Mass;
	CUDA_REAL * Mdot;
	CUDA_REAL * Radius2;
	CUDA_REAL(*Position)[Dim];
	CUDA_REAL(*Velocity)[Dim];
	//int size = NumberOfParticle;
	int size=0, j=0;
	
	// allocate memory to the temporary variables
	Mass     = new CUDA_REAL[NumberOfParticle];
	Mdot     = new CUDA_REAL[NumberOfParticle];
	Radius2  = new CUDA_REAL[NumberOfParticle];
	Position = new CUDA_REAL[NumberOfParticle][Dim];
	Velocity = new CUDA_REAL[NumberOfParticle][Dim];

	Particle *ptcl;

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
#else
	// variables for saving variables to send to GPU
	double * Mass;
	double * Mdot;
	double * Radius2;
	double(*Position)[Dim];
	double(*Velocity)[Dim];
	//int size = NumberOfParticle;
	int size=0, j=0;


	// allocate memory to the temporary variables
	Mass     = new double[NumberOfParticle];
	Mdot     = new double[NumberOfParticle];
	Radius2  = new double[NumberOfParticle];
	Position = new double[NumberOfParticle][Dim];
	Velocity = new double[NumberOfParticle][Dim];

	Particle *ptcl;

	// copy the data of particles to the arrays to be sent
		
	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl = &particles[i];

		if (!ptcl->isActive) {
			// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
			continue;
		}

		if (RegularList.find(i) != RegularList.end()) {
			IndexList[j] = size;
			j++;
		}


		Mass[size]    = ptcl->Mass;
		Mdot[size]    = 0; //particle[i]->Mass;
		Radius2[size] = ptcl->RadiusOfNeighbor; // mass weight?

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position[size], Velocity[size]);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, Position[size], Velocity[size]);

		ActiveIndexToOriginalIndex[size] = i;
		// std::cout << "(size , i) = "  << size << " " << i << std::endl;
		size++;
	} 
#endif

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