#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "../global.h"
#include "../QueueScheduler.h"
#include "cuda_functions.h"

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

#ifdef PerformanceTrace
	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
#endif


	// regIds are the list of positions of particles subject to regular force calculation in std::vector list particle

	// variables for opening GPU
	// const int buffer = 10;
	// int numGpuOpen = NNB+buffer;
	//const int NumPtclPerEachCalMax = 2048; // this also caps the number of particles computed each iteration
	int NeighborIndex; // this size should coincide with number of threads
	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];

	//int NumGpuCal;

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
	//int (*ACListReceive)[MaxNumNeighbor];

	//double* PotSend;
	// int **ACListReceive;
	int *ACListReceive;
	int *NumNeighborReceive;
	int MassFlag;


	double a_tmp[Dim]{0}, adot_tmp[Dim]{0};
	double da, dadot;
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;


	double DFR, FRD, SUM, AT3, BT2;
	double DTR, DTSQ, DT2, DT6,DTSQ12, DTR13;

	Particle *ptcl;


	double new_time = NextRegTimeBlock*time_step;  // next regular time


	// need to make array to send to GPU
	// allocate memory to the temporary variables
	//PotSend         = new double[ListSize];
	AccRegReceive    = new double[ListSize][Dim];
	AccRegDotReceive = new double[ListSize][Dim];
	AccIrr           = new double[ListSize][Dim];
	AccIrrDot        = new double[ListSize][Dim];

#ifdef CUDA_FLOAT
	AccRegReceive_f		= new CUDA_REAL[ListSize][Dim];
	AccRegDotReceive_f	= new CUDA_REAL[ListSize][Dim];
#endif 
	NumNeighborReceive  = new int[ListSize];

	// ACListReceive      = new int*[ListSize];
	ACListReceive = new int[ListSize * MaxNumNeighbor];

	for (int i=0; i<ListSize; i++) {
		// ACListReceive[i] = new int[MaxNumNeighbor];
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

#ifdef PerformanceTrace
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

#ifdef PerformanceTrace
	end_point_routine = std::chrono::high_resolution_clock::now();
	performance.RegularSendAllParticlesToGPU +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
	
	
	/*
	for (int i=0; i<ListSize; i++) {
		IndexList[i] = RegularList[i];
	} // endfor copy info
	*/


	//std::cout <<  "Starting Calculation On Device ..." << std::endl;
	// send information of all the particles to GPU
	// includes prediction
/*
#ifdef time_trace
	_time.reg_sendall.markStart();
#endif

	// Particles have been already at T_new through irregular time step

#ifdef time_trace
	_time.reg_sendall.markEnd();
	_time.reg_sendall.getDuration();
#endif

#ifdef time_trace
	_time.reg_gpu.markStart();
#endif
*/

#ifdef PerformanceTrace
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG
	std::cout << "CalculateAccelerationOnDevice starts" << std::endl;
#endif
  
#ifdef NSIGHT
	nvtxRangePushA("CalculateAccelerationOnDevice");
#endif
  
#ifdef CUDA_FLOAT
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive_f, AccRegDotReceive_f, NumNeighborReceive, ACListReceive);
#else
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive, AccRegDotReceive, NumNeighborReceive, ACListReceive);
#endif
  
#ifdef NSIGHT
	nvtxRangePop();
#endif
  
#ifdef DEBUG
	std::cout << "CalculateAccelerationOnDevice ended" << std::endl;
#endif

#ifdef PerformanceTrace
	end_point_routine = std::chrono::high_resolution_clock::now();
	performance.RegularGPU +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

	for (int i=0; i<ListSize; i++) {
		for (int dim=0; dim<Dim; dim++) {
			AccRegReceive[i][dim]    = (CUDA_REAL) AccRegReceive_f[i][dim];
			AccRegDotReceive[i][dim] = (CUDA_REAL) AccRegDotReceive_f[i][dim];
		}
	}
	/*
#ifdef time_trace
	_time.reg_gpu.markEnd();
	_time.reg_gpu.getDuration();

	_time.reg_cpu1.markStart();
#endif
*/



	/*
#ifdef time_trace
	_time.reg_cpu2.markEnd();
	_time.reg_cpu2.getDuration();

	_time.reg_cpu3.markStart();
#endif
*/



/*
	std::cout << "(REG_CUDA) RegularList, PID= ";
	for (int i=0; i<RegularList.size(); i++) {
		std::cout << RegularList[i]<< ", ";
	}
	std::cout << std::endl;
	*/

#ifdef PerformanceTrace
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG
	std::cout << "Adjust Regular Gravity starts" << std::endl;
#endif

#ifdef NSIGHT
	nvtxRangePushA("RegCuda");
#endif

	// Adjust Regular Gravity
	int i=0;
	TaskName task=RegCuda;
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
				MPI_Send(&task, 1, MPI_INT, (*worker)->MyRank, TASK_TAG, MPI_COMM_WORLD);
				MPI_Send(&ActiveIndexToOriginalIndex[IndexList[i]], 1, MPI_INT, (*worker)->MyRank, PTCL_TAG, MPI_COMM_WORLD);
				MPI_Send(&NumNeighborReceive[i], 1, MPI_INT, (*worker)->MyRank, 10, MPI_COMM_WORLD);
				MPI_Send(&ACListReceive[i * MaxNumNeighbor], NumNeighborReceive[i], MPI_INT, (*worker)->MyRank, 11, MPI_COMM_WORLD);
				MPI_Send(&AccRegReceive[i][0], 3, MPI_DOUBLE, (*worker)->MyRank, 12, MPI_COMM_WORLD);
				MPI_Send(&AccRegDotReceive[i][0], 3, MPI_DOUBLE, (*worker)->MyRank, 13, MPI_COMM_WORLD);
				((*worker))->onDuty = true;
				/*
				(*worker)->CurrentQueue++;
				(*worker)->CurrentQueue %= MAX_QUEUE;
				(*worker)->NumberOfQueues--;
				*/
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

#ifdef NSIGHT
	nvtxRangePop();
#endif

#ifdef DEBUG
	std::cout << "Adjust Regular Gravity ended" << std::endl;
#endif

#ifdef PerformanceTrace
	end_point_routine = std::chrono::high_resolution_clock::now();
	performance.RegularAdjust +=
		std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif


#ifdef nouse
	ws.initialize();
	int return_value, i = 0;
	ws._setTask(4);
	ws._total_tasks = RegularList.size();
	ws._completed_tasks = 0;

	do
	{
		if (ws._FreeWorkers.size() == 0 || ws._assigned_tasks == ws._total_tasks)
		{
			// have to add check all the sends are recved.
			// MPI_Waitall(NumberOfCommunication, requests, statuses);
			// NumberOfCommunication = 0;
			ws._checkCompletion(return_value);
			/* we can do something here */
			ws._Callback();
			ws._completed_tasks++;
		}
		if (ws._FreeWorkers.size() != 0 && ws._assigned_tasks < ws._total_tasks)
		{
			ws._WorkerTmp = ws._FreeWorkers.back();
			ws._FreeWorkers.pop_back();
			MPI_Send(&ws._WorkerTmp->task, 1, MPI_INT, ws._WorkerTmp->MyRank, TASK_TAG, MPI_COMM_WORLD);
			MPI_Send(&RegularList[i], 1, MPI_INT, ws._WorkerTmp->MyRank, PTCL_TAG, MPI_COMM_WORLD);
			MPI_Send(&NumNeighborReceive[i], 1, MPI_INT, ws._WorkerTmp->MyRank, 10, MPI_COMM_WORLD);
			MPI_Send(&ACListReceive[i * MaxNumNeighbor], NumNeighborReceive[i], MPI_INT, ws._WorkerTmp->MyRank, 11, MPI_COMM_WORLD);
			MPI_Send(&AccRegReceive[i][0], 3, MPI_DOUBLE, ws._WorkerTmp->MyRank, 12, MPI_COMM_WORLD);
			MPI_Send(&AccRegDotReceive[i][0], 3, MPI_DOUBLE, ws._WorkerTmp->MyRank, 13, MPI_COMM_WORLD);
			ws._WorkerTmp->onDuty = true;
			ws._assigned_tasks++;
			//fprintf(stdout, "assigned_tasks = %d/%d, number of free worker = %d pid = %d rank = %d\n",
					//ws._assigned_tasks, ws._total_tasks, ws._FreeWorkers.size(), RegularList[i], ws._WorkerTmp->MyRank);
		}
		i++;
	} while (ws._completed_tasks < ws._total_tasks);
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


#ifdef time_trace
	_time.reg_cpu3.markEnd();
	_time.reg_cpu3.getDuration();
#endif
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


