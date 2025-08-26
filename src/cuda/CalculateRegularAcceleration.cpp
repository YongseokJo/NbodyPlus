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
#include "cuda_defs.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

void sendAllParticlesToGPU(double new_time, const int& RegularListSize, std::vector<int>& RegularListIndices);

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
	double new_time = NextRegTimeBlock*time_step;  // next regular time

	Particle *ptcl;

#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG
	std::cout << "sendAllParticlesToGPU starts" << std::endl;
#endif

#ifdef NSIGHT
	nvtxRangePushA("sendAllParticlesToGPU");
#endif
	int RegularListSize = RegularList.size();
	std::vector<int> RegularListIndices;
	sendAllParticlesToGPU(new_time, RegularListSize, RegularListIndices);  // needs to be updated
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

	CalculateAccelerationOnDevice(&ListSize, RegularListIndices); //RegularListIndices
  
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
	/*
	for (int i=0; i<ListSize; i++) {
		ptcl = &particles[RegularListIndices[i]];

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
	*/
	/*
	fprintf(stderr, "DEBUGGING STARTS...\n");
	fprintf(stderr, "ListSize: %d\n", ListSize);
	for (int i=0; i<ListSize; i++) {
		ptcl = &particles[RegularListIndices[i]];
		fprintf(stderr, "PID: %d, NumNeighbor: %d, NewNumNeighbor: %d\n", ptcl->PID, ptcl->NumberOfNeighbor, ptcl->NewNumberOfNeighbor);
		for (int j=0; j<ptcl->NewNumberOfNeighbor; j++) {
			fprintf(stderr, "%d ", ptcl->NewNeighbors[j]);
		}
		fprintf(stderr, "\n");
	}
	*/

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

	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends




// (Query MY) Let's optimize this function later. Copying data to h_ptcl in _ReceiveFromHost of cuda_my_acceleation.cpp seems super inefficient. 2025.5.24
void sendAllParticlesToGPU(double new_time, const int& RegularListSize, std::vector<int>& RegularListIndices) {

/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/

	Queue queue = {PrepareGPUCalc, -1, new_time};
	MPI_Request requests[NumberOfWorker];
	for (int i = 0; i < NumberOfWorker; i++)
		MPI_Isend(&queue, 1, QueueType, i+1, QUEUE_TAG, MPI_COMM_WORLD, &requests[i]);

	std::vector<Jparticle> Jparticles;
	Jparticles.resize(NumberOfParticle);

	std::vector<Iparticle> Iparticles;
	Iparticles.resize(RegularListSize);

	RegularListIndices.resize(RegularListSize);

	std::vector<int> counts;
	counts.resize(NumberOfProcessor * 2);
	int send_buf[2] = {0, 0};

	std::vector<int> Jcounts;
	Jcounts.resize(NumberOfProcessor);
	Jcounts[0] = 0;

	std::vector<int> Icounts;
	Icounts.resize(NumberOfProcessor);
	Icounts[0] = 0;

	std::vector<int> Jdispls;
	Jdispls.resize(NumberOfProcessor);
	Jdispls[0] = 0;

	std::vector<int> Idispls;
	Idispls.resize(NumberOfProcessor);
	Idispls[0] = 0;

	MPI_Waitall(NumberOfWorker, requests, MPI_STATUSES_IGNORE);

/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "Send job took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/
/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/

	MPI_Gather(send_buf, 2, MPI_INT, counts.data(), 2, MPI_INT, ROOT, MPI_COMM_WORLD);

	for (int rank = 1; rank < NumberOfProcessor; rank++) {
	
		Jcounts[rank] = counts[rank * 2 + 0];
		Icounts[rank] = counts[rank * 2 + 1];
	
		Jdispls[rank] = Jdispls[rank - 1] + Jcounts[rank - 1];
		Idispls[rank] = Idispls[rank - 1] + Icounts[rank - 1];
	}
	assert(Jdispls[NumberOfProcessor - 1] + Jcounts[NumberOfProcessor - 1] == NumberOfParticle);
	assert(Idispls[NumberOfProcessor - 1] + Icounts[NumberOfProcessor - 1] == RegularListSize);

	MPI_Gatherv(nullptr, 0, JparticleType, 
				Jparticles.data(), Jcounts.data(), Jdispls.data(), JparticleType, ROOT, MPI_COMM_WORLD);
	MPI_Gatherv(nullptr, 0, IparticleType, 
				Iparticles.data(), Icounts.data(), Idispls.data(), IparticleType, ROOT, MPI_COMM_WORLD);
	MPI_Gatherv(nullptr, 0, MPI_INT,
				RegularListIndices.data(), Icounts.data(), Idispls.data(), MPI_INT, ROOT, MPI_COMM_WORLD);

/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "Entire loop took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/
	// send the arrays to GPU
/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/
	SendToDevice(Jparticles, Iparticles);
/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "SendToDevice took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/

}