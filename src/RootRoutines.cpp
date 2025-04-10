#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cassert>
#include <mpi.h>
#include "global.h"
#include "SkipList.h"
#include "Worker.h"
#include "QueueScheduler.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

#define noDEBUG

// void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG);
// void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG);
// void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
// void InitialAssignmentOfTasks(int* data, int NumTask, int TAG);
void InitialAssignmentOfTasks(Queue queue, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void ParticleSynchronization();

bool createSkipList(SkipList *skiplist);
bool updateSkipList(SkipList *skiplist, int ptcl_id);
int writeParticle(double current_time, int outputNum);
void calculateRegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueScheduler &queue_scheduler);

void formPrimordialBinaries(int beforeLastParticleIndex);
void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls, std::unordered_map<int, int>& existing, std::unordered_map<int, int>& terminated);
void FBTermination(Particle* ptclCM);
void Merge(Particle* p1, Particle* p2);

#ifdef MULTIMAP
void createRegularMap(std::multimap<ULL,int>& RegularMap);
#ifdef CUDA
void getRegularList(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList);
void updateRegularMap(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList);
#else
void getRegularList(std::multimap<ULL,int>& RegularMap, std::vector<int>& RegularList);
void updateRegularMap(std::multimap<ULL,int>& RegularMap, std::vector<int>& RegularList);
#endif
#else // no multimap
// #ifdef CUDA
void updateNextRegTime(std::unordered_set<int>& RegularList);
// #else
// void updateNextRegTime(std::vector<int>& RegularList);
// #endif
#endif // multimap

#ifdef SEVN
void StellarEvolution();
#endif

Worker* workers;

void RootRoutines() {

	std::cout << "Root processor is ready." << std::endl;

	Particle* ptcl;
	int min_time_level=0;
	//int worker_rank;
	TaskName task;
	int total_tasks;
	int remaining_tasks=0, completed_tasks=0, completed_rank;

	std::unordered_map<int, int> CMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11
	std::unordered_map<int, int> PrevCMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11
	std::vector<int> newCMptcls; // by EW 2025.1.6 // unordered_set? by EW 2025.1.11
	std::vector<int> EmptyIndex; // by EW 2025.1.7  empty slots in particles e.g., due to mergers
	// unordered_set? by EW 2025.1.11
	// merged particles & PISN will be contained here
	// new single Particle formed in Enzo can be formed in ParticleIndex of these ptcls
	// if empty, LastParticleIndex++

	bool bin_termination = false;
	bool new_binaries = false;
#ifdef MULTIMAP
	std::multimap<ULL,int> RegularMap;
#endif

// #ifdef CUDA
	std::unordered_set<int> RegularList;
// #else
	// std::vector<int> RegularList;
// #endif

	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object

	int ptcl_id;

	workers = new Worker[NumberOfWorker+1];

	for (int i=0; i<=NumberOfWorker; i++) {
		workers[i].initialize(i);
	}

	QueueScheduler queue_scheduler;

#ifdef PERFORMANCETRACE
	std::chrono::high_resolution_clock::time_point start_point_whole;
	std::chrono::high_resolution_clock::time_point end_point_whole;

	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
#ifdef MULTIMAP
	std::chrono::high_resolution_clock::time_point start_point_map;
	std::chrono::high_resolution_clock::time_point end_point_map;
#endif
#endif // performance

	/* Initialization */
	{
		std::vector<int> PIDs;
		PIDs.reserve(NumberOfParticle);
		PIDs.resize(NumberOfParticle);
		for (int i = 0; i <= LastParticleIndex; i++)
		{
			PIDs[i] = i;
		}

		std::cout << "Initialization of particles starts." << std::endl;
		queue_scheduler.initialize(InitAcc1);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			//queue_scheduler.printFreeWorker();
			//queue_scheduler.printWorkerToGo();
			queue_scheduler.assignQueueAuto();
			//queue_scheduler.printFreeWorker();
			//queue_scheduler.printWorkerToGo();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());
		std::cout << "Init 01 done" << std::endl;

		queue_scheduler.initialize(InitAcc2);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			queue_scheduler.assignQueueAuto();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());
		std::cout << "Init 02 done" << std::endl;

#ifdef FEWBODY
		// Primordial binary search
		queue_scheduler.initialize(SearchPrimordialGroup);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			queue_scheduler.assignQueueAuto();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());
		std::cout << "Primordial binary search done" << std::endl;

		Queue queue;
		int rank;
		int OriginalLastParticleIndex = LastParticleIndex;
		formPrimordialBinaries(OriginalLastParticleIndex);
		assert(OriginalLastParticleIndex <= LastParticleIndex); // for debugging by EW 2025.1.4
		assert(CMPtclWorker.empty()); // for debugging by EW 2025.1.4
		if (OriginalLastParticleIndex != LastParticleIndex) {
			std::cout << "In total, " << LastParticleIndex - OriginalLastParticleIndex
					  << " primordial binaries are created." << std::endl;
			queue_scheduler.initialize(MakePrimordialGroup);
			for (int i=OriginalLastParticleIndex+1; i<=LastParticleIndex; i++) {
				std::cout << "New Primordial Binary of PID="
						  << i << " is created with being assigned to a worker of rank "
						  << rank << "." << std::endl;
				ptcl = &particles[i];
				CMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker.size() % NumberOfWorker + 1});
				PIDs.push_back(ptcl->ParticleIndex);
				rank = CMPtclWorker[ptcl->ParticleIndex];

				queue.task = MakePrimordialGroup;
				queue.pid = ptcl->ParticleIndex;
				workers[rank].addQueue(queue);
				queue_scheduler.assignWorker(&workers[rank]);
			}
			queue_scheduler.setTotalQueue(CMPtclWorker.size());
			do {
				queue_scheduler.runQueueAuto();
				queue_scheduler.waitQueue(0);
			} while(queue_scheduler.isComplete());
		}
		else {
			std::cout << "There is no primordial binary." << std::endl;
		}
		fprintf(stdout, "PrimordialBinariesRoutine has ended...\n"
						"The total number of particles is  %d\n",
				NumberOfParticle);
		fflush(stdout);
#endif

		// Initialize Time Step
		queue_scheduler.initialize(InitTime);
		queue_scheduler.takeQueue(PIDs);
		do
		{
			queue_scheduler.assignQueueAuto();
			queue_scheduler.runQueueAuto();
			queue_scheduler.waitQueue(0); //blocking wait
		} while(queue_scheduler.isComplete());

		/*
		for (int i=0; i<=LastParticleIndex; i++) {
			ptcl = &particles[i];
			if (ptcl->isActive)
				fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
								"NumNeighbor= %d\n",
						ptcl->PID,
						ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockIrr,
						ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockReg,
						// NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
						// NextRegTimeBlock,
						ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeBlockIrr,
						ptcl->TimeLevelIrr,
						ptcl->TimeBlockReg,
						ptcl->TimeLevelReg,
						ptcl->NumberOfNeighbor);
		}
		*/
	} // Initialization ends



	/* synchronization */
	//ParticleSynchronization();

	/* timestep correction */
	{
		std::cout << "Time Step correction." << std::endl;
		for (int i=0; i<=LastParticleIndex; i++) {
			ptcl = &particles[i];

			if (!ptcl->isActive)
				continue;

			if (ptcl->NumberOfNeighbor != 0) {
				while (ptcl->TimeLevelIrr >= ptcl->TimeLevelReg) {
					ptcl->TimeStepIrr *= 0.5;
					ptcl->TimeBlockIrr *= 0.5;
					ptcl->TimeLevelIrr--;
				}
			}
			if (ptcl->TimeLevelIrr < min_time_level) {
				min_time_level = ptcl->TimeLevelIrr;
			}
		}

		// resetting time_block based on the system
		time_block = std::max(-60, min_time_level-MIN_LEVEL_BUFFER);
		block_max = static_cast<ULL>(pow(2, -time_block));
		time_step = pow(2,time_block);

		for (int i=0; i<=LastParticleIndex; i++) {
			ptcl = &particles[i];

			if (!ptcl->isActive)
				continue;

			ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr-time_block));
			ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block));
#ifdef IRR_TEST
			ptcl->TimeStepReg = 1;
			ptcl->TimeLevelReg = 0;
			ptcl->TimeBlockReg = block_max;
#endif
			ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
		}
		std::cout << "Time Step done." << std::endl;
	}

	/*
	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl = &particles[i];
		if (ptcl->isActive)
			fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
							"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
							"NumNeighbor= %d\n",
					ptcl->PID,
					ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
					ptcl->CurrentBlockIrr,
					ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
					ptcl->CurrentBlockReg,
					// NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
					// NextRegTimeBlock,
					ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
					ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
					ptcl->TimeBlockIrr,
					ptcl->TimeLevelIrr,
					ptcl->TimeBlockReg,
					ptcl->TimeLevelReg,
					ptcl->NumberOfNeighbor);
	}
	*/

	/* timestep variable synchronization */
	{
		std::cout << "Time Step synchronization." << std::endl;
		task=TimeSync;
		completed_tasks = 0; total_tasks = NumberOfWorker;
		Queue queue = {task, -1, -1.0};
		InitialAssignmentOfTasks(queue, NumberOfWorker, QUEUE_TAG);
		//MPI_Waitall(NumberOfCommunication, requests, statuses);
		//NumberOfCommunication = 0;
		broadcastFromRoot(time_block);
		broadcastFromRoot(block_max);
		broadcastFromRoot(time_step);
		fprintf(stdout, "TimeSync broadcasd done.\n");
		fflush(stdout);
		//MPI_Win_sync(win);  // Synchronize memory
		//MPI_Barrier(shared_comm);
		while (completed_tasks < total_tasks) {
			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			completed_tasks++;
		}
		fprintf(stderr, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);
		fflush(stderr);
	}

	/* Particle Initialization Check */
	// /*
	{
		for (int i=0; i<=LastParticleIndex; i++) {
			ptcl = &particles[i];
			if (ptcl->isActive)
				fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
								"NumNeighbor= %d\n",
						ptcl->PID,
						ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockIrr,
						ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->CurrentBlockReg,
						// NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
						// NextRegTimeBlock,
						ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
						ptcl->TimeBlockIrr,
						ptcl->TimeLevelIrr,
						ptcl->TimeBlockReg,
						ptcl->TimeLevelReg,
						ptcl->NumberOfNeighbor);
				/*
				fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
					a1_reg = (%.4e,%.4e,%.4e), a2_reg = (%.4e,%.4e,%.4e), a3_reg = (%.4e,%.4e,%.4e)\n\
					a1_irr = (%.4e,%.4e,%.4e), a2_irr = (%.4e,%.4e,%.4e), a3_irr = (%.4e,%.4e,%.4e)\n", 
					ptcl->a_tot[0][0],
					ptcl->a_tot[1][0],
					ptcl->a_tot[2][0],
					ptcl->a_reg[0][0],
					ptcl->a_reg[1][0],
					ptcl->a_reg[2][0],
					ptcl->a_irr[0][0],
					ptcl->a_irr[1][0],
					ptcl->a_irr[2][0],
					ptcl->NumberOfNeighbor,
					ptcl->RadiusOfNeighbor,
					ptcl->a_reg[0][1],
					ptcl->a_reg[1][1],
					ptcl->a_reg[2][1],
					ptcl->a_reg[0][2],
					ptcl->a_reg[1][2],
					ptcl->a_reg[2][2],
					ptcl->a_reg[0][3],
					ptcl->a_reg[1][3],
					ptcl->a_reg[2][3],
					ptcl->a_irr[0][1],
					ptcl->a_irr[1][1],
					ptcl->a_irr[2][1],
					ptcl->a_irr[0][2],
					ptcl->a_irr[1][2],
					ptcl->a_irr[2][2],
					ptcl->a_irr[0][3],
					ptcl->a_irr[1][3],
					ptcl->a_irr[2][3]);
				*/
		}
	} // Particle Initialization Check ends
	// */


	/* Actual Loop */
	{
		int max_level = 5;
		double prob = 0.5;
		SkipList *skiplist;
		Node* ThisLevelNode;
		double current_time_irr=0;
		double next_time=0;
		int ptcl_id_return;
		Worker* worker;
		
#ifdef MULTIMAP
		createRegularMap(RegularMap);
#endif

		//ParticleSynchronization();
		while (1) {

			// create output at appropriate time intervals
			if (global_time >= outputTime) {
				writeParticle(global_time, outNum++);
				outputTime += outputTimeStep;
			}

			// end if the global time exceeds the end time
			if (global_time >= 1) {
				task=Ends;
				Queue queue = {task, -1, -1.0};
				InitialAssignmentOfTasks(queue, NumberOfWorker, QUEUE_TAG);
				MPI_Type_free(&QueueType);
				//MPI_Waitall(NumberOfCommunication, requests, statuses);
				//NumberOfCommunication = 0;
				std::cout << EnzoTimeStep << std::endl;
				std::cout << "Simulation Done!" << std::endl;
				return;
			}

#ifdef PERFORMANCETRACE
			start_point_whole = std::chrono::high_resolution_clock::now();
#endif

#ifndef MULTIMAP

#ifdef PERFORMANCETRACE
			start_point_routine = std::chrono::high_resolution_clock::now();
#endif
#ifdef NSIGHT
			nvtxRangePushA("updateNextRegTime");
#endif
			updateNextRegTime(RegularList);
#ifdef NSIGHT
			nvtxRangePop();
#endif
#ifdef PERFORMANCETRACE
			end_point_routine = std::chrono::high_resolution_clock::now();
			performance.UpdateNextRegTime +=
				std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#endif // no multimap

			bin_termination = false;
			new_binaries = false;
				
			/*
			std::cout << "NextRegTimeBlock=" << NextRegTimeBlock << std::endl;
			std::cout << "PID= ";
			for (int i : RegularList)
				std::cout << i<< ", ";
			std::cout << std::endl;
			std::cout << "size of regularlist= " << RegularList.size() << std::endl;
			*/
#ifdef PERFORMANCETRACE
			start_point_routine = std::chrono::high_resolution_clock::now();
#endif		
#ifdef NSIGHT
			nvtxRangePushA("createSkipList");
#endif
			skiplist = new SkipList(max_level, prob);
			if (createSkipList(skiplist) == FAIL)
				fprintf(stderr, "There are no irregular particles!\nBut is it really happening? check skiplist->display()\n");
#ifdef NSIGHT
			nvtxRangePop();
#endif
#ifdef PERFORMANCETRACE
			end_point_routine = std::chrono::high_resolution_clock::now();
			performance.SkipListCreate +=
				std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif


			// Irregular
			while ( skiplist->getFirstNode() != nullptr) {

				ThisLevelNode = skiplist->getFirstNode();
				ThisLevelNode->ParticleList.erase(
					std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
						[](int i) {
						return !particles[i].isActive;
						}
					),
					ThisLevelNode->ParticleList.end()
				);
				if (ThisLevelNode->ParticleList.size() == 0) {
					skiplist->deleteFirstNode();
					continue;
				}
				// fprintf(stdout, "N: %d\n", ThisLevelNode->ParticleList.size());
				next_time     = particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr\
									 	    + particles[ThisLevelNode->ParticleList[0]].TimeStepIrr;

#ifdef PERFORMANCETRACE
				start_point_routine = std::chrono::high_resolution_clock::now();
#endif
			
#ifdef DEBUG
				// print out particlelist
				fprintf(stdout, "(IRR_FORCE) next_time: %e Myr\n", next_time*EnzoTimeStep*1e4);
				// /*
				fprintf(stdout, "PID: %d. CurrentTimeIrr: %e Myr, TimeStepIrr: %e Myr\n", 
							particles[ThisLevelNode->ParticleList[0]].PID, 
							particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr*EnzoTimeStep*1e4, 
							particles[ThisLevelNode->ParticleList[0]].TimeStepIrr*EnzoTimeStep*1e4);

				// fprintf(stdout, "PID (%d) = ", ThisLevelNode->ParticleList.size());
				for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
					ptcl = &particles[ThisLevelNode->ParticleList[i]];
					// fprintf(stdout, "%d, ", ptcl->PID);
					fprintf(stdout, "PID: %d. %e Myr, %e Myr\n", 
							ptcl->PID,
							ptcl->CurrentTimeIrr*EnzoTimeStep*1e4,
							ptcl->TimeStepIrr*EnzoTimeStep*1e4);
				}
				fprintf(stdout, "\n");
				// fflush(stdout);
				// */
#endif

				// Irregular Force
#ifdef FEWBODY
#ifdef NSIGHT
				nvtxRangePushA("IrregularForce");
#endif
#ifdef DEBUG
				std::cout << "Irr force starts" << std::endl;
#endif
// /*
				int cm_pid;
				Queue queue;
				queue_scheduler.initializeIrr(IrrForce, next_time, ThisLevelNode->ParticleList);
				auto iter = queue_scheduler.CMPtcls.begin();
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					do {
						worker = queue_scheduler.waitQueue(1); // non-blocking wait
						// if there's any CMPtcl
						if (queue_scheduler.CMPtcls.size() > 0)
						{
							// check if there's any CM ptcl ready to go for SDAR
							if (iter == queue_scheduler.CMPtcls.end())
								iter = queue_scheduler.CMPtcls.begin();
							cm_pid = *(iter);
							ptcl = &particles[cm_pid];
							/*
							for (int j = 0; j < ptcl->NumberOfNeighbor; j++)
							{
								// if (particles[ptcl->Neighbors[j]].isUpdateToDate == false) // original code
								if (particles[ptcl->Neighbors[j]].isActive && !particles[ptcl->Neighbors[j]].isUpdateToDate) // modified by EW 2025.2.26
								{
									iter++;
									goto skip_to_next;
								}
							}
							*/
							queue.task = ARIntegration;
							queue.pid = cm_pid;
							queue.next_time = next_time;
							workers[CMPtclWorker[cm_pid]].addQueue(queue);
							queue_scheduler.assignWorker(&workers[CMPtclWorker[cm_pid]]);
							iter = queue_scheduler.CMPtcls.erase(iter);
						// skip_to_next:;
						}
						if (worker != nullptr) 
						{

						}
					} while (worker == nullptr);
					queue_scheduler.callback(worker);
				} while (queue_scheduler.isComplete());
// */
/*
				queue_scheduler.initialize(IrrForce, next_time);
				queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());

				Queue queue;
				for (int ptcl_id : ThisLevelNode->ParticleList)
				{
					ptcl = &particles[ptcl_id];

					if (ptcl->isCMptcl) {
						int rank = CMPtclWorker[ptcl->ParticleIndex];
						queue.task = ARIntegration;
						queue.pid = ptcl->ParticleIndex;
						queue.next_time = next_time;
						workers[rank].addQueue(queue);
						workers[rank].runQueue();
						workers[rank].callback();
					}
				}
*/
#ifdef DEBUG
				 std::cout << "Irregular Force done" << std::endl;
#endif
#ifdef NSIGHT
				nvtxRangePop();
#endif
#else
				queue_scheduler.initialize(IrrForce, next_time);
				queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
#endif

#ifdef PERFORMANCETRACE
				end_point_routine = std::chrono::high_resolution_clock::now();
				performance.IrregularForce +=
					std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

				//ParticleSynchronization();


#ifdef PERFORMANCETRACE
				start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("IrregularUpdate");
#endif
/*
				// Irregular Update
				queue_scheduler.initialize(IrrUpdate);
				queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
*/
				for (int ptcl_id : ThisLevelNode->ParticleList)
				{
					ptcl = &particles[ptcl_id];

					if (ptcl->NumberOfNeighbor != 0) // IAR modified
						ptcl->updateParticle();
					ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
					ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
				}
#ifdef DEBUG
				for (int i: ThisLevelNode->ParticleList) {
					ptcl = &particles[i];
					if (ptcl->CurrentTimeIrr != next_time) {
						fprintf(stdout, "Error! PID: %d, CurrentTimeIrr: %e Myr, next_time: %e Myr\n", ptcl->PID, ptcl->CurrentTimeIrr*EnzoTimeStep*1e4, next_time*EnzoTimeStep*1e4);
						assert(ptcl->CurrentTimeIrr == next_time);
					}
				}
				 std::cout << "Irregular update done" << std::endl;
#endif
#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.IrregularUpdate +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#ifdef FEWBODY

#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("FewBodyTermination");
#endif
				int OriginalSize = ThisLevelNode->ParticleList.size();
				for (int i=0; i<OriginalSize; i++ ){
					ptcl = &particles[ThisLevelNode->ParticleList[i]];
					if (ptcl->getBinaryInterruptState() == BinaryInterruptState::merger ||
						ptcl->getBinaryInterruptState() == BinaryInterruptState::terminated) {

						assert(ptcl->isCMptcl); // for debugging by EW 2025.1.20
					
						if (ptcl->getBinaryInterruptState() == BinaryInterruptState::merger) {
							if (ptcl->NumberOfMember == 2) { // binary merger

								Particle* donor = &particles[ptcl->Members[0]];
								Particle* accretor = &particles[ptcl->Members[1]];

								Merge(donor, accretor);
#ifdef SEVN
								if (donor->StellarEvolution != nullptr && accretor->StellarEvolution != nullptr) {
									fprintf(stdout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->PID, donor->Mass*mass_unit, donor->StellarEvolution->get_zams());
									fprintf(stdout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->PID, accretor->Mass*mass_unit, accretor->StellarEvolution->get_zams());
									fprintf(stdout, "Donor: amiempty(): %d\n", donor->StellarEvolution->amiempty());
									fprintf(stdout, "Accretor: amiempty(): %d\n", accretor->StellarEvolution->amiempty());
									fflush(stdout);
								}
#endif
							}
							else { // from NewFBInitialization3

								assert(ptcl->NumberOfMember > 2); // for debugging by EW 2025.1.20

								Particle* donor;
								Particle* accretor;

								for (int j=0; j<ptcl->NumberOfMember; j++) {
									if (particles[ptcl->Members[j]].getBinaryInterruptState() == BinaryInterruptState::collision) {
										donor = &particles[ptcl->Members[j]];
										accretor = &particles[donor->getBinaryPairID()];

										assert(accretor->getBinaryInterruptState() == BinaryInterruptState::collision);
										assert(accretor->getBinaryPairID() == donor->ParticleIndex);
										break;
									}
								}

								Merge(donor, accretor);
#ifdef SEVN
								if (donor->StellarEvolution != nullptr && accretor->StellarEvolution != nullptr) {
									fprintf(stdout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->PID, donor->Mass*mass_unit, donor->StellarEvolution->get_zams());
									fprintf(stdout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->PID, accretor->Mass*mass_unit, accretor->StellarEvolution->get_zams());
									fprintf(stdout, "Donor: amiempty(): %d\n", donor->StellarEvolution->amiempty());
									fprintf(stdout, "Accretor: amiempty(): %d\n", accretor->StellarEvolution->amiempty());
									fflush(stdout);
								}
#endif
								int rank = CMPtclWorker[ptcl->ParticleIndex];
								queue.task = MergeManyBody;
								queue.pid = ptcl->ParticleIndex;
								workers[rank].addQueue(queue);
								workers[rank].runQueue();
								workers[rank].callback();
								
								continue;
							}
						}

						bin_termination = true;
						ptcl->isActive = false;

						if (ptcl->ParticleIndex == LastParticleIndex) {
							LastParticleIndex--;
							global_variable->LastParticleIndex == LastParticleIndex;
						}
						else
							PrevCMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker[ptcl->ParticleIndex]});
						CMPtclWorker.erase(ptcl->ParticleIndex);

						for (int j=0; j < ptcl->NumberOfMember; j++) {
							particles[ptcl->Members[j]].CMPtclIndex = -1;
							if (particles[ptcl->Members[j]].Mass < 0.0) {
#ifdef SEVN
								Particle* ptcl_erased = &particles[ptcl->Members[j]];
								fprintf(stdout, "ptcl_erased... PID: %d\n", ptcl_erased->PID);
								if (ptcl_erased->StellarEvolution != nullptr) {

									auto it = SEVNList.begin();
									while (it != SEVNList.end()) {
										if (it->second == ptcl_erased->ParticleIndex) {
											it = SEVNList.erase(it);
											fprintf(stdout, "Merger induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl_erased->PID);
											break;
										}
										else
											it++;
									}

									delete ptcl_erased->StellarEvolution;
									ptcl_erased->StellarEvolution = nullptr;
									fprintf(stdout, "Merger induced zero mass particle (PID: %d) SEVN memory is free now\n", ptcl_erased->PID);
								}
								fflush(stdout);
#endif
								continue;
							}
							ThisLevelNode->ParticleList.push_back(ptcl->Members[j]);
							particles[ptcl->Members[j]].isActive = true;
						}
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
						start_point_map = std::chrono::high_resolution_clock::now();
#endif
						auto range = RegularMap.equal_range(ptcl->CurrentBlockReg + ptcl->TimeBlockReg);
						for (auto it = range.first; it != range.second; ++it) {
							if (ptcl->PID == particles[it->second].PID) {
								RegularMap.erase(it);
								break;
							}
						}
#ifdef PERFORMANCETRACE
						end_point_map = std::chrono::high_resolution_clock::now();
						performance.RegularMap +=
							std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#endif // multimap
						FBTermination(ptcl);
					}
				}

				if (bin_termination) {
					for (int i=OriginalSize; i<ThisLevelNode->ParticleList.size(); i++) {
						ptcl = &particles[ThisLevelNode->ParticleList[i]];
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
						start_point_map = std::chrono::high_resolution_clock::now();
#endif

						RegularMap.insert({ptcl->CurrentBlockReg + ptcl->TimeBlockReg, ptcl->ParticleIndex});

#ifdef PERFORMANCETRACE
						end_point_map = std::chrono::high_resolution_clock::now();
						performance.RegularMap +=
							std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#else // no multimap
						if (ptcl->CurrentBlockReg + ptcl->TimeBlockReg == NextRegTimeBlock) {
// #ifdef CUDA
							RegularList.insert(ptcl->ParticleIndex);
// #else
							// RegularList.push_back(ptcl->ParticleIndex);
// #endif
						}
#endif // multimap

						ptcl->NewNumberOfNeighbor = 0;
						if (ptcl->TimeStepIrr * EnzoTimeStep * 1e4 < TSEARCH)
							ptcl->checkNewGroup4();
					}

					// Erase terminated CM particles by EW 2025.1.6
					ThisLevelNode->ParticleList.erase(
						std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
							[](int i) {
							return !particles[i].isActive;
							}
						),
						ThisLevelNode->ParticleList.end()
					);
				}

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.FewBodyTermination +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("FewBodySearch");
#endif

#ifdef DEBUG
				std::cout << "FB search starts" << std::endl;
#endif
/*
				// std::cerr << "FB search starts" << std::endl;
				// Few-body group search
				queue_scheduler.initialize(SearchGroup);
				queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
				do
				{
					queue_scheduler.assignQueueAuto();
					queue_scheduler.runQueueAuto();
					// queue_scheduler.printStatus();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());

				// std::cerr << "FB search ended" << std::endl;
*/
				for (int ptcl_id : ThisLevelNode->ParticleList)
				{
					ptcl = &particles[ptcl_id];
					
					if (ptcl->getBinaryInterruptState() == BinaryInterruptState::threebody) {
						ptcl->setBinaryInterruptState(BinaryInterruptState::none);
					}
					else if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody)
					{
						ptcl->NewNumberOfNeighbor = 0;
						ptcl->checkNewGroup2();
						ptcl->setBinaryInterruptState(BinaryInterruptState::none);
					}
					else if (ptcl->NewNumberOfNeighbor != 0)
					{	
						ptcl->checkNewGroup3();
					}
				}
#ifdef DEBUG
				std::cout << "FB search ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.FewBodySearch +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("FormBinaries");
#endif
				int OriginalParticleListSize = ThisLevelNode->ParticleList.size();
				int rank_delete, rank_new;
#ifdef DEBUG
				std::cout << "formBinaries starts" << std::endl;
#endif
				formBinaries(ThisLevelNode->ParticleList, newCMptcls, CMPtclWorker, PrevCMPtclWorker);
#ifdef DEBUG
				std::cout << "formBinaries ended" << std::endl;
#endif
				if (OriginalParticleListSize != ThisLevelNode->ParticleList.size()) {
#ifdef DEBUG
					std::cout << "New Binary!" << std::endl;
#endif
					new_binaries = true;

					// new code by EW 2025.1.26
					Particle* ptclCM;
					Particle* mem_ptclCM;
					for (int i=0; i<newCMptcls.size(); i++) {
						ptclCM = &particles[newCMptcls[i]]; // 2025.01.10 edited to newCMptcls[i] by YS

						for (int j=0; j<ptclCM->NewNumberOfNeighbor; j++) {
							mem_ptclCM = &particles[ptclCM->NewNeighbors[j]];
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
							start_point_map = std::chrono::high_resolution_clock::now();
#endif

							auto range = RegularMap.equal_range(mem_ptclCM->CurrentBlockReg + mem_ptclCM->TimeBlockReg);
							for (auto it = range.first; it != range.second; ++it) {
								if (mem_ptclCM->PID == particles[it->second].PID) {
									RegularMap.erase(it);
									break;
								}
							}
#ifdef PERFORMANCETRACE
							end_point_map = std::chrono::high_resolution_clock::now();
							performance.RegularMap +=
								std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#endif // multimap
							if (mem_ptclCM->isCMptcl) {
								fprintf(stdout, "manybody group detected; PID %d should be deleted first\n", mem_ptclCM->PID);
								rank_delete = CMPtclWorker[mem_ptclCM->ParticleIndex];
								fprintf(stdout, "Rank of CM ptcl %d: %d\n", mem_ptclCM->PID, rank_delete);
								queue.task = DeleteGroup;
								queue.pid = mem_ptclCM->ParticleIndex;
								workers[rank_delete].addQueue(queue);
								workers[rank_delete].runQueue();
								workers[rank_delete].callback();

								PrevCMPtclWorker.insert({mem_ptclCM->ParticleIndex, CMPtclWorker[mem_ptclCM->ParticleIndex]});
								CMPtclWorker.erase(mem_ptclCM->ParticleIndex);
							}
						}

						rank_new = CMPtclWorker[ptclCM->ParticleIndex];
#ifdef DEBUG
						fprintf(stdout, "Rank of CM ptcl %d: %d\n", ptclCM->PID, rank_new);
#endif
						queue.task = MakeGroup;
						queue.pid = ptclCM->ParticleIndex;
						workers[rank_new].addQueue(queue);
						workers[rank_new].runQueue();
						workers[rank_new].callback();
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
						start_point_map = std::chrono::high_resolution_clock::now();
#endif
						RegularMap.insert({ptclCM->CurrentBlockReg + ptclCM->TimeBlockReg, ptclCM->ParticleIndex});
#ifdef PERFORMANCETRACE
						end_point_map = std::chrono::high_resolution_clock::now();
						performance.RegularMap +=
							std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#else // no multimap
						if (ptclCM->CurrentBlockReg + ptclCM->TimeBlockReg == NextRegTimeBlock) {
// #ifdef CUDA
							RegularList.insert(ptcl->ParticleIndex);
// #else
							// RegularList.push_back(ptcl->ParticleIndex);
// #endif
						}
#endif
					}
#ifdef DEBUG
					std::cout << "All new fewbody objects are initialized." << std::endl;
#endif

					ThisLevelNode->ParticleList.erase(
						std::remove_if(
								ThisLevelNode->ParticleList.begin(), 
								ThisLevelNode->ParticleList.end(),
								[](int i) { return !particles[i].isActive; }
						),
						ThisLevelNode->ParticleList.end()
					);
				}
				newCMptcls.clear();

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.FewBodyInitialization +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#endif

#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif
#ifdef DEBUG
				std::cout << "updateSkipList starts" << std::endl;
#endif
				for (int i=0; i<ThisLevelNode->ParticleList.size(); i++)
					updateSkipList(skiplist, ThisLevelNode->ParticleList[i]);
#ifdef DEBUG
				std::cout << "updateSkipList ended" << std::endl;
#endif
#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.SkipListUpdate +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

				current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr*time_step;
#ifdef DEBUG
				std::cout << "skiplist->deleteFirstNode() starts" << std::endl;
#endif
				skiplist->deleteFirstNode();
#ifdef DEBUG
				std::cout << "skiplist->deleteFirstNode() ended" << std::endl;
#endif

#ifdef IRR_TEST
				std::cout << "current_time_irr=" << current_time_irr<< std::endl;
				// create output at appropriate time intervals
				if (current_time_irr >= outputTime) {
					writeParticle(current_time_irr, outNum++);
					outputTime += outputTimeStep;
				}

				// end if the global time exceeds the end time
				if (current_time_irr >= 1) {
					task=-100;
					Queue queue = {task, -1, -1.0};
					InitialAssignmentOfTasks(Queue, NumberOfWorker, QUEUE_TAG);
					MPI_Waitall(NumberOfCommunication, requests, statuses);
					NumberOfCommunication = 0;
					std::cout << EnzoTimeStep << std::endl;
					std::cout << "Simulation Done!" << std::endl;
					return;
				}
#endif
			} // Irr
#ifdef DEBUG
			std::cout << "delete skiplist" << std::endl;
#endif
			delete skiplist;
			skiplist = nullptr;
			//exit(SUCCESS);

#ifdef FEWBODY
			if (bin_termination || new_binaries) {
#ifdef MULTIMAP
				if (NextRegTimeBlock != RegularMap.begin()->first) {
					NextRegTimeBlock = RegularMap.begin()->first;
					global_variable->NextRegTimeBlock = NextRegTimeBlock;
					continue;
				}
#else // no multimap
				for (auto it = RegularList.begin(); it != RegularList.end(); ) {
					if (!particles[*it].isActive)
						it = RegularList.erase(it);
					else
						++it;
				}

				if (RegularList.empty())
					continue;
#endif
			}
#endif

#ifdef CUDA
			{
				next_time = NextRegTimeBlock*time_step;
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif
				getRegularList(RegularMap, RegularList);

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
#endif // multimap

#ifdef NSIGHT
				nvtxRangePushA("calculateRegAccelerationOnGPU");
#endif

#ifdef DEBUG
				std::cout << "calculateRegAccelerationOnGPU starts" << std::endl;
				std::cout << "RegularList size: " << RegularList.size() << std::endl;
#endif

				calculateRegAccelerationOnGPU(RegularList, queue_scheduler);

#ifdef DEBUG
				std::cout << "calculateRegAccelerationOnGPU ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("RegCudaUpdate");
#endif

#ifdef DEBUG
				std::cout << "update regular starts" << std::endl;
#endif
				// Update Regular
				queue_scheduler.initialize(RegCudaUpdate);
				queue_scheduler.takeQueueRegularList(RegularList);
				do
				{
					queue_scheduler.assignQueueAutoRegularList();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
#ifdef DEBUG
				std::cout << "update regular ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularUpdate +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif
				updateRegularMap(RegularMap, RegularList);

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
#endif // multimap
			}

#else // no cuda
			{
				next_time = NextRegTimeBlock*time_step;
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
				start_point_routine = std::chrono::high_resolution_clock::now();
#endif
				getRegularList(RegularMap, RegularList);

#ifdef PERFORMANCETRACE
				end_point_routine = std::chrono::high_resolution_clock::now();
				performance.RegularMap +=
					std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
#endif // multimap

#ifdef PERFORMANCETRACE
				start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("RegForce");
#endif

#ifdef DEBUG
				std::cout << "Regular force starts" << std::endl;
#endif
// /*
				// Regular force
				queue_scheduler.initialize(RegForce);
				queue_scheduler.takeQueueRegularList(RegularList);
				do
				{
					queue_scheduler.assignQueueAutoRegularList();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
// */
/*
				// Regular Gravity
				task = RegForce;
				completed_tasks = 0;
				total_tasks = RegularList.size();
				next_time = NextRegTimeBlock*time_step;

				//std::cout << "TotalTask=" << total_tasks << std::endl;

				//std::cout << std::endl;

				InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
				InitialAssignmentOfTasks(RegularList, next_time, total_tasks, PTCL_TAG);
				// MPI_Waitall(NumberOfCommunication, requests, statuses);
				// NumberOfCommunication = 0;

				// further assignments
				remaining_tasks = total_tasks-NumberOfWorker;
				while (completed_tasks < total_tasks) {
					// Check which worker is done
					MPI_Irecv(&ptcl_id_return, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
					completed_rank = status.MPI_SOURCE;

					if (remaining_tasks > 0) {
						ptcl_id = RegularList[NumberOfWorker + completed_tasks];
						MPI_Send(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD);
						MPI_Send(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD);
						MPI_Send(&next_time, 1, MPI_DOUBLE, completed_rank, TIME_TAG, MPI_COMM_WORLD);
						remaining_tasks--;
					} else {
						//printf("Rank %d: No more tasks to assign\n", completed_rank);
					}
					//updateSkipList(skiplist, ptcl_id_return);
					completed_tasks++;
				}
*/
#ifdef DEBUG
				std::cout << "Regular force ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularForce +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif


#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("RegUpdate");
#endif

#ifdef DEBUG
				std::cout << "update regular starts" << std::endl;
#endif
// /*
				// Update Regular
				queue_scheduler.initialize(RegUpdate);
				queue_scheduler.takeQueueRegularList(RegularList);
				do
				{
					queue_scheduler.assignQueueAutoRegularList();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
// */
/*
				// Regular Update
				//std::cout<< "Reg Acc Done." <<std::endl;
				task = RegUpdate;
				completed_tasks = 0;

				InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
				InitialAssignmentOfTasks(RegularList, total_tasks, PTCL_TAG);
				// MPI_Waitall(NumberOfCommunication, requests, statuses);
				// NumberOfCommunication = 0;

				// further assignments
				remaining_tasks = total_tasks-NumberOfWorker;
				while (completed_tasks < total_tasks) {
					// Check which worker is done
					MPI_Irecv(&ptcl_id_return, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
					completed_rank = status.MPI_SOURCE;

					if (remaining_tasks > 0) {
						ptcl_id = RegularList[NumberOfWorker + completed_tasks];
						//MPI_Isend(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD, &request);
						//MPI_Isend(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD, &request);
						MPI_Send(&task,      1, MPI_INT, completed_rank, TASK_TAG, MPI_COMM_WORLD);
						MPI_Send(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, MPI_COMM_WORLD);
						remaining_tasks--;
					} else {
						//printf("Rank %d: No more tasks to assign\n", completed_rank);
					}
					completed_tasks++;
				}
*/
#ifdef DEBUG
				std::cout << "update regular ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularUpdate +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif
				updateRegularMap(RegularMap, RegularList);

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
#endif // multimap
			} // no cuda regular routine ends
#endif // cuda

			global_time = NextRegTimeBlock*time_step;

#ifdef SEVN

#ifdef PERFORMANCETRACE
			start_point_routine = std::chrono::high_resolution_clock::now();
#endif

			if (!SEVNList.empty() && SEVNList.begin()->first <= global_time*EnzoTimeStep*1e4)
				StellarEvolution(); // Currently, evolving all the particles upto global_time

#ifdef PERFORMANCETRACE
			end_point_routine = std::chrono::high_resolution_clock::now();
			performance.StellarEvolution +=
				std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#endif

#ifdef PERFORMANCETRACE
			end_point_whole = std::chrono::high_resolution_clock::now();
			performance.WholeRoutine +=
				std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_whole - start_point_whole).count();
#endif
			//exit(SUCCESS);
		} // While(1)
	} // Actual Loop
}

#ifdef MULTIMAP
void createRegularMap(std::multimap<ULL,int>& RegularMap) {

	assert(RegularMap.empty());

	Particle* ptcl;
	for (int i=0; i<=LastParticleIndex; i++)
	{
		ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;

		RegularMap.insert({ptcl->CurrentBlockReg + ptcl->TimeBlockReg, ptcl->ParticleIndex});
	}
	NextRegTimeBlock = RegularMap.begin()->first;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
#ifdef CUDA
void getRegularList(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList) {

	// assert(RegularMap.size() == NumberOfParticle);
	if (RegularMap.size() != NumberOfParticle) { // PISN case

		int num_erased = RegularMap.size() - NumberOfParticle;
		int num = 0;

		fprintf(stdout, "PISN search (number: %d) in RegularMap...\n", num_erased);

		auto it = RegularMap.begin();
		while (it != RegularMap.end()) {
			if (particles[it->second].Mass > 0) {
				it++;
			}
			else {
				fprintf(stdout, "PISN (PID: %d) is erased in RegularMap\n", particles[it->second].PID);
				it = RegularMap.erase(it);
				num++;
				if (num == num_erased)
					break;
			}
		}
	}
	assert(RegularMap.begin()->first == NextRegTimeBlock);
	assert(RegularList.empty());

	auto it = RegularMap.begin();
	while (it->first == NextRegTimeBlock) {
		RegularList.insert(it->second);
		it = RegularMap.erase(it);
	}
	assert(RegularMap.size() + RegularList.size() == NumberOfParticle);
}

void updateRegularMap(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList) {

	assert(RegularMap.size() + RegularList.size() == NumberOfParticle);

	Particle* ptcl;
	while (!RegularList.empty()) {
		ptcl = &particles[*RegularList.begin()];
		RegularMap.insert({ptcl->CurrentBlockReg + ptcl->TimeBlockReg, ptcl->ParticleIndex});
		RegularList.erase(RegularList.begin());
	}
	assert(RegularList.empty());
	assert(RegularMap.size() == NumberOfParticle);

	NextRegTimeBlock = RegularMap.begin()->first;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
#else // no cuda
void getRegularList(std::multimap<ULL,int>& RegularMap, std::vector<int>& RegularList) {

	// assert(RegularMap.size() == NumberOfParticle);
	if (RegularMap.size() != NumberOfParticle) { // PISN case

		int num_erased = RegularMap.size() - NumberOfParticle;
		int num = 0;

		fprintf(stdout, "PISN search (number: %d) in RegularMap...\n", num_erased);

		auto it = RegularMap.begin();
		while (it != RegularMap.end()) {
			if (particles[it->second].Mass > 0) {
				it++;
			}
			else {
				fprintf(stdout, "PISN (PID: %d) is erased in RegularMap\n", particles[it->second].PID);
				it = RegularMap.erase(it);
				num++;
				if (num == num_erased)
					break;
			}
		}
	}
	assert(RegularMap.begin()->first == NextRegTimeBlock);
	assert(RegularList.empty());

	auto it = RegularMap.begin();
	while (it->first == NextRegTimeBlock) {
		RegularList.push_back(it->second);
		it = RegularMap.erase(it);
	}
	assert(RegularMap.size() + RegularList.size() == NumberOfParticle);
}

void updateRegularMap(std::multimap<ULL,int>& RegularMap, std::vector<int>& RegularList) {

	assert(RegularMap.size() + RegularList.size() == NumberOfParticle);

	Particle* ptcl;
	while (!RegularList.empty()) {
		ptcl = &particles[RegularList.back()];
		RegularMap.insert({ptcl->CurrentBlockReg + ptcl->TimeBlockReg, ptcl->ParticleIndex});
		RegularList.pop_back();
	}
	assert(RegularList.empty());
	assert(RegularMap.size() == NumberOfParticle);

	NextRegTimeBlock = RegularMap.begin()->first;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
#endif // cuda
#else // no multimap
// #ifdef CUDA
void updateNextRegTime(std::unordered_set<int>& RegularList) {

	ULL time_tmp=0, time=block_max;
	Particle *ptcl;

	RegularList.clear();

	for (int i=0; i<=LastParticleIndex; i++)
	{
		//std::cout << i << std::endl;
		ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;
		// Next regular time step
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;

		// Find the minum regular time step
		if (time_tmp <= time) {
			//fprintf(stderr, "PID=%d, time_tme=%llu\n", ptcl->PID, time_tmp);
			if (time_tmp < time) {
				RegularList.clear();
				time = time_tmp;
			}
			//RegularList.push_back(ptcl->ParticleIndex);
			RegularList.insert(ptcl->ParticleIndex);
		}
	}
	NextRegTimeBlock = time;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
// #else // no cuda
void updateNextRegTime(std::vector<int>& RegularList) {

	ULL time_tmp=0, time=block_max;
	Particle *ptcl;

	RegularList.clear();

	for (int i=0; i<=LastParticleIndex; i++)
	{
		//std::cout << i << std::endl;
		ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;
		// Next regular time step
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;

		// Find the minum regular time step
		if (time_tmp <= time) {
			//fprintf(stderr, "PID=%d, time_tme=%llu\n", ptcl->PID, time_tmp);
			if (time_tmp < time) {
				RegularList.clear();
				time = time_tmp;
			}
			RegularList.push_back(ptcl->ParticleIndex);
			// RegularList.insert(ptcl->ParticleIndex);
		}
	}
	NextRegTimeBlock = time;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
// #endif // cuda
#endif // multimap

bool createSkipList(SkipList *skiplist) {

	bool debug = false;
	//fprintf(stdout, "create level starts!\n");
	//fflush(stdout);

	if (debug) {
		fprintf(stdout, "create level starts!\n");
		fflush(stdout);
	}

	Particle* ptcl;

	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl =  &particles[i];

		// if ((ptcl->NumberOfNeighbor != 0) && (ptcl->NextBlockIrr <= NextRegTimeBlock)) { // IAR original
		if (ptcl->isActive && ptcl->NextBlockIrr <= NextRegTimeBlock) {	// IAR modified
			//fprintf(stdout, "PID=%d, NBI=%llu\n", ptcl->PID, ptcl->NextBlockIrr);
			if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
				skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);
		}
	}

	/*
	if (debug) {
		skiplist->display();
	}
	*/

	//fprintf(stdout, "create level ends!\n");
	//fflush(stdout);

	if (skiplist->getFirstNode() == nullptr)
		return FAIL;
	else
		return SUCCESS;
}



bool updateSkipList(SkipList *skiplist, int ptcl_id) {
	bool debug = false;
	if (debug) {
		fprintf(stdout, "update level starts!\n");
		fflush(stdout);
	}

	/* Update New Time Steps */
	//Node* ThisLevelNode = skiplist->getFirstNode();

	/*
	if (this->debug) {
	fprintf(stdout, "PID=%d, NBI=%llu, size=%lu\n", ptcl->PID, ptcl->NextBlockIrr, ThisLevelNode->particle_list.size());
	fprintf(stdout, "NextBlockIrr=%llu\n",ptcl->NextBlockIrr);
	fflush(stdout);
	}
	*/

	Particle * ptcl = &particles[ptcl_id];

	//std::cout << "NextBlockIrr of "<< ptcl_id<<" = " << ptcl->NextBlockIrr << std::endl;
	if (ptcl->NextBlockIrr > NextRegTimeBlock)
		return true;

	if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
		skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);

	if (debug) {
	}


	if (debug) {
		//skiplist->display();
		//fprintf(stderr, "This is it.\n\n\n");
	}

	return true;
}