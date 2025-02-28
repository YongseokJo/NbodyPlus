#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cassert>
#include <mpi.h>
#include <cuda_runtime.h>

#include "global.h"
#include "SkipList.h"
#include "Worker.h"
#include "performance.h"
//#include "QueueScheduler.h"
#include "QueueManager.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

#define noDEBUG

void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int* data, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void ParticleSynchronization();
void updateNextRegTime(std::unordered_set<int>& RegularList);
bool createSkipList(SkipList *skiplist);
bool updateSkipList(SkipList *skiplist, int ptcl_id);
int writeParticle(double current_time, int outputNum);
void calculateRegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueManager &qmanager);

void formPrimordialBinaries(int beforeLastParticleIndex);
void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls, std::unordered_map<int, int>& existing, std::unordered_map<int, int>& terminated);
void FBTermination(Particle* ptclCM);
void Merge(Particle* p1, Particle* p2);

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
	


	std::unordered_set<int> RegularList;
	//MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	//MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object

	//int sender_rank, sender_tag;
	int ptcl_id;

	workers = new Worker[NumberOfWorker+1];

	for (int i=0; i<=NumberOfWorker; i++) {
		workers[i].initialize(i);
	}

	QueueManager qmanager;

	/* Particle loading Check */
	/*
	{
		//, NextRegTime= %.3e Myr(%llu),
		for (int i=0; i<=LastParticleIndex; i++) {
			ptcl = &particles[i];
			fprintf(stdout, "PID=%d, pos=(%lf, %lf, %lf), vel=(%lf, %lf, %lf)\n",
					ptcl->PID,
					ptcl->Position[0],
					ptcl->Position[1],
					ptcl->Position[2],
					ptcl->Velocity[0],
					ptcl->Velocity[1],
					ptcl->Velocity[2]
					);
		}
		fflush(stdout);
	}*/

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
#ifdef PerformanceTrace
		performance.start(InitAcc1_Root);
#endif
		qmanager.initialize(PIDs, InitAcc1);
		qmanager.signalWorkers();
		qmanager.reportProgress();
#ifdef PerformanceTrace
		performance.end(InitAcc1_Root);
#endif

		std::cout << "InitAcc1 done" << std::endl;
		std::cerr << "Init1 Time for Root (ns) = " << performance.get(InitAcc1_Root) << std::endl;

		qmanager.initialize(PIDs, InitAcc2);
		qmanager.signalWorkers();
		qmanager.reportProgress();

		std::cout << "InitAcc2 done" << std::endl;

#ifdef _FEWBODY
		// Primordial binary search

		qmanager.initialize(PIDs, SearchPrimordialGroup);
		qmanager.signalWorkers();
		qmanager.reportProgress();

		std::cout << "Primordial binary search done" << std::endl;

		// example code by EW 2025.1.7
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

		qmanager.initialize(PIDs, InitTime);
		qmanager.signalWorkers();
		qmanager.reportProgress();
		std::cout << "InitTime done" << std::endl;
	}

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

	/* timestep variable synchronization */
	{
		std::cout << "Time Step synchronization." << std::endl;
		tasks[0] = TimeSync;
		global_variable->QueueSize = NumberOfWorker;
		qmanager.completed_queues = 0;
		qmanager.signalWorkers();
		broadcastFromRoot(time_block);
		broadcastFromRoot(block_max);
		broadcastFromRoot(time_step);
		//MPI_Barrier(MPI_COMM_WORLD);
		//qmanager.reportProgress();
		MPI_Win_fence(0, win);
        //MPI_Barrier(MPI_COMM_WORLD);
		//MPI_Win_fence(0, win);
		fprintf(stderr, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);
		fflush(stderr);
	}

	/* Particle Initialization Check */
	/*
	{
		//, NextRegTime= %.3e Myr(%llu),
		for (int i = 0; i <= LastParticleIndex; i++)
		{
			ptcl = &particles[i];
			fprintf(stdout, "%d(%d)=", ptcl->PID, ptcl->NumberOfNeighbor);
			for (int j = 0; j < ptcl->NumberOfNeighbor; j++)
			{
				fprintf(stdout, "%d, ", ptcl->Neighbors[j]);
			}
			fprintf(stdout, "\n");
		}
	}

	for (int i = 0; i <= LastParticleIndex; i++)
	{
		ptcl = &particles[i];
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
	}
	fflush(stdout);
	*/


	/* Actual Loop */
	{
		int max_level = 5;
		double prob = 0.5;
		SkipList *skiplist;
		Node* ThisLevelNode;
		double current_time_irr=0;
		double next_time=0;
		Worker* worker;
		Queue queue;
		std::unordered_set<int> CMPtclsForComputation;


		//ParticleSynchronization();
		while (1) {

			global_time = NextRegTimeBlock*time_step;

			// create output at appropriate time intervals
			if (global_time >= outputTime) {
				writeParticle(global_time, outNum++);
				outputTime += outputTimeStep;
			}

			// end if the global time exceeds the end time
			if (global_time >= 1) {
				tasks[0] = Ends;
				InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
				//MPI_Waitall(NumberOfCommunication, requests, statuses);
				//NumberOfCommunication = 0;
				std::cout << EnzoTimeStep << std::endl;
				std::cout << "Simulation Done!" << std::endl;
				return;
			}

			updateNextRegTime(RegularList);
			/*
			std::cout << "NextRegTimeBlock=" << NextRegTimeBlock << std::endl;
			std::cout << "PID= ";
			for (int i=0; i<RegularList.size(); i++)
				std::cout << RegularList[i]<< ", ";
			std::cout << std::endl;
			std::cout << "size of regularlist= " << RegularList.size() << std::endl;
			*/

			skiplist = new SkipList(max_level, prob);
			if (createSkipList(skiplist) == FAIL)
				fprintf(stderr, "There are no irregular particles!\nBut is it really happening? check skiplist->display()\n");

			bool bin_termination = false;
			bool new_binaries = false;


			// Irregular
			while ( skiplist->getFirstNode() != nullptr) {
				// update_idx=0; // commented out by EW 2025.1.6
				ThisLevelNode = skiplist->getFirstNode();
				next_time     = particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr\
									 	    + particles[ThisLevelNode->ParticleList[0]].TimeStepIrr;
				global_variable->next_time = next_time;
#ifdef DEBUG
				// print out particlelist
				fprintf(stdout, "(IRR_FORCE) next_time: %e Myr\n", next_time*EnzoTimeStep*1e4);
				/*
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
				*/
#endif

				// Irregular Force
#ifdef FEWBODY
#ifdef NSIGHT
				nvtxRangePushA("IrregularForce");
#endif

#ifdef PerformanceTrace
				performance.start(IrrForce_Root);
#endif
#ifdef DEBUG
				std::cout << "Irr force starts" << std::endl;
#endif

#ifdef IRR_CM_SEPARATE

				CMPtclsForComputation.clear();
				global_variable->QueueSize = ThisLevelNode->ParticleList.size();

				//std::cout << "Queues = ";
				for (int i = 0; i < ThisLevelNode->ParticleList.size(); i++)
				{
					if (!particles[ThisLevelNode->ParticleList[i]].isActive)
					{
						//std::cout << "PID "<< ThisLevelNode->ParticleList[i] <<" is inactive." << std::endl;
						assert(particles[ThisLevelNode->ParticleList[i]].isActive);
					}
					if (particles[ThisLevelNode->ParticleList[i]].isCMptcl)
					{
						CMPtclsForComputation.insert(ThisLevelNode->ParticleList[i]);
					}
					queues[i] = ThisLevelNode->ParticleList[i];
					//std::cout << queues[i] << ", ";
				}
				//std::cout << std::endl;
				tasks[0] = IrrForce;
				qmanager.completed_queues = 0;
				qmanager.NumberOfCommunication = 0;
				MPI_Win_sync(win4);
				MPI_Win_sync(win5);
				MPI_Win_flush_all(win4);
				MPI_Win_flush_all(win5);
				//qmanager.initialize(ThisLevelNode->ParticleList, IrrForce);
				qmanager.signalWorkers();
				qmanager.reportProgress();

				if (CMPtclsForComputation.size() > 0)
				{
					{
						int i = 0;
						for (int cm_pid : CMPtclsForComputation)
						{
							queues[i] = cm_pid;
							std::cout << "cm_pid = " << cm_pid << std::endl;
						}
					}
					global_variable->QueueSize = CMPtclsForComputation.size();
					tasks[0] = ARIntegration;
					qmanager.completed_queues = 0;
					qmanager.NumberOfCommunication = 0;
					qmanager.signalWorkers();
					qmanager.reportProgress();
				}

#else // IRR_CM_SEPARATE
 
#endif // IRR_CM_SEPARATE


#ifdef DEBUG
				std::cout << "Irregular Force done" << std::endl;
#endif
#ifdef PerformanceTrace
				performance.end(IrrForce_Root);
#endif

		//std::cerr << "IrrForce Time for Root (ns) = " << performance.get(IrrForce_Root) << std::endl;
#ifdef NSIGHT
				nvtxRangePop();
#endif 
#else // FEWBODY
#ifdef NSIGHT
				nvtxRangePushA("IrregularForce");
#endif
#ifdef PerformanceTrace
				performance.start(IrrForce_Root);
#endif
				qmanager.initialize(ThisLevelNode->ParticleList, IrrForce);
				qmanager.signalWorkers();
				qmanager.reportProgress();
#ifdef PerformanceTrace
				performance.end(IrrForce_Root);
#endif
#ifdef NSIGHT
				nvtxRangePop();
#endif 
#endif // FEWBODY

				// if remaining, further sorting
				/*
				if (update_idx < total_tasks) {
						updateSkipList(skiplist, ThisLevelNode->ParticleList[update_idx]);
						update_idx++;
					}				updateSkipList(skiplist, ptcl_id_return);
					*/

#ifdef NSIGHT
				nvtxRangePushA("IrregularUpdate");
#endif
#ifdef PerformanceTrace
				performance.start(IrrUpdate1);
#endif
				// Irregular Update
				qmanager.initialize(ThisLevelNode->ParticleList, IrrUpdate);
				/*
				std::cout << "Queues = ";
				for (int i=0; i<global_variable->QueueSize;i++)
					std::cout << queues[i] << " ";
				std::cout << std::endl; 
				*/
				qmanager.signalWorkers();
				qmanager.reportProgress();
#ifdef PerformanceTrace
				performance.end(IrrUpdate1);
#endif
#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef NSIGHT
				nvtxRangePushA("Irregular update on root");
#endif

#ifdef PerformanceTrace
				performance.start(IrrUpdate2);
#endif
				for (int ptcl_id : ThisLevelNode->ParticleList)
				{
					ptcl = &particles[ptcl_id];
					if (ptcl->NumberOfNeighbor != 0) // IAR modified
						ptcl->updateParticle();
					ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockIrr * time_step;
				}

#ifdef PerformanceTrace
				performance.end(IrrUpdate2);
#endif
#ifdef NSIGHT
				nvtxRangePop();
#endif
#ifdef DEBUG
				std::cerr <<  "ParticleList Size=" << ThisLevelNode->ParticleList.size() << std::endl;
				for (int i : ThisLevelNode->ParticleList)
				{
					std::cerr << i << " ";
					ptcl = &particles[i];
					if (ptcl->CurrentTimeIrr != next_time)
					{
						fprintf(stdout, "Error! PID: %d, CurrentTimeIrr: %e Myr, next_time: %e Myr\n", ptcl->PID, ptcl->CurrentTimeIrr * EnzoTimeStep * 1e4, next_time * EnzoTimeStep * 1e4);
						assert(ptcl->CurrentTimeIrr == next_time);
					}
				}
				std::cerr << std::endl;
				std::cout << "Irregular update done" << std::endl;
#endif


				//std::cerr << "IrrUpdate1 (ns) = " << performance.get(IrrUpdate1) << std::endl;
				//std::cerr << "IrrUpdate2 (ns) = " << performance.get(IrrUpdate2) << std::endl;
				//exit(1);

#define noTEST_1
#ifdef TEST_1
#ifdef NSIGHT
				nvtxRangePushA("FewBodySearch");
#endif

#ifdef DEBUG
				std::cout << "FB search starts" << std::endl;
#endif
#ifdef PerformanceTrace
				//performance.start(IrrUpdate1);
#endif
				qmanager.initialize(ThisLevelNode->ParticleList, SearchGroup);
				qmanager.signalWorkers();
				qmanager.reportProgress();
#ifdef PerformanceTrace
				//performance.end(IrrUpdate1);
#endif
#ifdef DEBUG
				std::cout << "FB search ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif
#endif





#ifdef _FEWBODY	

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
							if (ptcl->NewNumberOfNeighbor == 2) { // binary merger

								Particle* donor = &particles[ptcl->NewNeighbors[0]];
								Particle* accretor = &particles[ptcl->NewNeighbors[1]];

								Merge(donor, accretor);
								
							}
							else { // from NewFBInitialization3

								assert(ptcl->NewNumberOfNeighbor > 2); // for debugging by EW 2025.1.20

								Particle* donor;
								Particle* accretor;

								for (int j=0; j<ptcl->NewNumberOfNeighbor; j++) {
									if (particles[ptcl->NewNeighbors[j]].getBinaryInterruptState() == BinaryInterruptState::collision) {
										donor = &particles[ptcl->NewNeighbors[j]];
										accretor = &particles[donor->getBinaryPairID()];

										assert(accretor->getBinaryInterruptState() == BinaryInterruptState::collision);
										assert(accretor->getBinaryPairID() == donor->ParticleIndex);
										break;
									}
								}

								Merge(donor, accretor);

								// (Query to myself) this can be parallelized 
								int rank = CMPtclWorker[ptcl->ParticleIndex];
								qmanager.completed_queues = 0;
								tasks[0] = MergeManyBody;
								global_variable->QueueSize = 1;
								queues[rank-1] = ptcl->ParticleIndex;
								qmanager.signalWorker(rank);
								qmanager.reportProgress();

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

						for (int j=0; j < ptcl->NewNumberOfNeighbor; j++) {
							if (particles[ptcl->NewNeighbors[j]].Mass == 0.0)
								continue;
							ThisLevelNode->ParticleList.push_back(ptcl->NewNeighbors[j]);
							particles[ptcl->NewNeighbors[j]].isActive = true;
						}

						FBTermination(ptcl);
						ptcl->isActive = false;
					}
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
#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef NSIGHT
				nvtxRangePushA("FewBodySearch");
#endif

#ifdef DEBUG
				std::cout << "FB search starts" << std::endl;
#endif
#ifdef PerformanceTrace
				//performance.start(IrrUpdate1);
#endif
				qmanager.initialize(ThisLevelNode->ParticleList, SearchGroup);
				qmanager.signalWorkers();
				qmanager.reportProgress();
#ifdef PerformanceTrace
				//performance.end(IrrUpdate1);
#endif
#ifdef DEBUG
				std::cout << "FB search ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
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
					for (int i=OriginalParticleListSize; i<ThisLevelNode->ParticleList.size(); i++) {
						ptcl = &particles[ThisLevelNode->ParticleList[i]];
						fprintf(stdout, "New CM Particle PID: %d\n", ptcl->PID);
						fflush(stdout);
					}

					// new code by EW 2025.1.26
					Particle* ptclCM;
					Particle* mem_ptclCM;
					for (int i=0; i<newCMptcls.size(); i++) {
						ptclCM = &particles[newCMptcls[i]]; // 2025.01.10 edited to newCMptcls[i] by YS

						for (int j=0; j<ptclCM->NewNumberOfNeighbor; j++) {
							mem_ptclCM = &particles[ptclCM->NewNeighbors[j]];
							if (mem_ptclCM->isCMptcl) {
								fprintf(stdout, "manybody group detected; PID %d should be deleted first\n", mem_ptclCM->PID);
								rank_delete = CMPtclWorker[mem_ptclCM->ParticleIndex];
								fprintf(stdout, "Rank of CM ptcl %d: %d\n", mem_ptclCM->PID, rank_delete);

								queue.task = DeleteGroup;
								queue.pid = mem_ptclCM->ParticleIndex;
								workers[rank_delete].addQueue(queue);
								workers[rank_delete].runQueue();
								workers[rank_delete].callback();

								// (Query to myself) this can be parallelized
								qmanager.completed_queues = 0;
								tasks[0] = DeleteGroup;
								global_variable->QueueSize = 1;
								queues[rank_delete-1] = mem_ptclCM->ParticleIndex;
								qmanager.signalWorker(rank_delete);
								qmanager.reportProgress();

								PrevCMPtclWorker.insert({mem_ptclCM->ParticleIndex, CMPtclWorker[mem_ptclCM->ParticleIndex]});
								CMPtclWorker.erase(mem_ptclCM->ParticleIndex);
							}
						}

						rank_new = CMPtclWorker[ptclCM->ParticleIndex];
						// fprintf(stdout, "New CM ptcl is %d\n", newCMptcls[0]);
						fprintf(stdout, "Rank of CM ptcl %d: %d\n", ptclCM->PID, rank_new);
						// (Query to myself) this can be parallelized
						qmanager.completed_queues = 0;
						tasks[0] = MakeGroup;
						global_variable->QueueSize = 1;
						queues[rank_new-1] = ptclCM->ParticleIndex;
						qmanager.signalWorker(rank_new);
						qmanager.reportProgress();
					}
					std::cout << "All new fewbody objects are initialized." << std::endl;

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

				// std::cout << "erase success" << std::endl;
#ifdef NSIGHT
				nvtxRangePop();
#endif
#endif

#ifdef DEBUG
				std::cout << "updateSkipList starts" << std::endl;
#endif
				for (int i=0; i<ThisLevelNode->ParticleList.size(); i++)
					updateSkipList(skiplist, ThisLevelNode->ParticleList[i]);
#ifdef DEBUG
				std::cout << "updateSkipList ended" << std::endl;
#endif

				//std::cout << "update success" << std::endl;
				//skiplist->display();
				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
						ptcl = &particles_original[ThisLevelNode->ParticleList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e10/1e6,
								ptcl->NextBlockIrr,
								ptcl->NumberOfNeighbor
								);

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
								ptcl->a_irr[2][3]
									);
					}
					fflush(stdout);
				}
				*/
				current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr*time_step;
#ifdef DEBUG
				std::cout << "skiplist->deleteFirstNode() starts" << std::endl;
#endif
				skiplist->deleteFirstNode();
#ifdef DEBUG
				std::cout << "skiplist->deleteFirstNode() ended" << std::endl;
#endif

				//std::cout << "deleteFirstNode success" << std::endl;
				//ThisLevelNode = skiplist->getFirstNode();

				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
						ptcl = &particles[ThisLevelNode->ParticleList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NumberOfNeighbor
								);
					}
				}
				*/
				//exit(SUCCESS);
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
					InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
					MPI_Waitall(NumberOfCommunication, requests, statuses);
					NumberOfCommunication = 0;
					std::cout << EnzoTimeStep << std::endl;
					std::cout << "Simulation Done!" << std::endl;
					return;
				}
#endif
			} // Irr
#ifdef DEBUG
			delete skiplist;
			std::cout << "delete skiplist" << std::endl;
#endif
			skiplist = nullptr;
			//exit(SUCCESS);

#ifdef FEWBODY
			if (bin_termination || new_binaries) {
#ifdef DEBUG
			std::cout << "(FB) updateNextRegTime starts" << std::endl;
#endif
				updateNextRegTime(RegularList);
#ifdef DEBUG
			std::cout << "(FB) updateNextRegTime done" << std::endl;
			std::cout << "(FB) RegularList size: " << RegularList.size() << std::endl;
#endif
			}
#endif

#ifdef CUDA
			{
				//total_tasks = RegularList.size();
				next_time = NextRegTimeBlock*time_step;

				//fprintf(stdout, "Regular starts\n");
#ifdef NSIGHT
				nvtxRangePushA("calculateRegAccelerationOnGPU");
#endif

#ifdef DEBUG
				std::cout << "calculateRegAccelerationOnGPU starts" << std::endl;
				std::cout << "RegularList size: " << RegularList.size() << std::endl;
#endif

				calculateRegAccelerationOnGPU(RegularList, qmanager);

#ifdef DEBUG
				std::cout << "calculateRegAccelerationOnGPU ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef NSIGHT
				nvtxRangePushA("UpdateRegular");
#endif

#ifdef DEBUG
				std::cout << "update regular starts" << std::endl;
#endif

				// Update Regular
				qmanager.initialize(RegularList, RegCudaUpdate);
				qmanager.signalWorkers();
				qmanager.reportProgress();
#ifdef DEBUG
				std::cout << "update regular ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif
			}
				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<RegularList.size(); i++) {
						ptcl = &particles[RegularList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e10/1e6,
								ptcl->NextBlockIrr,
								ptcl->NumberOfNeighbor
								);

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
								ptcl->a_irr[2][3]
									);
					}
					//fflush(stdout); 
				}
				*/
#else
			//std::cout << "Regular Routine Starts." << std::endl;
			// Regular
			{
				// Regular Gravity
				task = RegForce;
				completed_tasks = 0;
				total_tasks = RegularList.size();
				next_time = NextRegTimeBlock*time_step;

				//std::cout << "TotalTask=" << total_tasks << std::endl;

				/*
				std::cout << "RegularList, PID= ";
				for (int i=0; i<total_tasks; i++) {
					std::cout << RegularList[i]<< ", ";
				}*/
				//std::cout << std::endl;

				InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
				InitialAssignmentOfTasks(RegularList, next_time, total_tasks, PTCL_TAG);
				MPI_Waitall(NumberOfCommunication, requests, statuses);
				NumberOfCommunication = 0;


				// further assignments
				remaining_tasks = total_tasks-NumberOfWorker;
				while (completed_tasks < total_tasks) {
					// Check which worker is done
					MPI_Irecv(&ptcl_id_return, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
					// Poll until send completes
					/*
						 flag=0;
						 while (!flag) {
						 MPI_Test(&request, &flag, &status);
					// Perform other work while waiting
					}
					*/
					MPI_Wait(&request, &status);
					completed_rank = status.MPI_SOURCE;
					//printf("Rank %d: Send operation completed (%d).\n",completed_rank, ptcl_id_return);

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


				//ParticleSynchronization();


				// Regular Update
				//std::cout<< "Reg Acc Done." <<std::endl;
				task = RegUpdate;
				completed_tasks = 0;

				InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
				InitialAssignmentOfTasks(RegularList, total_tasks, PTCL_TAG);
				MPI_Waitall(NumberOfCommunication, requests, statuses);
				NumberOfCommunication = 0;


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

				//ParticleSynchronization();
				/*
				{
					//, NextRegTime= %.3e Myr(%llu),
					for (int i=0; i<total_tasks; i++) {
						ptcl = &particles[RegularList[i]];
						fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
								"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
								"NumNeighbor= %d\n",
								ptcl->PID,
								ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockIrr,
								ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
								ptcl->CurrentBlockReg,
								NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
								NextRegTimeBlock,
								ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
								ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
								ptcl->TimeBlockIrr,
								ptcl->TimeLevelIrr,
								ptcl->TimeBlockReg,
								ptcl->TimeLevelReg,
								ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e10/1e6,
								ptcl->NextBlockIrr,
								ptcl->NumberOfNeighbor
								);

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
								ptcl->NewNumberOfNeighbor,
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
								ptcl->a_irr[2][3]
									);
					}
					//fflush(stdout); 
				}
				*/
				//current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr;
			} // Regular Done.
#endif

#ifdef SEVN
			StellarEvolution(); // How about evolving particles inside RegularList only? by EW 2025.1.19
								// Currently, evolving all the particles upto global_time
#endif
			//exit(SUCCESS);
		} // While(1)
	} // Actual Loop
}




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
}







bool createSkipList(SkipList *skiplist) {

	bool debug = false;
	//fprintf(stdout, "create level starts!\n");
	//fflush(stdout);

	if (debug) {
		fprintf(stdout, "create level starts!\n");
		fflush(stdout);
	}

	/*
#ifdef time_trace
	_time.irr_chain.markStart();
#endif
*/

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
#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif
*/

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

	/*
#ifdef time_trace
	_time.irr_sort.markStart();
#endif
*/

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

	/*
#ifdef time_trace
	_time.irr_sort.markEnd();
	_time.irr_sort.getDuration();
#endif
*/

	if (debug) {
		//skiplist->display();
		//fprintf(stderr, "This is it.\n\n\n");
	}

	return true;
}




