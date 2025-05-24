#include <iostream>
#include <vector>
#include <errno.h>
#include "global.h"
#include "Queue.h"

void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void makePrimordialGroup(Particle* ptclCM);
void NewFBInitialization(Particle* ptclCM);
void deleteGroup(Particle* ptclCM);
void NewFBInitialization3(Group* group);
#ifdef MultiNode
void updateInitAcc01(int update_count, int* update_count_list, int* displs);
void updateInitAcc23(int update_count, int* update_count_list, int* displs);
void sendNewNeighbors(int update_count, std::vector<int>& neighbors);
void updateTimeVariables(int update_count, int* update_count_list, int* displs);
void updateTimeCorrection(int update_count);
#endif

void WorkerRoutines() {

	//std::cout << "Processor " << MyRank << " is ready." << std::endl;

	TaskName task = Error;
	MPI_Status status;
	MPI_Request request;
	int ptcl_id;
	double next_time;
	Particle *ptcl;
	Queue queue;
#ifdef MultiNode
	int* update_count_list = new int[NumberOfNode];
	int* displs = new int[NumberOfNode];
	std::vector<int> neighbors;
#endif

	while (true) {

		MPI_Recv(&queue, 1, QueueType, ROOT, QUEUE_TAG, MPI_COMM_WORLD, &status);
		task = queue.task;
		ptcl_id = queue.pid;
		next_time = queue.next_time;

		switch (task) {
			case IrrForce: // Irregular Acceleration

				ptcl = &particles[ptcl_id];
				ptcl->computeAccelerationIrr();

				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->calculateTimeStepIrr();
				ptcl->NextBlockIrr = ptcl->NewCurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->isUpdateToDate = true;
				break;

			case RegForce: // Regular Acceleration

				ptcl = &particles[ptcl_id];
				ptcl->computeAccelerationReg();
				break;

			case IrrUpdate: // Irregular Update Particle

				ptcl = &particles[ptcl_id];
				if (ptcl->NumberOfNeighbor != 0) // IAR modified
					ptcl->updateParticle();
				ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
				ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
				break;

			case RegUpdate: // Regular Update Particle

				ptcl = &particles[ptcl_id];
				ptcl->updateParticle();

				for (int i=0; i<ptcl->NewNumberOfNeighbor; i++)
					ptcl->Neighbors[i] = ptcl->NewNeighbors[i];
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->CurrentBlockReg += ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg   = ptcl->CurrentBlockReg*time_step;
				ptcl->calculateTimeStepReg();
				// ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockReg; // commented out by EW 2025.3.3 to match with RegCudaUpdate task
				ptcl->calculateTimeStepIrr();
				ptcl->updateRadius();
				if (ptcl->NumberOfNeighbor == 0) {
					ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*time_step;
				}
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				break;

			case RegCuda: // Update Regular Particle CUDA

				ptcl = &particles[ptcl_id];
				ptcl->updateRegularParticleCuda();
				break;

			case RegCudaUpdate: // Update Regular Particle CUDA II

				ptcl = &particles[ptcl_id];

				for (int j = 0; j < ptcl->NewNumberOfNeighbor; j++)
					ptcl->Neighbors[j] = ptcl->NewNeighbors[j];
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->updateParticle();
				ptcl->CurrentBlockReg = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg = ptcl->CurrentBlockReg * time_step;
				ptcl->calculateTimeStepReg();
				ptcl->calculateTimeStepIrr();
				// /*
				if (ptcl->NumberOfNeighbor == 0) {
					/*
					if (ptcl->CurrentBlockIrr != ptcl->CurrentBlockReg || ptcl->CurrentTimeIrr != ptcl->CurrentBlockReg*time_step) {
						fprintf(stderr, "PID: %d\n", ptcl->PID);
						fprintf(stderr, "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n", ptcl->CurrentBlockIrr, ptcl->CurrentBlockReg);
						fprintf(stderr, "CurrentBlockIrr * time_step: %e, CurrentBlockReg * time_step: %e\n", ptcl->CurrentBlockIrr*time_step, ptcl->CurrentBlockReg*time_step);
						fprintf(stderr, "CurrentTimeIrr: %e, CurrentTimeReg: %e\n", ptcl->CurrentTimeIrr, ptcl->CurrentTimeReg);
						fprintf(stderr, "NextRegTimeBlock: %llu\n", global_variable->NextRegTimeBlock);
						fflush(stderr);
						assert(ptcl->CurrentBlockIrr == ptcl->CurrentBlockReg);
						assert(ptcl->CurrentTimeIrr == ptcl->CurrentBlockReg*time_step);
					}
					*/
					ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*time_step;
				}
				// */
				ptcl->updateRadius();
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of ptcl particle
				break;

			case InitAcc1: // Initialize Acceleration(01)

				ptcl = &particles[ptcl_id];
				CalculateAcceleration01(ptcl);
				break;
#ifdef MultiNode
			case UpdateInitAcc01:

				updateInitAcc01(ptcl_id, update_count_list, displs);
				break;

			case UpdateInitAcc23:

				updateInitAcc23(ptcl_id, update_count_list, displs);
				break;

			case SendNewNeighbors:

				sendNewNeighbors(ptcl_id, neighbors);
				break;

			case UpdateLastParticleIndex:

				LastParticleIndex = ptcl_id;
				global_variable->LastParticleIndex = LastParticleIndex;
				break;

			case UpdateTimeVariables:

				updateTimeVariables(ptcl_id, update_count_list, displs);
				break;

			case UpdateTimeCorrection:

				updateTimeCorrection(ptcl_id);
				break;
#endif

			case InitAcc2: // Initialize Acceleration(23)

				ptcl = &particles[ptcl_id];
				CalculateAcceleration23(ptcl);
				break;

			case InitTime: // Initialize Time Step

				ptcl = &particles[ptcl_id];
				if (ptcl->isActive)
					ptcl->initializeTimeStep();
				break;

			case TimeSync: // Initialize Timestep variables
				broadcastFromRoot(time_block);
				broadcastFromRoot(block_max);
				broadcastFromRoot(time_step);
				//MPI_Win_sync(win);  // Synchronize memory
				//MPI_Barrier(shared_comm);
				//MPI_Win_fence(0, win);
				fprintf(stderr, "(%d) nbody+:time_block = %d, EnzoTimeStep=%e\n", MyRank, time_block, EnzoTimeStep);
				fflush(stderr);
				break;


#ifdef FEWBODY
			case SearchPrimordialGroup: // Primordial binary search

				ptcl = &particles[ptcl_id];

				ptcl->NewNumberOfNeighbor = 0;
				ptcl->checkNewGroup2();
				break;

			case SearchGroup: // Few-body group search

				ptcl = &particles[ptcl_id];

				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::threebody) {
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
				}
				else if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->NewNumberOfNeighbor = 0;
					ptcl->checkNewGroup2();
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
				}
				else {
					ptcl->NewNumberOfNeighbor = 0;
					if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < TSEARCH)
						ptcl->checkNewGroup();
				}
				/*
				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
					std::cout << "ptcl PID: " << ptcl->PID << ", ptcl NewNumberOfNeighbor: " << ptcl->NewNumberOfNeighbor << std::endl;
				}
				else {
					ptcl->NewNumberOfNeighbor = 0;
					if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < TSEARCH)
						ptcl->checkNewGroup();
				}
				*/
				break;

			case MakePrimordialGroup: // Make a primordial group

				ptcl = &particles[ptcl_id];
				makePrimordialGroup(ptcl);
				break;

			case MakeGroup: // Make a group

				ptcl = &particles[ptcl_id];

				NewFBInitialization(ptcl);
#ifdef DEBUG
				std::cout << "FewBody object of particle " << ptcl->PID
						  << " is successfully initialized on rank " << MyRank << "." <<std::endl;
#endif
				break;

			case DeleteGroup: // Delete a Group struct

				ptcl = &particles[ptcl_id];
				deleteGroup(ptcl);
				break;

			case ARIntegration: // SDAR for few body encounters
				
				ptcl = &particles[ptcl_id];

				if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
					exit(EXIT_FAILURE);
				}
				
				ptcl->GroupInfo->ARIntegration(next_time);
				if (!ptcl->GroupInfo->isMerger && !ptcl->GroupInfo->isTerminate)
					ptcl->GroupInfo->isTerminate = ptcl->GroupInfo->CheckBreak();

				if (ptcl->GroupInfo->isTerminate) {
					if (ptcl->getBinaryInterruptState() == BinaryInterruptState::none)
						ptcl->setBinaryInterruptState(BinaryInterruptState::terminated);

					delete ptcl->GroupInfo;
#ifdef DEBUG
					std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << " deleted!" <<std::endl;
#endif
				}
#ifdef DEBUG
				else
					std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << " done!" <<std::endl;
#endif
				break;
			
			case MergeManyBody: // Merger insided many-body (>2) group
				
				ptcl = &particles[ptcl_id];
				std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << std::endl;

				if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
					exit(EXIT_FAILURE);
				}

				NewFBInitialization3(ptcl->GroupInfo);

				ptcl->GroupInfo->isMerger = false;
				ptcl->setBinaryInterruptState(BinaryInterruptState::none);

				std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << " NewFBInitialization3 done!" <<std::endl;
				break;
#endif 

			case Synchronize: // Synchronize
				MPI_Win_sync(win);  // Synchronize memory
				MPI_Barrier(shared_comm);
				break;

			case Ends: // Simualtion ends
				std::cout << "Processor " << MyRank<< " returns." << std::endl;
				return;
				break;

			case Error:
				perror("Error task assignments");
				exit(EXIT_FAILURE);
				break;
			default:
				break;
		}

		// return that it's over // (Query) EW: I think there is no need to use MPI_Isend here. Let's use MPI_Send instead. 2025.5.20
#ifdef MultiNode
		if (task == SearchPrimordialGroup) {
			int return_value[2] = {ptcl_id, ptcl->NewNumberOfNeighbor};
			MPI_Isend(return_value, 2, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);
		}
		else if (task == SendNewNeighbors) {
			MPI_Isend(neighbors.data(), neighbors.size(), MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);
		}
		else {
			MPI_Isend(&ptcl_id, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);
		}
#else
		MPI_Isend(&ptcl_id, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);
#endif

		MPI_Wait(&request, &status);
#ifdef MultiNode
		neighbors.clear();
#endif
		//std::cerr << "Processor " << MyRank << " done." << std::endl;
	}
}