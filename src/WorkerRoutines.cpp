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
#ifdef CUDA
void sendAllParticlesToGPU_Worker(double new_time);
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

				std::memcpy(Neighbors + ptcl->NeighborsOffset, NewNeighbors + ptcl->NeighborsOffset, sizeof(int) * ptcl->NewNumberOfNeighbor);
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

				ptcl->CurrentBlockReg += ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg = ptcl->CurrentBlockReg * time_step;

				ptcl->updateParticle();
				std::memcpy(Neighbors + ptcl->NeighborsOffset, NewNeighbors + ptcl->NeighborsOffset, sizeof(int) * ptcl->NewNumberOfNeighbor);
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->calculateTimeStepReg();
				ptcl->calculateTimeStepIrr();
				// /*
				if (ptcl->CurrentBlockIrr != ptcl->CurrentBlockReg || ptcl->CurrentTimeIrr != ptcl->CurrentTimeReg) {
					fprintf(stderr, "PID: %d\n", ptcl->PID);
					fprintf(stderr, "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n", ptcl->CurrentBlockIrr, ptcl->CurrentBlockReg);
					fprintf(stderr, "CurrentBlockIrr * time_step: %.17g, CurrentBlockReg * time_step: %.17g\n", ptcl->CurrentBlockIrr*time_step, ptcl->CurrentBlockReg*time_step);
					fprintf(stderr, "CurrentTimeIrr: %.17g, CurrentTimeReg: %.17g\n", ptcl->CurrentTimeIrr, ptcl->CurrentTimeReg);
					fprintf(stderr, "NextRegTimeBlock: %llu\n", global_variable->NextRegTimeBlock);
					fflush(stderr);
					assert(ptcl->CurrentTimeIrr == ptcl->CurrentTimeReg);
					assert(ptcl->CurrentBlockIrr == ptcl->CurrentBlockReg);
				}
				// */
				ptcl->updateRadius();
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of ptcl particle
				break;

			case InitAcc1: // Initialize Acceleration(01)

				ptcl = &particles[ptcl_id];
				CalculateAcceleration01(ptcl);
				break;

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
				fprintf(workerout, "MyRank = %d time_block = %d, EnzoTimeStep = %e\n\n", MyRank, time_block, EnzoTimeStep);
				break;

#ifdef FEWBODY
			case SearchPrimordialGroup: // Primordial binary search

				ptcl = &particles[ptcl_id];

				ptcl->NewNumberOfMember = 0;
				ptcl->checkNewGroup2();
				break;

			case SearchGroup: // Few-body group search

				ptcl = &particles[ptcl_id];

				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::threebody) {
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
				}
				else if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->NewNumberOfMember = 0;
					ptcl->checkNewGroup2();
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
				}
				else {
					ptcl->NewNumberOfMember = 0;
					if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < TSEARCH)
						ptcl->checkNewGroup();
				}
				/*
				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
					std::cout << "ptcl PID: " << ptcl->PID << ", ptcl NewNumberOfMember: " << ptcl->NewNumberOfMember << std::endl;
				}
				else {
					ptcl->NewNumberOfMember = 0;
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
					ptcl->GroupInfo->isTerminate = ptcl->GroupInfo->CheckBreak2();

				if (ptcl->GroupInfo->isTerminate) {
					if (ptcl->getBinaryInterruptState() == BinaryInterruptState::none)
						ptcl->setBinaryInterruptState(BinaryInterruptState::terminated);

					E_binary -= ptcl->GroupInfo->sym_int.getEtot();
					E_binary_SD -= ptcl->GroupInfo->sym_int.getEtotSlowDown();

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

			case GetTotalEnergy:

				MPI_Reduce(&E_binary,		nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&E_binary_SD,	nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&E_merger,		nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&E_PN,			nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				continue;
#ifdef CUDA
			case PrepareGPUCalc:

				sendAllParticlesToGPU_Worker(next_time);
				continue;
#endif
			case Synchronize: // Synchronize
				MPI_Win_sync(win);  // Synchronize memory
				MPI_Barrier(shared_comm);
				break;

			case Ends: // Simualtion ends
				fprintf(workerout, "Processor %d returns.\n", MyRank);
				return;

			case Error:
				perror("Error task assignments");
				exit(EXIT_FAILURE);
				break;

			default:
				break;
		}

		// return that it's over
		//task = -1;
		if (task == IrrForce || task == RegForce || task == IrrUpdate || task == RegUpdate)
			MPI_Isend(&ptcl_id, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD,&request);
		else
			MPI_Isend(&task, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD,&request);

		MPI_Wait(&request, &status);
		//std::cerr << "Processor " << MyRank << " done." << std::endl;
	}
}