#include <iostream>
#include <vector>
#include <errno.h>
#include "global.h"
#include "Queue.h"
#include <cstring>
// #include <unordered_set> // for sendAllParticlesToGPU
#include "cuda/cuda_functions.h"

void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void makePrimordialGroup(Particle* ptclCM);
void NewFBInitialization(Particle* ptclCM);
void deleteGroup(Particle* ptclCM);
void NewFBInitialization3(Group* group);

void RegularWorker(int NumTargetTotal, int Jstart, int Jend, int gpu_id);
// void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList, int N_start, int N_end, int gpu_id);
void AllocateDeviceMemory(int N_i,int N_j,int gpu_id);
void SendToDeviceMPI(int N_i, int N_j, int gpu_id);

#ifdef MultiNode
void updateInitAcc01(int update_count, int* update_count_list, int* displs);
void updateInitAcc23(int update_count, int* update_count_list, int* displs);
void sendNewNeighbors(int update_count, std::vector<int>& neighbors);
void updateTimeVariables(int update_count, int* update_count_list, int* displs);
void updateTimeCorrection(int update_count);
void updateIrregularForce(int update_count, int* update_count_list, int* displs);
void updateFBTermination(int update_count);
void updateNewGroup(int ptcl_id);
void updateAfterRegCuda(int update_count, int* update_count_list, int* displs);
void updateAfterRegCudaUpdate(int update_count, int* update_count_list, int* displs);
#ifdef SEVN
void updateStellarEvolution(int update_count);
#endif
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
	UpdateBinary* updateBinary_list;
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
			
			// In development by Minyong Jung
			case RegSend:
				int N_j, N_i;
				
				MPI_Recv(&N_i, 1, MPI_INT, ROOT, 1010, MPI_COMM_DEVICE, &status);
				MPI_Recv(&N_j, 1, MPI_INT, ROOT, 1011, MPI_COMM_DEVICE, &status);
				AllocateDeviceMemory(N_i, N_j, ptcl_id); // this part cannot be parallelized
				SendToDeviceMPI(N_i, N_j, ptcl_id); //this part should be parallelized
				// sendAllParticlesToGPU(next_time, RegularList, IndexList, N_start, N_end, ptcl_id);
				// predict and send particles to device
				break;
				
			case RegCal: {
				int* int_lists = new int[3];
				MPI_Recv(int_lists, 3, MPI_INT, ROOT, 1, MPI_COMM_DEVICE, &status);
				RegularWorker(int_lists[0], int_lists[1], int_lists[2], ptcl_id);
				delete int_lists;
				break;
			}

			case IrrUpdate: // Irregular Update Particle

				ptcl = &particles[ptcl_id];
				if (ptcl->NumberOfNeighbor != 0) // IAR modified
					ptcl->updateParticle();
				ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
				ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
				break;

			case RegUpdate: // Regular Update Particle

				ptcl = &particles[ptcl_id];
#ifndef MultiNode		
				std::memcpy(ptcl->Neighbors, ptcl->NewNeighbors, sizeof(int) * ptcl->NewNumberOfNeighbor);
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;
				ptcl->updateParticle();
#endif
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
#ifndef MultiNode
				std::memcpy(ptcl->Neighbors, ptcl->NewNeighbors, sizeof(int) * ptcl->NewNumberOfNeighbor);
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;
				ptcl->updateParticle();
#endif
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

			case UpdateIrregularForce:

				updateIrregularForce(ptcl_id, update_count_list, displs);
				break;

			case UpdateBinaryMerger:

				ptcl = &particles[ptcl_id];
				ptcl->Mass = -1.0;
				ptcl->CMPtclIndex = -1;
				break;

			case UpdateFBTermination:

				updateFBTermination(ptcl_id);
				break;

			case UpdateNewGroup:

				updateNewGroup(ptcl_id);
				break;

			case UpdateManybodyMerger: {

				ptcl = &particles[ptcl_id]; // This is zero mass particle
				Particle* ptclCM = &particles[ptcl->CMPtclIndex];
				for (int i = 0; i < ptclCM->NumberOfMember; i++) {
					Particle* member_ptcl = &particles[ptclCM->Members[i]];
					if (member_ptcl->PID == ptcl->PID) {
						member_ptcl->Mass = -1.0;
						member_ptcl->CMPtclIndex = -1;
						ptclCM->Members[i] = ptclCM->Members[ptclCM->NumberOfMember - 1];
						ptclCM->NumberOfMember--;
						break;
					}
				}
				break;
			}

			case UpdateActiveIndexToOriginalIndex:

				MPI_Recv(ActiveIndexToOriginalIndex, ptcl_id, MPI_INT, ROOT, 1, MPI_COMM_WORLD, &status);
				break;

			case UpdateBeforeRegCuda: {

				int *IndexList = new int[ptcl_id];
				int *NumNeighborReceive = new int[ptcl_id];
				int *ACListReceive = new int[ptcl_id * MaxNumNeighbor];
				CUDA_REAL (*AccRegReceive_f)[Dim] = new CUDA_REAL[ptcl_id][Dim];
				CUDA_REAL (*AccRegDotReceive_f)[Dim] = new CUDA_REAL[ptcl_id][Dim];
				MPI_Recv(IndexList, 			ptcl_id, 					MPI_INT, 	ROOT, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(NumNeighborReceive, 	ptcl_id, 					MPI_INT, 	ROOT, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(ACListReceive, 		ptcl_id * MaxNumNeighbor, 	MPI_INT, 	ROOT, 3, MPI_COMM_WORLD, &status);
				MPI_Recv(AccRegReceive_f, 		ptcl_id * Dim, 				MPI_FLOAT,	ROOT, 4, MPI_COMM_WORLD, &status);
				MPI_Recv(AccRegDotReceive_f, 	ptcl_id * Dim, 				MPI_FLOAT,	ROOT, 5, MPI_COMM_WORLD, &status);

				for (int i=0; i<ptcl_id; i++) {
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
				break;
			}

			case UpdateAfterRegCuda:

				updateAfterRegCuda(ptcl_id, update_count_list, displs);
				break;

			case UpdateAfterRegCudaUpdate:

				updateAfterRegCudaUpdate(ptcl_id, update_count_list, displs);
				break;
#ifdef SEVN
			case UpdateStellarEvolution1:

				updateStellarEvolution(ptcl_id);
				break;
#endif
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

#ifdef MultiNode
				updateBinary_list = new UpdateBinary[ptcl->NumberOfMember + 1];
				updateBinary_list[0].pid = ptcl->ParticleIndex;
				updateBinary_list[0].binary_state = ptcl->binary_state;
				std::memcpy(updateBinary_list[0].position, ptcl->Position, sizeof(double) * Dim);
				std::memcpy(updateBinary_list[0].velocity, ptcl->Velocity, sizeof(double) * Dim);
				updateBinary_list[0].mass = ptcl->Mass;
				updateBinary_list[0].currenttime_irr = ptcl->CurrentTimeIrr;

				for (int i = 0; i < ptcl->NumberOfMember; i++) {
					int member_index = ptcl->Members[i];
					Particle* member_ptcl = &particles[member_index];
					updateBinary_list[i + 1].pid = member_ptcl->ParticleIndex;
					updateBinary_list[i + 1].binary_state = member_ptcl->binary_state;
					std::memcpy(updateBinary_list[i + 1].position, member_ptcl->Position, sizeof(double) * Dim);
					std::memcpy(updateBinary_list[i + 1].velocity, member_ptcl->Velocity, sizeof(double) * Dim);
					updateBinary_list[i + 1].mass = member_ptcl->Mass;
					updateBinary_list[i + 1].currenttime_irr = member_ptcl->CurrentTimeIrr;
				}
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
		else if (task == ARIntegration) {
			MPI_Isend(updateBinary_list, ptcl->NumberOfMember + 1, UpdateBinaryType, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);
		}
		else if (task == UpdateLastParticleIndex || task == UpdateBinaryMerger || task == UpdateManybodyMerger)
			continue;
		else {
			MPI_Isend(&ptcl_id, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);
		}
#else
		MPI_Isend(&ptcl_id, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);
#endif

		MPI_Wait(&request, &status);
#ifdef MultiNode
		if (task == SendNewNeighbors) {
			neighbors.clear();
		}
		else if (task == ARIntegration) {
			delete[] updateBinary_list;
			updateBinary_list = nullptr;
		}
#endif
		//std::cerr << "Processor " << MyRank << " done." << std::endl;
	}
}