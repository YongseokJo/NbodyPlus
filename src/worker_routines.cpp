#include <iostream>
#include <vector>
#include <errno.h>
#include "global.h"
#include "queue.h"
#include "profiler.h"

void broadcastFromRoot(double &data);
void broadcastFromRoot(ull_t &data);
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

	//std::cout << "Processor " << my_rank << " is ready." << std::endl;

#ifdef PERFORMANCETRACE
	// Phase 21.5: Initialize worker tracking so compute times can be recorded
	// num_workers is a global variable set during MPI initialization
	if (num_workers > 0) {
		Profiler::instance().initializeWorkerTracking(num_workers);
	}
#endif

	task_name_t task = TASK_ERROR;
	MPI_Status status;
	MPI_Request request;
	int ptcl_id;
	double next_time;
	Particle *ptcl;
	Queue queue;

	while (true) {

		PROFILE_START(TimerID::WorkerRecvWait);
		MPI_Recv(&queue, 1, queue_type_mpi, ROOT, QUEUE_TAG, MPI_COMM_WORLD, &status);
		PROFILE_STOP(TimerID::WorkerRecvWait);
		task = queue.task;
		ptcl_id = queue.pid;
		next_time = queue.next_time;

		PROFILE_START(TimerID::WorkerTaskDispatch);
		switch (task) {
			case TASK_IRR_FORCE: // Irregular Acceleration

				ptcl = &particles[ptcl_id];
				ptcl->compute_acceleration_irr();

				ptcl->new_current_block_irr = ptcl->current_block_irr + ptcl->time_block_irr; // of this particle
				ptcl->calculate_time_step_irr();
				ptcl->next_block_irr = ptcl->new_current_block_irr + ptcl->time_block_irr; // of this particle
				ptcl->is_up_to_date = true;
				break;

			case TASK_REG_FORCE: // Regular Acceleration

				ptcl = &particles[ptcl_id];
				ptcl->compute_acceleration_reg();
				break;

			case TASK_IRR_UPDATE: // Irregular Update Particle

				ptcl = &particles[ptcl_id];
				if (ptcl->num_neighbors != 0) // IAR modified
					ptcl->update_particle();
				ptcl->current_block_irr = ptcl->new_current_block_irr;
				ptcl->current_time_irr  = ptcl->current_block_irr*time_step;
				break;

			case TASK_REG_UPDATE: // Regular Update Particle

				ptcl = &particles[ptcl_id];

				ptcl->current_block_reg += ptcl->time_block_reg;
				ptcl->current_time_reg = ptcl->current_block_reg * time_step;

				ptcl->update_particle();
				std::memcpy(neighbors + ptcl->neighbors_offset, new_neighbors + ptcl->neighbors_offset, sizeof(int) * ptcl->new_num_neighbors);
				ptcl->num_neighbors = ptcl->new_num_neighbors;

				ptcl->calculate_time_step_reg();
				ptcl->calculate_time_step_irr();
				// /*
				if (ptcl->current_block_irr != ptcl->current_block_reg || ptcl->current_time_irr != ptcl->current_time_reg) {
					fprintf(stderr, "WARNING!!! In RegCudaUpdate...\n");
					fprintf(stderr, "PID: %d\n", ptcl->pid);
					fprintf(stderr, "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n", ptcl->current_block_irr, ptcl->current_block_reg);
					fprintf(stderr, "TimeBlockIrr: %llu, TimeBlockReg: %llu\n", ptcl->time_block_irr, ptcl->time_block_reg);
					fprintf(stderr, "CurrentBlockIrr * time_step: %.17g, CurrentBlockReg * time_step: %.17g\n", ptcl->current_block_irr*time_step, ptcl->current_block_reg*time_step);
					fprintf(stderr, "CurrentTimeIrr: %.17g, CurrentTimeReg: %.17g\n", ptcl->current_time_irr, ptcl->current_time_reg);
					fprintf(stderr, "next_reg_time_block: %llu\n", g_state->next_reg_time_block);
					fflush(stderr);
					assert(ptcl->current_block_irr + ptcl->time_block_irr > ptcl->current_block_reg);
					// assert(ptcl->current_time_irr == ptcl->current_time_reg);
					// assert(ptcl->current_block_irr == ptcl->current_block_reg);
				}
				// */
				ptcl->update_radius();
				ptcl->next_block_irr = ptcl->current_block_irr + ptcl->time_block_irr; // of ptcl particle
				break;

			case TASK_REG_CUDA: // Regular Correct Particle After CUDA

				ptcl = &particles[ptcl_id];
				ptcl->update_regular_particle_cuda();
				break;

			case TASK_INIT_ACC_1: // Initialize Acceleration(01)

				ptcl = &particles[ptcl_id];
				CalculateAcceleration01(ptcl);
				break;

			case TASK_INIT_ACC_2: // Initialize Acceleration(23)

				ptcl = &particles[ptcl_id];
				CalculateAcceleration23(ptcl);
				break;

			case TASK_INIT_TIME: // Initialize Time Step

				ptcl = &particles[ptcl_id];
				if (ptcl->is_active)
					ptcl->initialize_time_step();
				break;

			case TASK_TIME_SYNC: // Initialize Timestep variables
				broadcastFromRoot(time_block);
				broadcastFromRoot(block_max);
				broadcastFromRoot(time_step);
				//MPI_Win_sync(win);  // TASK_SYNCHRONIZE memory
				//MPI_Barrier(shared_comm);
				//MPI_Win_fence(0, win);
				fprintf(worker_output_file, "my_rank = %d time_block = %d, enzo_time_step = %e\n\n", my_rank, time_block, enzo_time_step);
				break;

#ifdef FEWBODY
			case TASK_SEARCH_PRIMORDIAL_GROUP: // Primordial binary search

				ptcl = &particles[ptcl_id];

				ptcl->new_num_members = 0;
				ptcl->check_new_group_v2();
				break;

			case TASK_SEARCH_GROUP: // Few-body group search

				ptcl = &particles[ptcl_id];

				if (ptcl->get_binary_interrupt_state()==BinaryInterruptState::threebody) {
					ptcl->set_binary_interrupt_state(BinaryInterruptState::none);
				}
				else if (ptcl->get_binary_interrupt_state()==BinaryInterruptState::manybody) {
					ptcl->new_num_members = 0;
					ptcl->check_new_group_v2();
					ptcl->set_binary_interrupt_state(BinaryInterruptState::none);
				}
				else {
					ptcl->new_num_members = 0;
					if (ptcl->time_step_irr < t_search)
						ptcl->check_new_group();
				}
				/*
				if (ptcl->get_binary_interrupt_state()==BinaryInterruptState::manybody) {
					ptcl->set_binary_interrupt_state(BinaryInterruptState::none);
					std::cout << "ptcl PID: " << ptcl->pid << ", ptcl NewNumberOfMember: " << ptcl->new_num_members << std::endl;
				}
				else {
					ptcl->new_num_members = 0;
					if (ptcl->time_step_irr < t_search)
						ptcl->check_new_group();
				}
				*/
				break;

			case TASK_MAKE_PRIMORDIAL_GROUP: // Make a primordial group

				ptcl = &particles[ptcl_id];
				makePrimordialGroup(ptcl);
				break;

			case TASK_MAKE_GROUP: // Make a group

				ptcl = &particles[ptcl_id];

				NewFBInitialization(ptcl);
#ifdef DEBUG
				std::cout << "FewBody object of particle " << ptcl->pid
						  << " is successfully initialized on rank " << my_rank << "." <<std::endl;
#endif
				break;

			case TASK_DELETE_GROUP: // Delete a Group struct

				ptcl = &particles[ptcl_id];
				deleteGroup(ptcl);
				break;

			case TASK_AR_INTEGRATION: // SDAR for few body encounters
				
				ptcl = &particles[ptcl_id];

				if (!ptcl->is_cm_particle || ptcl->group_info == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->is_cm_particle=%d ptcl->group_info=%p\n", ptcl->is_cm_particle, ptcl->group_info);
					exit(EXIT_FAILURE);
				}
				
				ptcl->group_info->ARIntegration(next_time);
				if (!ptcl->group_info->isMerger && !ptcl->group_info->isTerminate)
					ptcl->group_info->isTerminate = ptcl->group_info->CheckBreak2();

				if (ptcl->group_info->isTerminate) {
					if (ptcl->get_binary_interrupt_state() == BinaryInterruptState::none)
						ptcl->set_binary_interrupt_state(BinaryInterruptState::terminated);

					energy_binary -= ptcl->group_info->sym_int.getEtot();
					energy_binary_sd -= ptcl->group_info->sym_int.getEtotSlowDown();

					delete ptcl->group_info;
#ifdef DEBUG
					std::cout << "(SDAR) Processor " << my_rank<< ": PID= "<<ptcl->pid << " deleted!" <<std::endl;
#endif
				}
#ifdef DEBUG
				else
					std::cout << "(SDAR) Processor " << my_rank<< ": PID= "<<ptcl->pid << " done!" <<std::endl;
#endif
				break;
			
			case TASK_MERGE_MANYBODY: // Merger insided many-body (>2) group
				
				ptcl = &particles[ptcl_id];
				std::cout << "(SDAR) Processor " << my_rank<< ": PID= "<<ptcl->pid << std::endl;

				if (!ptcl->is_cm_particle || ptcl->group_info == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->is_cm_particle=%d ptcl->group_info=%p\n", ptcl->is_cm_particle, ptcl->group_info);
					exit(EXIT_FAILURE);
				}

				NewFBInitialization3(ptcl->group_info);

				ptcl->group_info->isMerger = false;
				ptcl->set_binary_interrupt_state(BinaryInterruptState::none);

				std::cout << "(SDAR) Processor " << my_rank<< ": PID= "<<ptcl->pid << " NewFBInitialization3 done!" <<std::endl;
				break;
#endif 

			case TASK_GET_TOTAL_ENERGY:

				MPI_Reduce(&energy_binary,		nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&energy_binary_sd,	nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&energy_merger,		nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				MPI_Reduce(&energy_pn,			nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
				continue;
#ifdef CUDA
			case TASK_PREPARE_GPU_CALC:

				sendAllParticlesToGPU_Worker(next_time);
				continue;
#endif
			case TASK_SEND_PROFILING: // Phase 21.5: Send profiling data to root
				{
					ProfilerTransferData pdata = profiler().packTransferData(my_rank);
					MPI_Send(&pdata, sizeof(ProfilerTransferData), MPI_BYTE, ROOT, PROFILING_TAG, MPI_COMM_WORLD);
				}
				continue;  // Don't send TERMINATE_TAG, root is waiting for profiling data

			case TASK_SYNCHRONIZE: // TASK_SYNCHRONIZE
				MPI_Win_sync(win);  // TASK_SYNCHRONIZE memory
				MPI_Barrier(shared_comm);
				break;

			case TASK_END: // Simualtion ends
				fprintf(worker_output_file, "Processor %d returns.\n", my_rank);
				return;

			case TASK_ERROR:
				perror("TASK_ERROR task assignments");
				exit(EXIT_FAILURE);
				break;

			default:
				break;
		}
		PROFILE_STOP(TimerID::WorkerTaskDispatch);

		// return that it's over
		//task = -1;
		PROFILE_START(TimerID::WorkerSendComplete);
		int send_value = ptcl_id;
		if (!(task == TASK_IRR_FORCE || task == TASK_REG_FORCE || task == TASK_IRR_UPDATE || task == TASK_REG_UPDATE))
			send_value = static_cast<int>(task);
		MPI_Isend(&send_value, 1, MPI_INT, ROOT, TERMINATE_TAG, MPI_COMM_WORLD, &request);

		MPI_Wait(&request, &status);
		PROFILE_STOP(TimerID::WorkerSendComplete);
		//std::cerr << "Processor " << my_rank << " done." << std::endl;
	}
}
