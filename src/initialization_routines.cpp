#include "queue_scheduler.h"

void formPrimordialBinaries(int beforeLastParticleIndex);
#ifdef SEVN_BINARY
bool makeSEVNBinary(Particle* ptclCM);
#endif
void InitialAssignmentOfTasks(Queue queue, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ull_t &data);
void broadcastFromRoot(int &data);

void InitializationRoutines(QueueScheduler &queue_scheduler, Worker *workers) {

    std::cout << "Initialization of particles starts." << std::endl;

    Particle* ptcl;
    Queue queue;
    task_name_t task;
    int total_tasks;
	int completed_tasks;
    MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object

#ifdef SEVN_BINARY
	std::chrono::high_resolution_clock::time_point start_point_BSE;
	std::chrono::high_resolution_clock::time_point end_point_BSE;
#endif

    std::vector<int> PIDs;
    PIDs.reserve(num_particles);
    PIDs.resize(num_particles);
    for (int i = 0; i <= last_particle_index; i++)
        PIDs[i] = i;

    queue_scheduler.initialize(TASK_INIT_ACC_1);
    queue_scheduler.takeQueue(PIDs);
    do {
        //queue_scheduler.printFreeWorker();
        //queue_scheduler.printWorkerToGo();
        queue_scheduler.assignQueueAuto();
        //queue_scheduler.printFreeWorker();
        //queue_scheduler.printWorkerToGo();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); //blocking wait
    } while(queue_scheduler.isComplete());
    std::cout << "Init 01 done" << std::endl;

    queue_scheduler.initialize(TASK_INIT_ACC_2);
    queue_scheduler.takeQueue(PIDs);
    do {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); //blocking wait
    } while(queue_scheduler.isComplete());
    std::cout << "Init 02 done" << std::endl;

#ifdef FEWBODY
    // Primordial binary search
    queue_scheduler.initialize(TASK_SEARCH_PRIMORDIAL_GROUP);
    queue_scheduler.takeQueue(PIDs);
    do {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); //blocking wait
    } while(queue_scheduler.isComplete());
    std::cout << "Primordial binary search done" << std::endl;

    int rank;
    int OriginalLastParticleIndex = last_particle_index;
    formPrimordialBinaries(OriginalLastParticleIndex);
    assert(OriginalLastParticleIndex <= last_particle_index); // for debugging by EW 2025.1.4
    assert(cm_particle_worker_map.empty()); // for debugging by EW 2025.1.4
    // Let's modify this primordial binary part later!!! by EW 2025.5.24
    if (OriginalLastParticleIndex != last_particle_index) {
        std::cout << "In total, " << last_particle_index - OriginalLastParticleIndex
                  << " primordial binaries are created." << std::endl;
        queue_scheduler.initialize(TASK_MAKE_PRIMORDIAL_GROUP);
        for (int i=OriginalLastParticleIndex+1; i<=last_particle_index; i++) {
            ptcl = &particles[i];
            cm_particle_worker_map.insert({ptcl->particle_index, cm_particle_worker_map.size() % num_workers + 1});
            PIDs.push_back(ptcl->particle_index);
            rank = cm_particle_worker_map[ptcl->particle_index];
            std::cout << "New Primordial Binary of PID="
                      << ptcl->pid << " is created with being assigned to a worker of rank "
                      << rank << "." << std::endl;
#ifdef SEVN_BINARY
#ifdef PERFORMANCETRACE
            start_point_BSE = std::chrono::high_resolution_clock::now();
#endif
            if (!makeSEVNBinary(ptcl)) {
                if (ptcl->particle_index == last_particle_index) {
                    last_particle_index--;
                    g_state->last_particle_index = last_particle_index;
                }
                else
                    prev_cm_particle_worker_map.insert({ptcl->particle_index, cm_particle_worker_map[ptcl->particle_index]});
                cm_particle_worker_map.erase(ptcl->particle_index);
                continue;
            }
#ifdef PERFORMANCETRACE
            end_point_BSE = std::chrono::high_resolution_clock::now();
            performance.BinaryStellarEvolution +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_BSE - start_point_BSE).count();
#endif
#endif
            queue.task = TASK_MAKE_PRIMORDIAL_GROUP;
            queue.pid = ptcl->particle_index;
            workers[rank].add_queue(queue);
            queue_scheduler.assignWorker(&workers[rank]);
        }
        queue_scheduler.setTotalQueue(cm_particle_worker_map.size());
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
            num_particles);
    fflush(stdout);
#endif

    // Initialize Time Step
    queue_scheduler.initialize(TASK_INIT_TIME);
    queue_scheduler.takeQueue(PIDs);
    do {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); //blocking wait
    } while(queue_scheduler.isComplete());

    /* Synchronization */
    //ParticleSynchronization();

    /* Timestep correction */
    int min_time_level=0;
    std::cout << "Time Step correction." << std::endl;
    for (int i=0; i<=last_particle_index; i++) {
        ptcl = &particles[i];

        if (!ptcl->is_active)
            continue;

        if (ptcl->num_neighbors != 0) {
            while (ptcl->time_level_irr >= ptcl->time_level_reg) {
                ptcl->time_step_irr *= 0.5;
                ptcl->time_block_irr *= 0.5;
                ptcl->time_level_irr--;
            }
        }
        if (ptcl->time_level_irr < min_time_level) {
            min_time_level = ptcl->time_level_irr;
        }
    }

    // resetting time_block based on the system
    time_block = std::max(-60, min_time_level-MIN_LEVEL_BUFFER);
    block_max = static_cast<ull_t>(pow(2, -time_block));
    time_step = pow(2,time_block);

    for (int i=0; i<=last_particle_index; i++) {
        ptcl = &particles[i];

        if (!ptcl->is_active)
            continue;

        ptcl->time_block_irr = static_cast<ull_t>(pow(2, ptcl->time_level_irr-time_block));
        ptcl->time_block_reg = static_cast<ull_t>(pow(2, ptcl->time_level_reg-time_block));
#ifdef IRR_TEST
        ptcl->time_step_reg = 1;
        ptcl->time_level_reg = 0;
        ptcl->time_block_reg = block_max;
#endif
        ptcl->next_block_irr = ptcl->current_block_irr + ptcl->time_block_irr; // of this particle
    }
    std::cout << "Time Step done." << std::endl;
    // Timestep correction ends

    /* Timestep variable synchronization */
    std::cout << "Time Step synchronization." << std::endl;
    task=TASK_TIME_SYNC;
    completed_tasks = 0; total_tasks = num_workers;
    queue = {task, -1, -1.0};
    InitialAssignmentOfTasks(queue, num_workers, QUEUE_TAG);
    broadcastFromRoot(time_block);
    broadcastFromRoot(block_max);
    broadcastFromRoot(time_step);
    fprintf(stdout, "TASK_TIME_SYNC broadcast done.\n");
    fflush(stdout);
    //MPI_Win_sync(win);  // TASK_SYNCHRONIZE memory
    //MPI_Barrier(shared_comm);
    int task_signal = 0; // Use int buffer for MPI_INT payloads.
    while (completed_tasks < total_tasks) {
        MPI_Irecv(&task_signal, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        completed_tasks++;
    }
    fprintf(stdout, "my_rank = %d time_block = %d, enzo_time_step = %e\n", my_rank, time_block, enzo_time_step);
    fflush(stdout);

    /* Particle Initialization Check */
    // /*
    for (int i=0; i<=last_particle_index; i++) {
        ptcl = &particles[i];
        if (ptcl->is_active)
            fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
                            "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
                            "NumNeighbor= %d\n",
                    ptcl->pid,
                    ptcl->current_time_irr * enzo_time_step * 1e10 / 1e6,
                    ptcl->current_block_irr,
                    ptcl->current_time_reg * enzo_time_step * 1e10 / 1e6,
                    ptcl->current_block_reg,
                    ptcl->time_step_irr * enzo_time_step * 1e10 / 1e6,
                    ptcl->time_step_reg * enzo_time_step * 1e10 / 1e6,
                    ptcl->time_block_irr,
                    ptcl->time_level_irr,
                    ptcl->time_block_reg,
                    ptcl->time_level_reg,
                    ptcl->num_neighbors);
            /*
            fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
                a1_reg = (%.4e,%.4e,%.4e), a2_reg = (%.4e,%.4e,%.4e), a3_reg = (%.4e,%.4e,%.4e)\n\
                a1_irr = (%.4e,%.4e,%.4e), a2_irr = (%.4e,%.4e,%.4e), a3_irr = (%.4e,%.4e,%.4e)\n", 
                ptcl->acc_total[0][0],	ptcl->acc_total[1][0],	ptcl->acc_total[2][0],
                ptcl->acc_regular[0][0],	ptcl->acc_regular[1][0],	ptcl->acc_regular[2][0],
                ptcl->acc_irregular[0][0],	ptcl->acc_irregular[1][0],	ptcl->acc_irregular[2][0],
                ptcl->num_neighbors,	ptcl->neighbor_radius_sq,
                ptcl->acc_regular[0][1],	ptcl->acc_regular[1][1],	ptcl->acc_regular[2][1],	
                ptcl->acc_regular[0][2],	ptcl->acc_regular[1][2],	ptcl->acc_regular[2][2],	
                ptcl->acc_regular[0][3],	ptcl->acc_regular[1][3],	ptcl->acc_regular[2][3],	
                ptcl->acc_irregular[0][1],	ptcl->acc_irregular[1][1],	ptcl->acc_irregular[2][1],
                ptcl->acc_irregular[0][2],	ptcl->acc_irregular[1][2],	ptcl->acc_irregular[2][2],
                ptcl->acc_irregular[0][3],	ptcl->acc_irregular[1][3],	ptcl->acc_irregular[2][3]);
            */
    }
}
