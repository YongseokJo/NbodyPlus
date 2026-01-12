#include "QueueScheduler.h"

void formPrimordialBinaries(int beforeLastParticleIndex);
#ifdef SEVN_BINARY
bool makeSEVNBinary(Particle* ptclCM);
#endif
void InitialAssignmentOfTasks(Queue queue, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);

void InitializationRoutines(QueueScheduler &queue_scheduler, Worker *workers) {

    std::cout << "Initialization of particles starts." << std::endl;

    Particle* ptcl;
    Queue queue;
    TaskName task;
    int total_tasks;
	int completed_tasks;
    MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object

#ifdef SEVN_BINARY
	std::chrono::high_resolution_clock::time_point start_point_BSE;
	std::chrono::high_resolution_clock::time_point end_point_BSE;
#endif

    std::vector<int> PIDs;
    PIDs.reserve(NumberOfParticle);
    PIDs.resize(NumberOfParticle);
    for (int i = 0; i <= LastParticleIndex; i++)
        PIDs[i] = i;

    queue_scheduler.initialize(InitAcc1);
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

    queue_scheduler.initialize(InitAcc2);
    queue_scheduler.takeQueue(PIDs);
    do {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); //blocking wait
    } while(queue_scheduler.isComplete());
    std::cout << "Init 02 done" << std::endl;

#ifdef FEWBODY
    // Primordial binary search
    queue_scheduler.initialize(SearchPrimordialGroup);
    queue_scheduler.takeQueue(PIDs);
    do {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); //blocking wait
    } while(queue_scheduler.isComplete());
    std::cout << "Primordial binary search done" << std::endl;

    int rank;
    int OriginalLastParticleIndex = LastParticleIndex;
    formPrimordialBinaries(OriginalLastParticleIndex);
    assert(OriginalLastParticleIndex <= LastParticleIndex); // for debugging by EW 2025.1.4
    assert(CMPtclWorker.empty()); // for debugging by EW 2025.1.4
    // Let's modify this primordial binary part later!!! by EW 2025.5.24
    if (OriginalLastParticleIndex != LastParticleIndex) {
        std::cout << "In total, " << LastParticleIndex - OriginalLastParticleIndex
                  << " primordial binaries are created." << std::endl;
        queue_scheduler.initialize(MakePrimordialGroup);
        for (int i=OriginalLastParticleIndex+1; i<=LastParticleIndex; i++) {
            ptcl = &particles[i];
            CMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker.size() % NumberOfWorker + 1});
            PIDs.push_back(ptcl->ParticleIndex);
            rank = CMPtclWorker[ptcl->ParticleIndex];
            std::cout << "New Primordial Binary of PID="
                      << ptcl->PID << " is created with being assigned to a worker of rank "
                      << rank << "." << std::endl;
#ifdef SEVN_BINARY
#ifdef PERFORMANCETRACE
            start_point_BSE = std::chrono::high_resolution_clock::now();
#endif
            if (!makeSEVNBinary(ptcl)) {
                if (ptcl->ParticleIndex == LastParticleIndex) {
                    LastParticleIndex--;
                    global_variable->LastParticleIndex = LastParticleIndex;
                }
                else
                    PrevCMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker[ptcl->ParticleIndex]});
                CMPtclWorker.erase(ptcl->ParticleIndex);
                continue;
            }
#ifdef PERFORMANCETRACE
            end_point_BSE = std::chrono::high_resolution_clock::now();
            performance.BinaryStellarEvolution +=
                std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_BSE - start_point_BSE).count();
#endif
#endif
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
    // Timestep correction ends

    /* Timestep variable synchronization */
    std::cout << "Time Step synchronization." << std::endl;
    task=TimeSync;
    completed_tasks = 0; total_tasks = NumberOfWorker;
    queue = {task, -1, -1.0};
    InitialAssignmentOfTasks(queue, NumberOfWorker, QUEUE_TAG);
    broadcastFromRoot(time_block);
    broadcastFromRoot(block_max);
    broadcastFromRoot(time_step);
    fprintf(stdout, "TimeSync broadcast done.\n");
    fflush(stdout);
    //MPI_Win_sync(win);  // Synchronize memory
    //MPI_Barrier(shared_comm);
    while (completed_tasks < total_tasks) {
        int recv_task = 0; // TaskName is int8_t; receive into int to match MPI_INT.
        int mpi_rc = MPI_Irecv(&recv_task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &request);
        if (mpi_rc != MPI_SUCCESS) {
            char err[MPI_MAX_ERROR_STRING];
            int err_len = 0;
            MPI_Error_string(mpi_rc, err, &err_len);
            fprintf(stderr, "MPI_Irecv failed in TimeSync (rc=%d): %s\n", mpi_rc, err);
            MPI_Abort(MPI_COMM_WORLD, mpi_rc);
        }
        mpi_rc = MPI_Wait(&request, &status);
        if (mpi_rc != MPI_SUCCESS || status.MPI_ERROR != MPI_SUCCESS) {
            char err[MPI_MAX_ERROR_STRING];
            int err_len = 0;
            int status_rc = (mpi_rc != MPI_SUCCESS) ? mpi_rc : status.MPI_ERROR;
            MPI_Error_string(status_rc, err, &err_len);
            fprintf(stderr, "MPI_Wait failed in TimeSync (rc=%d): %s\n", status_rc, err);
            MPI_Abort(MPI_COMM_WORLD, status_rc);
        }
        completed_tasks++;
    }
    fprintf(stdout, "MyRank = %d time_block = %d, EnzoTimeStep = %e\n", MyRank, time_block, EnzoTimeStep);
    fflush(stdout);

    /* Particle Initialization Check */
    // /*
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
                ptcl->a_tot[0][0],	ptcl->a_tot[1][0],	ptcl->a_tot[2][0],
                ptcl->a_reg[0][0],	ptcl->a_reg[1][0],	ptcl->a_reg[2][0],
                ptcl->a_irr[0][0],	ptcl->a_irr[1][0],	ptcl->a_irr[2][0],
                ptcl->NumberOfNeighbor,	ptcl->RadiusOfNeighbor,
                ptcl->a_reg[0][1],	ptcl->a_reg[1][1],	ptcl->a_reg[2][1],	
                ptcl->a_reg[0][2],	ptcl->a_reg[1][2],	ptcl->a_reg[2][2],	
                ptcl->a_reg[0][3],	ptcl->a_reg[1][3],	ptcl->a_reg[2][3],	
                ptcl->a_irr[0][1],	ptcl->a_irr[1][1],	ptcl->a_irr[2][1],
                ptcl->a_irr[0][2],	ptcl->a_irr[1][2],	ptcl->a_irr[2][2],
                ptcl->a_irr[0][3],	ptcl->a_irr[1][3],	ptcl->a_irr[2][3]);
            */
    }
}
