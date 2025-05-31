#ifdef MultiNode

#include "global.h"
#include "Queue.h"
#include <cstring>

void updateInitAcc01(int update_count, int* update_count_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(update_count);
    MPI_Recv(update_pid_list.data(), update_count, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&update_count, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);
    // fprintf(workerout, "1. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

    UpdateInitAcc1* update_list = new UpdateInitAcc1[update_count_list[update_rank]];
    std::vector<int> neighbors;
    neighbors.reserve(MaxNumNeighbor * update_count);

    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_list[i].pid];
        update_list[i].numberofneighbors = ptcl->NumberOfNeighbor;
        update_list[i].airr0[0] = ptcl->a_irr[0][0];
        update_list[i].airr0[1] = ptcl->a_irr[1][0];
        update_list[i].airr0[2] = ptcl->a_irr[2][0];
        update_list[i].airr1[0] = ptcl->a_irr[0][1];
        update_list[i].airr1[1] = ptcl->a_irr[1][1];
        update_list[i].airr1[2] = ptcl->a_irr[2][1];
        update_list[i].areg0[0] = ptcl->a_reg[0][0];
        update_list[i].areg0[1] = ptcl->a_reg[1][0];
        update_list[i].areg0[2] = ptcl->a_reg[2][0];
        update_list[i].areg1[0] = ptcl->a_reg[0][1];
        update_list[i].areg1[1] = ptcl->a_reg[1][1];
        update_list[i].areg1[2] = ptcl->a_reg[2][1];
        // fprintf(workerout, "1. Gathering... MyRank: %d. PID: %d, NN: %d. airr0: %e %e %e, airr1: %e %e %e, areg0: %e %e %e, areg1: %e %e %e\n",
        //     MyRank, ptcl->PID, ptcl->NumberOfNeighbor,
        //     ptcl->a_irr[0][0], ptcl->a_irr[1][0], ptcl->a_irr[2][0],
        //     ptcl->a_irr[0][1], ptcl->a_irr[1][1], ptcl->a_irr[2][1],
        //     ptcl->a_reg[0][0], ptcl->a_reg[1][0], ptcl->a_reg[2][0],
        //     ptcl->a_reg[0][1], ptcl->a_reg[1][1], ptcl->a_reg[2][1]);
        // fprintf(workerout, "(");
        // for (int j=0; j<ptcl->NumberOfNeighbor; j++)
        //     fprintf(workerout, "%d ", ptcl->Neighbors[j]);
        // fprintf(workerout, ")\n");

        neighbors.insert(neighbors.end(), ptcl->Neighbors, ptcl->Neighbors + ptcl->NumberOfNeighbor);
    }

    int neighbor_count = neighbors.size();
    int* neighbor_count_list = new int[NumberOfNode];
    MPI_Allgather(&neighbor_count, 1, MPI_INT, neighbor_count_list, 1, MPI_INT, update_comm);
    int total_neighbor_count = 0;
    int* displs_neighbors = new int[NumberOfNode];

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];

        displs_neighbors[i] = total_neighbor_count;
        total_neighbor_count += neighbor_count_list[i];
    }
    UpdateInitAcc1* total_update_list = new UpdateInitAcc1[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateInitAcc1Type, total_update_list, update_count_list, displs, UpdateInitAcc1Type, update_comm);
    int* total_neighbors = new int[total_neighbor_count];
    MPI_Allgatherv(neighbors.data(), neighbors.size(), MPI_INT, total_neighbors, neighbor_count_list, displs_neighbors, MPI_INT, update_comm);
    // fprintf(workerout, "1. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);

    int offset = 0;
    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank]) {
            ptcl = &particles[total_update_list[i].pid];
            offset += ptcl->NumberOfNeighbor;
            continue;
        }

        ptcl = &particles[total_update_list[i].pid];
        ptcl->NumberOfNeighbor = total_update_list[i].numberofneighbors;
        ptcl->a_irr[0][0] = total_update_list[i].airr0[0];
        ptcl->a_irr[1][0] = total_update_list[i].airr0[1];
        ptcl->a_irr[2][0] = total_update_list[i].airr0[2];
        ptcl->a_irr[0][1] = total_update_list[i].airr1[0];
        ptcl->a_irr[1][1] = total_update_list[i].airr1[1];
        ptcl->a_irr[2][1] = total_update_list[i].airr1[2];
        ptcl->a_reg[0][0] = total_update_list[i].areg0[0];
        ptcl->a_reg[1][0] = total_update_list[i].areg0[1];
        ptcl->a_reg[2][0] = total_update_list[i].areg0[2];
        ptcl->a_reg[0][1] = total_update_list[i].areg1[0];
        ptcl->a_reg[1][1] = total_update_list[i].areg1[1];
        ptcl->a_reg[2][1] = total_update_list[i].areg1[2];
        ptcl->a_tot[0][0] = ptcl->a_irr[0][0] + ptcl->a_reg[0][0];
        ptcl->a_tot[1][0] = ptcl->a_irr[1][0] + ptcl->a_reg[1][0];
        ptcl->a_tot[2][0] = ptcl->a_irr[2][0] + ptcl->a_reg[2][0];
        ptcl->a_tot[0][1] = ptcl->a_irr[0][1] + ptcl->a_reg[0][1];
        ptcl->a_tot[1][1] = ptcl->a_irr[1][1] + ptcl->a_reg[1][1];
        ptcl->a_tot[2][1] = ptcl->a_irr[2][1] + ptcl->a_reg[2][1];
        std::memcpy(ptcl->Neighbors, &total_neighbors[offset], sizeof(int) * ptcl->NumberOfNeighbor);
        offset += ptcl->NumberOfNeighbor;

        // fprintf(workerout, "1. Receiving... MyRank: %d. PID: %d, NN: %d, airr0: %e %e %e, airr1: %e %e %e, areg0: %e %e %e, areg1: %e %e %e\n",
        //     MyRank, ptcl->PID, ptcl->NumberOfNeighbor,
        //     ptcl->a_irr[0][0], ptcl->a_irr[1][0], ptcl->a_irr[2][0],
        //     ptcl->a_irr[0][1], ptcl->a_irr[1][1], ptcl->a_irr[2][1],
        //     ptcl->a_reg[0][0], ptcl->a_reg[1][0], ptcl->a_reg[2][0],
        //     ptcl->a_reg[0][1], ptcl->a_reg[1][1], ptcl->a_reg[2][1]);
        // fprintf(workerout, "(");
        // for (int j=0; j<ptcl->NumberOfNeighbor; j++)
        //     fprintf(workerout, "%d ", ptcl->Neighbors[j]);
        // fprintf(workerout, ")\n");
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
    delete [] neighbor_count_list;
    neighbor_count_list = nullptr;
    delete [] displs_neighbors;
    displs_neighbors = nullptr;
    delete [] total_neighbors;
    total_neighbors = nullptr;
}

void updateInitAcc23(int update_count, int* update_count_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(update_count);
    MPI_Recv(update_pid_list.data(), update_count, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&update_count, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);
    // fprintf(workerout, "2. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

    UpdateInitAcc2* update_list = new UpdateInitAcc2[update_count_list[update_rank]];
    
    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_list[i].pid];
        update_list[i].airr2[0] = ptcl->a_irr[0][2];
        update_list[i].airr2[1] = ptcl->a_irr[1][2];
        update_list[i].airr2[2] = ptcl->a_irr[2][2];
        update_list[i].airr3[0] = ptcl->a_irr[0][3];
        update_list[i].airr3[1] = ptcl->a_irr[1][3];
        update_list[i].airr3[2] = ptcl->a_irr[2][3];
        update_list[i].areg2[0] = ptcl->a_reg[0][2];
        update_list[i].areg2[1] = ptcl->a_reg[1][2];
        update_list[i].areg2[2] = ptcl->a_reg[2][2];
        update_list[i].areg3[0] = ptcl->a_reg[0][3];
        update_list[i].areg3[1] = ptcl->a_reg[1][3];
        update_list[i].areg3[2] = ptcl->a_reg[2][3];
        // fprintf(workerout, "2. Gathering... MyRank: %d. PID: %d, airr2: %e %e %e, airr3: %e %e %e, areg2: %e %e %e, areg3: %e %e %e\n",
        //     MyRank, ptcl->PID,
        //     ptcl->a_irr[0][2], ptcl->a_irr[1][2], ptcl->a_irr[2][2],
        //     ptcl->a_irr[0][3], ptcl->a_irr[1][3], ptcl->a_irr[2][3],
        //     ptcl->a_reg[0][2], ptcl->a_reg[1][2], ptcl->a_reg[2][2],
        //     ptcl->a_reg[0][3], ptcl->a_reg[1][3], ptcl->a_reg[2][3]);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    UpdateInitAcc2* total_update_list = new UpdateInitAcc2[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateInitAcc2Type, total_update_list, update_count_list, displs, UpdateInitAcc2Type, update_comm);
    // fprintf(workerout, "2. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank])
            continue;

        ptcl = &particles[total_update_list[i].pid];
        ptcl->a_irr[0][2] = total_update_list[i].airr2[0];
        ptcl->a_irr[1][2] = total_update_list[i].airr2[1];
        ptcl->a_irr[2][2] = total_update_list[i].airr2[2];
        ptcl->a_irr[0][3] = total_update_list[i].airr3[0];
        ptcl->a_irr[1][3] = total_update_list[i].airr3[1];
        ptcl->a_irr[2][3] = total_update_list[i].airr3[2];
        ptcl->a_reg[0][2] = total_update_list[i].areg2[0];
        ptcl->a_reg[1][2] = total_update_list[i].areg2[1];
        ptcl->a_reg[2][2] = total_update_list[i].areg2[2];
        ptcl->a_reg[0][3] = total_update_list[i].areg3[0];
        ptcl->a_reg[1][3] = total_update_list[i].areg3[1];
        ptcl->a_reg[2][3] = total_update_list[i].areg3[2];
        ptcl->a_tot[0][2] = ptcl->a_irr[0][2] + ptcl->a_reg[0][2];
        ptcl->a_tot[1][2] = ptcl->a_irr[1][2] + ptcl->a_reg[1][2];
        ptcl->a_tot[2][2] = ptcl->a_irr[2][2] + ptcl->a_reg[2][2];
        ptcl->a_tot[0][3] = ptcl->a_irr[0][3] + ptcl->a_reg[0][3];
        ptcl->a_tot[1][3] = ptcl->a_irr[1][3] + ptcl->a_reg[1][3];
        ptcl->a_tot[2][3] = ptcl->a_irr[2][3] + ptcl->a_reg[2][3];
        // fprintf(workerout, "2. Receiving... MyRank: %d. PID: %d, airr2: %e %e %e, airr3: %e %e %e, areg2: %e %e %e, areg3: %e %e %e\n",
        //     MyRank, ptcl->PID,
        //     ptcl->a_irr[0][2], ptcl->a_irr[1][2], ptcl->a_irr[2][2],
        //     ptcl->a_irr[0][3], ptcl->a_irr[1][3], ptcl->a_irr[2][3],
        //     ptcl->a_reg[0][2], ptcl->a_reg[1][2], ptcl->a_reg[2][2],
        //     ptcl->a_reg[0][3], ptcl->a_reg[1][3], ptcl->a_reg[2][3]);
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
}

void sendNewNeighbors(int update_count, std::vector<int>& neighbors) {
    
    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(update_count);
    MPI_Recv(update_pid_list.data(), update_count, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // fprintf(workerout, "3. After Recv... MyRank: %d. list_size: %d\n", MyRank, update_count);

    neighbors.reserve(10 * update_count);

    for (int i=0; i<update_count; i++) {
        ptcl = &particles[update_pid_list[i]];
        neighbors.insert(neighbors.end(), ptcl->NewNeighbors, ptcl->NewNeighbors + ptcl->NewNumberOfNeighbor);
    }
}

void updateTimeVariables(int update_count, int* update_count_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(update_count);
    MPI_Recv(update_pid_list.data(), update_count, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&update_count, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);
    // fprintf(workerout, "4. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

    UpdateTime* update_list = new UpdateTime[update_count_list[update_rank]];

    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_pid_list[i]];
        update_list[i].timestep_reg = ptcl->TimeStepReg;
        update_list[i].timeblock_reg = ptcl->TimeBlockReg;
        update_list[i].timelevel_reg = ptcl->TimeLevelReg;
        update_list[i].timestep_irr = ptcl->TimeStepIrr;
        update_list[i].timeblock_irr = ptcl->TimeBlockIrr;
        update_list[i].timelevel_irr = ptcl->TimeLevelIrr;
        update_list[i].currenttime_irr = ptcl->CurrentTimeIrr;
        update_list[i].currenttime_reg = ptcl->CurrentTimeReg;
        update_list[i].currentblock_irr = ptcl->CurrentBlockIrr;
        update_list[i].currentblock_reg = ptcl->CurrentBlockReg;
        // fprintf(workerout, "4. Gathering... MyRank: %d. PID: %d, TimeStepReg: %e, TimeBlockReg: %llu, TimeLevelReg: %d, "
        //     "TimeStepIrr: %e, TimeBlockIrr: %llu, TimeLevelIrr: %d, CurrentTimeIrr: %e, CurrentTimeReg: %e, "
        //     "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n",
        //     MyRank, ptcl->PID,
        //     ptcl->TimeStepReg, ptcl->TimeBlockReg, ptcl->TimeLevelReg,
        //     ptcl->TimeStepIrr, ptcl->TimeBlockIrr, ptcl->TimeLevelIrr,
        //     ptcl->CurrentTimeIrr, ptcl->CurrentTimeReg, ptcl->CurrentBlockIrr, ptcl->CurrentBlockReg);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    UpdateTime* total_update_list = new UpdateTime[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateTimeType, total_update_list, update_count_list, displs, UpdateTimeType, update_comm);
    // fprintf(workerout, "4. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank])
            continue;

        ptcl = &particles[total_update_list[i].pid];
        ptcl->TimeStepReg = total_update_list[i].timestep_reg;
        ptcl->TimeBlockReg = total_update_list[i].timeblock_reg;
        ptcl->TimeLevelReg = total_update_list[i].timelevel_reg;
        ptcl->TimeStepIrr = total_update_list[i].timestep_irr;
        ptcl->TimeBlockIrr = total_update_list[i].timeblock_irr;
        ptcl->TimeLevelIrr = total_update_list[i].timelevel_irr;
        ptcl->CurrentTimeIrr = total_update_list[i].currenttime_irr;
        ptcl->CurrentTimeReg = total_update_list[i].currenttime_reg;
        ptcl->CurrentBlockIrr = total_update_list[i].currentblock_irr;
        ptcl->CurrentBlockReg = total_update_list[i].currentblock_reg;
        // fprintf(workerout, "4. Receiving... MyRank: %d. PID: %d, TimeStepReg: %e, TimeBlockReg: %llu, TimeLevelReg: %d, "
        //     "TimeStepIrr: %e, TimeBlockIrr: %llu, TimeLevelIrr: %d, CurrentTimeIrr: %e, CurrentTimeReg: %e, "
        //     "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n",
        //     MyRank, ptcl->PID,
        //     ptcl->TimeStepReg, ptcl->TimeBlockReg, ptcl->TimeLevelReg,
        //     ptcl->TimeStepIrr, ptcl->TimeBlockIrr, ptcl->TimeLevelIrr,
        //     ptcl->CurrentTimeIrr, ptcl->CurrentTimeReg, ptcl->CurrentBlockIrr, ptcl->CurrentBlockReg);
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
}

void updateTimeCorrection(int update_count) {

    Particle* ptcl;

    UpdateTimeCorr* total_update_list = new UpdateTimeCorr[update_count];

    MPI_Recv(total_update_list, update_count, UpdateTimeCorrType, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i < update_count; ++i) {
        ptcl = &particles[total_update_list[i].pid];
        ptcl->TimeStepIrr = total_update_list[i].timestep_irr;
        ptcl->TimeBlockIrr = total_update_list[i].timeblock_irr;
        ptcl->TimeLevelIrr = total_update_list[i].timelevel_irr;

        ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg-time_block));
        ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr;
    }
    delete[] total_update_list;
    total_update_list = nullptr;
}

void updateIrregularForce(int update_count, int* update_count_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(update_count);
    MPI_Recv(update_pid_list.data(), update_count, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&update_count, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);
    // fprintf(workerout, "5. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);
    // fflush(workerout);

    UpdateIrrForce* update_list = new UpdateIrrForce[update_count_list[update_rank]];
    
    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_list[i].pid];
        update_list[i].newnumberofneighbors = ptcl->NewNumberOfNeighbor;
        std::memcpy(update_list[i].newneighbors, ptcl->NewNeighbors, sizeof(int) * ptcl->NewNumberOfNeighbor);
        std::memcpy(update_list[i].newposition, ptcl->NewPosition, sizeof(double) * 3);
        std::memcpy(update_list[i].newvelocity, ptcl->NewVelocity, sizeof(double) * 3);
        for (int dim = 0; dim < Dim; dim++) {
            for (int j = 0; j < HERMITE_ORDER; j++) {
                update_list[i].airr[dim][j] = ptcl->a_irr[dim][j];
                update_list[i].atot[dim][j] = ptcl->a_tot[dim][j];
            }
        }
        update_list[i].newcurrentblock_irr = ptcl->NewCurrentBlockIrr;
        update_list[i].timelevel_irr = ptcl->TimeLevelIrr;
        update_list[i].timestep_irr = ptcl->TimeStepIrr;
        update_list[i].timeblock_irr = ptcl->TimeBlockIrr;
        update_list[i].nextblock_irr = ptcl->NextBlockIrr;
        // fprintf(workerout, "5. Gathering... MyRank: %d. PID: %d, NewNN: %d, NewPosition: %e %e %e, NewVelocity: %e %e %e, "
        //     "NewCurrentBlockIrr: %llu, TimeLevelIrr: %d, TimeStepIrr: %e, TimeBlockIrr: %llu, NextBlockIrr: %llu\n",
        //     MyRank, ptcl->PID, ptcl->NewNumberOfNeighbor,
        //     ptcl->NewPosition[0], ptcl->NewPosition[1], ptcl->NewPosition[2],
        //     ptcl->NewVelocity[0], ptcl->NewVelocity[1], ptcl->NewVelocity[2],
        //     ptcl->NextBlockIrr, ptcl->TimeLevelIrr, ptcl->TimeStepIrr,
        //     ptcl->TimeBlockIrr, ptcl->NextBlockIrr);
        // fflush(workerout);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    UpdateIrrForce* total_update_list = new UpdateIrrForce[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateIrrForceType, total_update_list, update_count_list, displs, UpdateIrrForceType, update_comm);
    // fprintf(workerout, "5. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
    // fflush(workerout);

    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank]) {
            ptcl = &particles[total_update_list[i].pid];
            
            // UpdateIrr Routine Starts
            if (ptcl->NumberOfNeighbor != 0) // IAR modified
					ptcl->updateParticle();
            ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
            ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
            continue;
        }

        ptcl = &particles[total_update_list[i].pid];
        ptcl->NewNumberOfNeighbor = total_update_list[i].newnumberofneighbors;
        std::memcpy(ptcl->NewNeighbors, total_update_list[i].newneighbors, sizeof(int) * ptcl->NewNumberOfNeighbor);
        std::memcpy(ptcl->NewPosition, total_update_list[i].newposition, sizeof(double) * 3);
        std::memcpy(ptcl->NewVelocity, total_update_list[i].newvelocity, sizeof(double) * 3);
        for (int dim = 0; dim < Dim; dim++) {
            for (int j = 0; j < HERMITE_ORDER; j++) {
                ptcl->a_irr[dim][j] = total_update_list[i].airr[dim][j];
                ptcl->a_tot[dim][j] = total_update_list[i].atot[dim][j];
            }
        }
        ptcl->NewCurrentBlockIrr = total_update_list[i].newcurrentblock_irr;
        ptcl->TimeLevelIrr = total_update_list[i].timelevel_irr;
        ptcl->TimeStepIrr = total_update_list[i].timestep_irr;
        ptcl->TimeBlockIrr = total_update_list[i].timeblock_irr;
        ptcl->NextBlockIrr = total_update_list[i].nextblock_irr;

        // UpdateIrr Routine Starts
        if (ptcl->NumberOfNeighbor != 0) // IAR modified
            ptcl->updateParticle();
        ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
        ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;

        // fprintf(workerout, "5. Receiving... MyRank: %d. PID: %d, NewNN: %d, NewPosition: %e %e %e, NewVelocity: %e %e %e, "
        //     "NewCurrentBlockIrr: %llu, TimeLevelIrr: %d, TimeStepIrr: %e, TimeBlockIrr: %llu, NextBlockIrr: %llu\n",
        //     MyRank, ptcl->PID, ptcl->NewNumberOfNeighbor,
        //     ptcl->NewPosition[0], ptcl->NewPosition[1], ptcl->NewPosition[2],
        //     ptcl->NewVelocity[0], ptcl->NewVelocity[1], ptcl->NewVelocity[2],
        //     ptcl->NextBlockIrr, ptcl->TimeLevelIrr, ptcl->TimeStepIrr,
        //     ptcl->TimeBlockIrr, ptcl->NextBlockIrr);
        // fflush(workerout);
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
}

void updateFBTermination(int update_count) {

    Particle* ptcl;

    UpdateFBTerm* total_update_list = new UpdateFBTerm[update_count];

    MPI_Recv(total_update_list, update_count, UpdateFBTermType, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i < update_count; ++i) {
        ptcl = &particles[total_update_list[i].pid];
        ptcl->isActive = true;

        // Clear CM particle information
        particles[ptcl->CMPtclIndex].clear();

        ptcl->CMPtclIndex = -1;
        ptcl->binary_state = total_update_list[i].binary_state;
        ptcl->CurrentBlockIrr = total_update_list[i].currentblock_irr;
        ptcl->CurrentTimeIrr = total_update_list[i].currenttime_irr;
        ptcl->CurrentBlockReg = total_update_list[i].currentblock_reg;
        ptcl->CurrentTimeReg = total_update_list[i].currenttime_reg;
        ptcl->NewCurrentBlockIrr = total_update_list[i].newcurrentblock_irr;
        ptcl->NextBlockIrr = total_update_list[i].nextblock_irr;
        ptcl->TimeLevelIrr = total_update_list[i].timelevel_irr;
        ptcl->TimeStepIrr = total_update_list[i].timestep_irr;
        ptcl->TimeBlockIrr = total_update_list[i].timeblock_irr;
        ptcl->TimeLevelReg = total_update_list[i].timelevel_reg;
        ptcl->TimeStepReg = total_update_list[i].timestep_reg;
        ptcl->TimeBlockReg = total_update_list[i].timeblock_reg;
        ptcl->RadiusOfNeighbor = total_update_list[i].radiusofneighbor;
        for (int dim = 0; dim < Dim; dim++) {
            for (int j = 0; j < HERMITE_ORDER; j++) {
                ptcl->a_irr[dim][j] = total_update_list[i].airr[dim][j];
                ptcl->a_reg[dim][j] = total_update_list[i].areg[dim][j];
                ptcl->a_tot[dim][j] = ptcl->a_irr[dim][j] + ptcl->a_reg[dim][j];
            }
        }
        ptcl->NumberOfNeighbor = total_update_list[i].numberofneighbors;
        std::memcpy(ptcl->Neighbors, total_update_list[i].neighbors, sizeof(int) * ptcl->NumberOfNeighbor);
        std::memcpy(ptcl->Position, total_update_list[i].position, sizeof(double) * 3);
        std::memcpy(ptcl->Velocity, total_update_list[i].velocity, sizeof(double) * 3);
    }
    delete[] total_update_list;
    total_update_list = nullptr;
}

void updateNewGroup(int ptcl_id) {

    int nodenum;
    MPI_Recv(&nodenum, 1, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    Particle* ptcl = &particles[ptcl_id];
    UpdateNewCM* update_new_cm = new UpdateNewCM[1];

    if (update_rank == nodenum) {
        assert(ptcl->isActive);
        update_new_cm[0].PID = ptcl->PID;
        update_new_cm[0].numberofmembers = ptcl->NumberOfMember;
        std::memcpy(update_new_cm[0].members, ptcl->Members, sizeof(int) * ptcl->NumberOfMember);
        std::memcpy(update_new_cm[0].position, ptcl->Position, sizeof(double) * Dim);
        std::memcpy(update_new_cm[0].velocity, ptcl->Velocity, sizeof(double) * Dim);
        update_new_cm[0].mass = ptcl->Mass;
        update_new_cm[0].radiusofneighbor = ptcl->RadiusOfNeighbor;
        update_new_cm[0].currenttime_irr = ptcl->CurrentTimeIrr;
        update_new_cm[0].currenttime_reg = ptcl->CurrentTimeReg;
        update_new_cm[0].currentblock_irr = ptcl->CurrentBlockIrr;
        update_new_cm[0].currentblock_reg = ptcl->CurrentBlockReg;
        update_new_cm[0].newcurrentblock_irr = ptcl->NewCurrentBlockIrr;
        update_new_cm[0].nextblock_irr = ptcl->NextBlockIrr;
        update_new_cm[0].timestep_irr = ptcl->TimeStepIrr;
        update_new_cm[0].timeblock_irr = ptcl->TimeBlockIrr;
        update_new_cm[0].timelevel_irr = ptcl->TimeLevelIrr;
        update_new_cm[0].timestep_reg = ptcl->TimeStepReg;
        update_new_cm[0].timeblock_reg = ptcl->TimeBlockReg;
        update_new_cm[0].timelevel_reg = ptcl->TimeLevelReg;
        update_new_cm[0].numberofneighbors = ptcl->NumberOfNeighbor;
        std::memcpy(update_new_cm[0].neighbors, ptcl->Neighbors, sizeof(int) * ptcl->NumberOfNeighbor);
        for (int i = 0; i < HERMITE_ORDER; i++) {
            for (int dim = 0; dim < Dim; dim++) {
                update_new_cm[0].airr[dim][i] = ptcl->a_irr[dim][i];
                update_new_cm[0].areg[dim][i] = ptcl->a_reg[dim][i];
            }
        }
        MPI_Bcast(update_new_cm, 1, UpdateNewCMType, nodenum, update_comm);
    }
    else {
        MPI_Bcast(update_new_cm, 1, UpdateNewCMType, nodenum, update_comm);
        ptcl->PID = update_new_cm[0].PID;
        ptcl->ParticleIndex = ptcl_id;
        ptcl->isActive = true;
        ptcl->isCMptcl = true;
        ptcl->setBinaryInterruptState(BinaryInterruptState::none);
        ptcl->NumberOfMember = update_new_cm[0].numberofmembers;
        std::memcpy(ptcl->Members, update_new_cm[0].members, sizeof(int) * ptcl->NumberOfMember);
        for (int i = 0; i < ptcl->NumberOfMember; i++) {
            particles[ptcl->Members[i]].CMPtclIndex = ptcl_id;
            particles[ptcl->Members[i]].isActive = false;
        }
        std::memcpy(ptcl->Position, update_new_cm[0].position, sizeof(double) * Dim);
        std::memcpy(ptcl->Velocity, update_new_cm[0].velocity, sizeof(double) * Dim);
        ptcl->Mass = update_new_cm[0].mass;
        ptcl->RadiusOfNeighbor = update_new_cm[0].radiusofneighbor;
        ptcl->CurrentTimeIrr = update_new_cm[0].currenttime_irr;
        ptcl->CurrentTimeReg = update_new_cm[0].currenttime_reg;
        ptcl->CurrentBlockIrr = update_new_cm[0].currentblock_irr;
        ptcl->CurrentBlockReg = update_new_cm[0].currentblock_reg;
        ptcl->NewCurrentBlockIrr = update_new_cm[0].newcurrentblock_irr;
        ptcl->NextBlockIrr = update_new_cm[0].nextblock_irr;
        ptcl->TimeStepIrr = update_new_cm[0].timestep_irr;
        ptcl->TimeBlockIrr = update_new_cm[0].timeblock_irr;
        ptcl->TimeLevelIrr = update_new_cm[0].timelevel_irr;
        ptcl->TimeStepReg = update_new_cm[0].timestep_reg;
        ptcl->TimeBlockReg = update_new_cm[0].timeblock_reg;
        ptcl->TimeLevelReg = update_new_cm[0].timelevel_reg;
        ptcl->NumberOfNeighbor = update_new_cm[0].numberofneighbors;
        std::memcpy(ptcl->Neighbors, update_new_cm[0].neighbors, sizeof(int) * ptcl->NumberOfNeighbor);
        for (int i = 0; i < HERMITE_ORDER; i++) {
            for (int dim = 0; dim < Dim; dim++) {
                ptcl->a_irr[dim][i] = update_new_cm[0].airr[dim][i];
                ptcl->a_reg[dim][i] = update_new_cm[0].areg[dim][i];
                ptcl->a_tot[dim][i] = ptcl->a_irr[dim][i] + ptcl->a_reg[dim][i];
            }
        }
    }

    delete[] update_new_cm;
    update_new_cm = nullptr;
}

void updateAfterRegCuda(int update_count, int* update_count_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(update_count);
    MPI_Recv(update_pid_list.data(), update_count, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&update_count, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);
    // fprintf(workerout, "6. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

    UpdateRegCuda* update_list = new UpdateRegCuda[update_count_list[update_rank]];

    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_list[i].pid];
        update_list[i].newnumberofneighbors = ptcl->NewNumberOfNeighbor;
        std::memcpy(update_list[i].newneighbors, ptcl->NewNeighbors, sizeof(int) * ptcl->NewNumberOfNeighbor);
        std::memcpy(update_list[i].newposition, ptcl->NewPosition, sizeof(double) * 3);
        std::memcpy(update_list[i].newvelocity, ptcl->NewVelocity, sizeof(double) * 3);
        for (int dim = 0; dim < Dim; dim++) {
            for (int j = 0; j < 2; j++) {
                update_list[i].airr[dim][j] = ptcl->a_irr[dim][j];
                update_list[i].areg[dim][j] = ptcl->a_reg[dim][j];
                update_list[i].atot[dim][j] = ptcl->a_tot[dim][j];
            }
            for (int j = 2; j < HERMITE_ORDER; j++) {
                update_list[i].areg[dim][j] = ptcl->a_reg[dim][j];
                update_list[i].atot[dim][j] = ptcl->a_tot[dim][j];
            }
        }
        // fprintf(workerout, "6. Gathering... MyRank: %d. PID: %d, NewNN: %d, NewPosition: %e %e %e, NewVelocity: %e %e %e\n",
        //     MyRank, ptcl->PID, ptcl->NewNumberOfNeighbor,
        //     ptcl->NewPosition[0], ptcl->NewPosition[1], ptcl->NewPosition[2],
        //     ptcl->NewVelocity[0], ptcl->NewVelocity[1], ptcl->NewVelocity[2]);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    UpdateRegCuda* total_update_list = new UpdateRegCuda[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateRegCudaType, total_update_list, update_count_list, displs, UpdateRegCudaType, update_comm);
    // fprintf(workerout, "6. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank]) {

            ptcl = &particles[total_update_list[i].pid];

            ptcl->NumberOfNeighbor = total_update_list[i].newnumberofneighbors;
            std::memcpy(ptcl->Neighbors, total_update_list[i].newneighbors, sizeof(int) * ptcl->NumberOfNeighbor);
            std::memcpy(ptcl->Position, total_update_list[i].newposition, sizeof(double) * 3);
            std::memcpy(ptcl->Velocity, total_update_list[i].newvelocity, sizeof(double) * 3);
            continue;
        }

        ptcl = &particles[total_update_list[i].pid];
        ptcl->NumberOfNeighbor = total_update_list[i].newnumberofneighbors;
        std::memcpy(ptcl->Neighbors, total_update_list[i].newneighbors, sizeof(int) * ptcl->NumberOfNeighbor);
        std::memcpy(ptcl->Position, total_update_list[i].newposition, sizeof(double) * 3);
        std::memcpy(ptcl->Velocity, total_update_list[i].newvelocity, sizeof(double) * 3);
        for (int dim = 0; dim < Dim; dim++) {
            for (int j = 0; j < 2; j++) {
                ptcl->a_irr[dim][j] = total_update_list[i].airr[dim][j];
                ptcl->a_reg[dim][j] = total_update_list[i].areg[dim][j];
                ptcl->a_tot[dim][j] = total_update_list[i].atot[dim][j];
            }
            for (int j = 2; j < HERMITE_ORDER; j++) {
                ptcl->a_reg[dim][j] = total_update_list[i].areg[dim][j];
                ptcl->a_tot[dim][j] = total_update_list[i].atot[dim][j];
            }
        }

        // Binary members pos&vel adjustment after regular routine // error fixed by EW 2025.5.29
        if (ptcl->isCMptcl) {
            for (int k = 0; k < ptcl->NumberOfMember; k++) {
                for (int dim=0; dim<Dim; dim++) {
                    particles[ptcl->Members[k]].Position[dim] += ptcl->NewPosition[dim] - ptcl->Position[dim];
                    particles[ptcl->Members[k]].Velocity[dim] += ptcl->NewVelocity[dim] - ptcl->Velocity[dim];
                }
            }
        }
        // fprintf(workerout, "6. Receiving... MyRank: %d. PID: %d, NewNN: %d, NewPosition: %e %e %e, NewVelocity: %e %e %e\n",
        //     MyRank, ptcl->PID, ptcl->NewNumberOfNeighbor,
        //     ptcl->NewPosition[0], ptcl->NewPosition[1], ptcl->NewPosition[2],
        //     ptcl->NewVelocity[0], ptcl->NewVelocity[1], ptcl->NewVelocity[2]);
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
}

void updateAfterRegCudaUpdate(int update_count, int* update_count_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(update_count);
    MPI_Recv(update_pid_list.data(), update_count, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&update_count, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);
    // fprintf(workerout, "7. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

    UpdateRegCudaUpdate* update_list = new UpdateRegCudaUpdate[update_count_list[update_rank]];

    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_list[i].pid];
        update_list[i].currentblock_reg = ptcl->CurrentBlockReg;
        update_list[i].currenttime_reg = ptcl->CurrentTimeReg;
        update_list[i].timelevel_reg = ptcl->TimeLevelReg;
        update_list[i].timestep_reg = ptcl->TimeStepReg;
        update_list[i].timeblock_reg = ptcl->TimeBlockReg;
        update_list[i].timelevel_irr = ptcl->TimeLevelIrr;
        update_list[i].timestep_irr = ptcl->TimeStepIrr;
        update_list[i].timeblock_irr = ptcl->TimeBlockIrr;
        update_list[i].radiusofneighbor = ptcl->RadiusOfNeighbor;
        update_list[i].nextblock_irr = ptcl->NextBlockIrr;
        // fprintf(stderr, "7. Gathering... MyRank: %d. PID: %d, CurrentBlockReg: %llu, CurrentTimeReg: %e, "
        //     "TimeLevelReg: %d, TimeStepReg: %e, TimeBlockReg: %llu, TimeLevelIrr: %d, "
        //     "TimeStepIrr: %e, TimeBlockIrr: %llu, RadiusOfNeighbor: %e, NextBlockIrr: %llu\n",
        //     MyRank, ptcl->PID,
        //     ptcl->CurrentBlockReg, ptcl->CurrentTimeReg,
        //     ptcl->TimeLevelReg, ptcl->TimeStepReg, ptcl->TimeBlockReg,
        //     ptcl->TimeLevelIrr, ptcl->TimeStepIrr, ptcl->TimeBlockIrr,
        //     ptcl->RadiusOfNeighbor, ptcl->NextBlockIrr);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    UpdateRegCudaUpdate* total_update_list = new UpdateRegCudaUpdate[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateRegCudaUpdateType,
                   total_update_list, update_count_list, displs,
                   UpdateRegCudaUpdateType, update_comm);
    // fprintf(workerout, "7. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank])
            continue;

        ptcl = &particles[total_update_list[i].pid];
        ptcl->CurrentBlockReg = total_update_list[i].currentblock_reg;
        ptcl->CurrentTimeReg = total_update_list[i].currenttime_reg;
        ptcl->TimeLevelReg = total_update_list[i].timelevel_reg;
        ptcl->TimeStepReg = total_update_list[i].timestep_reg;
        ptcl->TimeBlockReg = total_update_list[i].timeblock_reg;
        ptcl->TimeLevelIrr = total_update_list[i].timelevel_irr;
        ptcl->TimeStepIrr = total_update_list[i].timestep_irr;
        ptcl->TimeBlockIrr = total_update_list[i].timeblock_irr;
        ptcl->RadiusOfNeighbor = total_update_list[i].radiusofneighbor;
        ptcl->NextBlockIrr = total_update_list[i].nextblock_irr;

        if (ptcl->NumberOfNeighbor == 0) {
            ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
            ptcl->CurrentTimeIrr = ptcl->CurrentTimeReg;
        }
        // fprintf(stderr, "7. Receiving... MyRank: %d. PID: %d, CurrentBlockReg: %llu, CurrentTimeReg: %e, "
        //     "TimeLevelReg: %d, TimeStepReg: %e, TimeBlockReg: %llu, TimeLevelIrr: %d, "
        //     "TimeStepIrr: %e, TimeBlockIrr: %llu, RadiusOfNeighbor: %e, NextBlockIrr: %llu\n",
        //     MyRank, ptcl->PID,
        //     ptcl->CurrentBlockReg, ptcl->CurrentTimeReg,
        //     ptcl->TimeLevelReg, ptcl->TimeStepReg, ptcl->TimeBlockReg,
        //     ptcl->TimeLevelIrr, ptcl->TimeStepIrr, ptcl->TimeBlockIrr,
        //     ptcl->RadiusOfNeighbor, ptcl->NextBlockIrr);
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
}
#endif