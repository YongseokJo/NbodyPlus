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
    fprintf(workerout, "1. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

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
        fprintf(workerout, "1. Gathering... MyRank: %d. PID: %d, NN: %d. airr0: %e %e %e, airr1: %e %e %e, areg0: %e %e %e, areg1: %e %e %e\n",
            MyRank, ptcl->PID, ptcl->NumberOfNeighbor,
            ptcl->a_irr[0][0], ptcl->a_irr[1][0], ptcl->a_irr[2][0],
            ptcl->a_irr[0][1], ptcl->a_irr[1][1], ptcl->a_irr[2][1],
            ptcl->a_reg[0][0], ptcl->a_reg[1][0], ptcl->a_reg[2][0],
            ptcl->a_reg[0][1], ptcl->a_reg[1][1], ptcl->a_reg[2][1]);
        fprintf(workerout, "(");
        for (int j=0; j<ptcl->NumberOfNeighbor; j++)
            fprintf(workerout, "%d ", ptcl->Neighbors[j]);
        fprintf(workerout, ")\n");

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
    fprintf(workerout, "1. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);

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

        fprintf(workerout, "1. Receiving... MyRank: %d. PID: %d, NN: %d, airr0: %e %e %e, airr1: %e %e %e, areg0: %e %e %e, areg1: %e %e %e\n",
            MyRank, ptcl->PID, ptcl->NumberOfNeighbor,
            ptcl->a_irr[0][0], ptcl->a_irr[1][0], ptcl->a_irr[2][0],
            ptcl->a_irr[0][1], ptcl->a_irr[1][1], ptcl->a_irr[2][1],
            ptcl->a_reg[0][0], ptcl->a_reg[1][0], ptcl->a_reg[2][0],
            ptcl->a_reg[0][1], ptcl->a_reg[1][1], ptcl->a_reg[2][1]);
        fprintf(workerout, "(");
        for (int j=0; j<ptcl->NumberOfNeighbor; j++)
            fprintf(workerout, "%d ", ptcl->Neighbors[j]);
        fprintf(workerout, ")\n");
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
    fprintf(workerout, "2. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

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
        fprintf(workerout, "2. Gathering... MyRank: %d. PID: %d, airr2: %e %e %e, airr3: %e %e %e, areg2: %e %e %e, areg3: %e %e %e\n",
            MyRank, ptcl->PID,
            ptcl->a_irr[0][2], ptcl->a_irr[1][2], ptcl->a_irr[2][2],
            ptcl->a_irr[0][3], ptcl->a_irr[1][3], ptcl->a_irr[2][3],
            ptcl->a_reg[0][2], ptcl->a_reg[1][2], ptcl->a_reg[2][2],
            ptcl->a_reg[0][3], ptcl->a_reg[1][3], ptcl->a_reg[2][3]);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    UpdateInitAcc2* total_update_list = new UpdateInitAcc2[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateInitAcc2Type, total_update_list, update_count_list, displs, UpdateInitAcc2Type, update_comm);
    fprintf(workerout, "2. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
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
        fprintf(workerout, "2. Receiving... MyRank: %d. PID: %d, airr2: %e %e %e, airr3: %e %e %e, areg2: %e %e %e, areg3: %e %e %e\n",
            MyRank, ptcl->PID,
            ptcl->a_irr[0][2], ptcl->a_irr[1][2], ptcl->a_irr[2][2],
            ptcl->a_irr[0][3], ptcl->a_irr[1][3], ptcl->a_irr[2][3],
            ptcl->a_reg[0][2], ptcl->a_reg[1][2], ptcl->a_reg[2][2],
            ptcl->a_reg[0][3], ptcl->a_reg[1][3], ptcl->a_reg[2][3]);
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
    fprintf(workerout, "3. After Recv... MyRank: %d. list_size: %d\n", MyRank, update_count);

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
    fprintf(workerout, "4. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);

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
        fprintf(workerout, "4. Gathering... MyRank: %d. PID: %d, TimeStepReg: %e, TimeBlockReg: %llu, TimeLevelReg: %d, "
            "TimeStepIrr: %e, TimeBlockIrr: %llu, TimeLevelIrr: %d, CurrentTimeIrr: %e, CurrentTimeReg: %e, "
            "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n",
            MyRank, ptcl->PID,
            ptcl->TimeStepReg, ptcl->TimeBlockReg, ptcl->TimeLevelReg,
            ptcl->TimeStepIrr, ptcl->TimeBlockIrr, ptcl->TimeLevelIrr,
            ptcl->CurrentTimeIrr, ptcl->CurrentTimeReg, ptcl->CurrentBlockIrr, ptcl->CurrentBlockReg);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    UpdateTime* total_update_list = new UpdateTime[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateTimeType, total_update_list, update_count_list, displs, UpdateTimeType, update_comm);
    fprintf(workerout, "4. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
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
        fprintf(workerout, "4. Receiving... MyRank: %d. PID: %d, TimeStepReg: %e, TimeBlockReg: %llu, TimeLevelReg: %d, "
            "TimeStepIrr: %e, TimeBlockIrr: %llu, TimeLevelIrr: %d, CurrentTimeIrr: %e, CurrentTimeReg: %e, "
            "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n",
            MyRank, ptcl->PID,
            ptcl->TimeStepReg, ptcl->TimeBlockReg, ptcl->TimeLevelReg,
            ptcl->TimeStepIrr, ptcl->TimeBlockIrr, ptcl->TimeLevelIrr,
            ptcl->CurrentTimeIrr, ptcl->CurrentTimeReg, ptcl->CurrentBlockIrr, ptcl->CurrentBlockReg);
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
#endif