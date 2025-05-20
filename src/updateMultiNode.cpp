#ifdef MultiNode

#include "global.h"
#include "Queue.h"

void updateInitAcc01(int ptcl_id, int* update_count_list, UpdateInitAcc* update_list, UpdateInitAcc* total_update_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(ptcl_id);
    MPI_Recv(update_pid_list.data(), ptcl_id, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&ptcl_id, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);

    update_list = new UpdateInitAcc[update_count_list[update_rank]];
    fprintf(stderr, "1. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);
    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_list[i].pid];
        update_list[i].acc1[0] = ptcl->a_tot[0][0];
        update_list[i].acc1[1] = ptcl->a_tot[1][0];
        update_list[i].acc1[2] = ptcl->a_tot[2][0];
        update_list[i].acc2[0] = ptcl->a_tot[0][1];
        update_list[i].acc2[1] = ptcl->a_tot[1][1];
        update_list[i].acc2[2] = ptcl->a_tot[2][1];
        fprintf(stderr, "1. Gathering... MyRank: %d. PID: %d, acc1: %e %e %e, acc2: %e %e %e\n",
            MyRank, ptcl->PID, ptcl->a_tot[0][0], ptcl->a_tot[1][0], ptcl->a_tot[2][0],
            ptcl->a_tot[0][1], ptcl->a_tot[1][1], ptcl->a_tot[2][1]);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    total_update_list = new UpdateInitAcc[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateInitAccType, total_update_list, update_count_list, displs, UpdateInitAccType, update_comm);
    fprintf(stderr, "1. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);

    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank])
            continue;

        ptcl = &particles[total_update_list[i].pid];
        ptcl->a_tot[0][0] = total_update_list[i].acc1[0];
        ptcl->a_tot[1][0] = total_update_list[i].acc1[1];
        ptcl->a_tot[2][0] = total_update_list[i].acc1[2];
        ptcl->a_tot[0][1] = total_update_list[i].acc2[0];
        ptcl->a_tot[1][1] = total_update_list[i].acc2[1];
        ptcl->a_tot[2][1] = total_update_list[i].acc2[2];
        fprintf(stderr, "1. Receiving... MyRank: %d. PID: %d, acc1: %e %e %e, acc2: %e %e %e\n",
            MyRank, ptcl->PID, ptcl->a_tot[0][0], ptcl->a_tot[1][0], ptcl->a_tot[2][0],
            ptcl->a_tot[0][1], ptcl->a_tot[1][1], ptcl->a_tot[2][1]);
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
}

void updateInitAcc23(int ptcl_id, int* update_count_list, UpdateInitAcc* update_list, UpdateInitAcc* total_update_list, int* displs) {

    std::vector<int> update_pid_list;
    Particle* ptcl;

    update_pid_list.resize(ptcl_id);
    MPI_Recv(update_pid_list.data(), ptcl_id, MPI_INT, ROOT, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Allgather(&ptcl_id, 1, MPI_INT, update_count_list, 1, MPI_INT, update_comm);

    update_list = new UpdateInitAcc[update_count_list[update_rank]];
    fprintf(stderr, "2. After Allgather... MyRank: %d. list_size: %d\n", MyRank, update_count_list[update_rank]);
    for (int i=0; i<update_count_list[update_rank]; i++) {
        update_list[i].pid = update_pid_list[i];
        ptcl = &particles[update_list[i].pid];
        update_list[i].acc1[0] = ptcl->a_tot[0][2];
        update_list[i].acc1[1] = ptcl->a_tot[1][2];
        update_list[i].acc1[2] = ptcl->a_tot[2][2];
        update_list[i].acc2[0] = ptcl->a_tot[0][3];
        update_list[i].acc2[1] = ptcl->a_tot[1][3];
        update_list[i].acc2[2] = ptcl->a_tot[2][3];
        fprintf(stderr, "2. Gathering... MyRank: %d. PID: %d, acc1: %e %e %e, acc2: %e %e %e\n",
            MyRank, ptcl->PID, ptcl->a_tot[0][2], ptcl->a_tot[1][2], ptcl->a_tot[2][2],
            ptcl->a_tot[0][3], ptcl->a_tot[1][3], ptcl->a_tot[2][3]);
    }

    int total_recv_count = 0;
    for (int i = 0; i < NumberOfNode; ++i) {
        displs[i] = total_recv_count;
        total_recv_count += update_count_list[i];
    }
    total_update_list = new UpdateInitAcc[total_recv_count];
    MPI_Allgatherv(update_list, update_count_list[update_rank], UpdateInitAccType, total_update_list, update_count_list, displs, UpdateInitAccType, update_comm);
    fprintf(stderr, "2. After Allgatherv... MyRank: %d. total_recv_count: %d\n", MyRank, total_recv_count);
    for (int i=0; i<total_recv_count; i++) {
        if (i >= displs[update_rank] && i < displs[update_rank] + update_count_list[update_rank])
            continue;

        ptcl = &particles[total_update_list[i].pid];
        ptcl->a_tot[0][2] = total_update_list[i].acc1[0];
        ptcl->a_tot[1][2] = total_update_list[i].acc1[1];
        ptcl->a_tot[2][2] = total_update_list[i].acc1[2];
        ptcl->a_tot[0][3] = total_update_list[i].acc2[0];
        ptcl->a_tot[1][3] = total_update_list[i].acc2[1];
        ptcl->a_tot[2][3] = total_update_list[i].acc2[2];
        fprintf(stderr, "2. Receiving... MyRank: %d. PID: %d, acc1: %e %e %e, acc2: %e %e %e\n",
            MyRank, ptcl->PID, ptcl->a_tot[0][2], ptcl->a_tot[1][2], ptcl->a_tot[2][2],
            ptcl->a_tot[0][3], ptcl->a_tot[1][3], ptcl->a_tot[2][3]);
    }
    delete [] update_list;
    update_list = nullptr;
    delete [] total_update_list;
    total_update_list = nullptr;
}

#endif