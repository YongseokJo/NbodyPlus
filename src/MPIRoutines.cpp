#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "global.h"
#include <mpi.h>
#include <cstddef>
#include "Queue.h"
#include "cuda/cuda_defs.h"


template <typename T>
void InitialAssignmentOfTasks(T data, int NumTask);


MPI_Datatype createQueueType();
MPI_Datatype createIparticleType();
MPI_Datatype createJparticleType();

void initializeMPI(int argc, char *argv[]) {
	/* MPI Initialization */
	//MPI_Win win;
	//int MyRank;
	//int NumberOfProcessor;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
	MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcessor);
	NumberOfWorker = NumberOfProcessor - 1;

	if (NumberOfProcessor < 2) {
		std::cerr << "This program requires at least 2 processes.\n";
		MPI_Finalize();
	}

	QueueType		= createQueueType();
	IparticleType	= createIparticleType();
	JparticleType	= createJparticleType();
	/*
	// comm for each node
	MPI_Comm shmcomm;
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
	int local_rank, local_size;
	MPI_Comm_rank(shmcomm, &local_rank);
	MPI_Comm_size(shmcomm, &local_size);
	*/


	// Create a shared memory communicator
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, MyRank, MPI_INFO_NULL, &shared_comm);

	// Get rank and size in the shared communicator
	int shared_rank, shared_size;
	MPI_Comm_rank(shared_comm, &shared_rank);
	MPI_Comm_size(shared_comm, &shared_size);
	if (MyRank == ROOT)
		fprintf(stdout, "MyRank = %d, NumberOfProcessor = %d : Shared Rank = %d, Shared size = %d\n", MyRank, NumberOfProcessor, shared_rank, shared_size);

	// Allocate shared memory
	if (shared_rank == 0) {
		//MPI_Win_allocate_shared(sizeof(int), sizeof(int), MPI_INFO_NULL, shared_comm, &shared_mem, &win);
		MPI_Win_allocate_shared(sizeof(Particle) * MaxNumberOfParticle, sizeof(Particle), 
		MPI_INFO_NULL, shared_comm, &particles_original, &win);
		MPI_Win_allocate_shared(sizeof(GlobalVariable), sizeof(GlobalVariable),
		 MPI_INFO_NULL, shared_comm, &global_variable_original, &win2);
	} else {
		MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles_original, &win);
		MPI_Win_allocate_shared(0, sizeof(GlobalVariable), MPI_INFO_NULL, shared_comm, &global_variable_original, &win2);
	}
	// Query shared memory of rank 0
	
	MPI_Aint size_bytes;
	int disp_unit;

	MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
	MPI_Win_shared_query(win2, 0, &size_bytes, &disp_unit, &global_variable);
}

void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i],   1, MPI_INT,    i+1, TAG, MPI_COMM_WORLD);
	}
}


void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		//MPI_Isend(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);

		MPI_Send(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD);
		MPI_Send(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD);
	}
}


void InitialAssignmentOfTasks(int* data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD);
	}
}

void InitialAssignmentOfTasks(int data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		//std::cerr << "InitialAssignmentOfTasks out of" << NumTask<< ": " << i << std::endl;
		if (i >= NumTask)
			break;
		//MPI_Isend(&data, 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data, 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD);
	}
	//fprintf(stderr, "Number of tasks assigned = %d\n", i);
	//fflush(stderr);
}

void InitialAssignmentOfTasks(Queue queue, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		MPI_Send(&queue, 1, QueueType, i+1, TAG, MPI_COMM_WORLD);
	}
}



void broadcastFromRoot(int &data) {
	MPI_Bcast(&data, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
}
void broadcastFromRoot(double &data) {
	MPI_Bcast(&data, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
}
void broadcastFromRoot(ULL &data) {
	MPI_Bcast(&data, 1, MPI_UNSIGNED_LONG_LONG, ROOT, MPI_COMM_WORLD);
}

void ParticleSynchronization() {
		//std::cout << "Particle synchronization." << std::endl;
		MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
		MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
		int task=100;
		int completed_tasks = 0;
		//int ptcl_id=104;
		InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
		//std::cerr << "before, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		MPI_Win_sync(win);  // Synchronize memory
		MPI_Barrier(shared_comm);
		//std::cerr << "after, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		for (int i=0; i<NumberOfWorker; i++) {

			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &requests[i]);
		}
		MPI_Waitall(NumberOfWorker, requests, statuses);
		//std::cerr << "Done" << std::endl;
}

void FurtherAssignmentOfTasks() {

}

MPI_Datatype createQueueType() {
    MPI_Datatype QueueType;
    int blocklen[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT8_T, MPI_INT, MPI_DOUBLE};
	MPI_Aint disp[3], base;

	Queue sample;
    MPI_Get_address(&sample,           &base);
    MPI_Get_address(&sample.task,      &disp[0]);
    MPI_Get_address(&sample.pid,       &disp[1]);
    MPI_Get_address(&sample.next_time, &disp[2]);

	for (int i = 0; i < 3; ++i) disp[i] -= base;

    MPI_Type_create_struct(3, blocklen, disp, types, &QueueType);
    MPI_Type_commit(&QueueType);

    return QueueType;
}

MPI_Datatype createIparticleType() {
	MPI_Datatype IparticleType;
	int blocklen[8] = {1,1,1,1,1,1,1,1};
#ifdef CUDA_FLOAT
	MPI_Datatype types[8] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
							MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
#else
	MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
							MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
#endif
	MPI_Aint disp[8], base;

	Iparticle sample;
	MPI_Get_address(&sample,		&base);
	MPI_Get_address(&sample.posx,	&disp[0]);
	MPI_Get_address(&sample.posy,	&disp[1]);
	MPI_Get_address(&sample.posz,	&disp[2]);
	MPI_Get_address(&sample.r2,		&disp[3]);
	MPI_Get_address(&sample.velx,	&disp[4]);
	MPI_Get_address(&sample.vely,	&disp[5]);
	MPI_Get_address(&sample.velz,	&disp[6]);
	MPI_Get_address(&sample.dtr,	&disp[7]);

	for (int i = 0; i < 8; ++i) disp[i] -= base;

	MPI_Type_create_struct(8, blocklen, disp, types, &IparticleType);
	MPI_Type_commit(&IparticleType);

	return IparticleType;
}

MPI_Datatype createJparticleType() {
	MPI_Datatype JparticleType;
	int blocklen[8] = {1,1,1,1,1,1,1,1};
#ifdef CUDA_FLOAT
	MPI_Datatype types[8] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
							MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
#else
	MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
							MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG_LONG};
#endif
	MPI_Aint disp[8], base;

	Jparticle sample;
	MPI_Get_address(&sample,		&base);
	MPI_Get_address(&sample.posx,	&disp[0]);
	MPI_Get_address(&sample.posy,	&disp[1]);
	MPI_Get_address(&sample.posz,	&disp[2]);
	MPI_Get_address(&sample.mass,	&disp[3]);
	MPI_Get_address(&sample.velx,	&disp[4]);
	MPI_Get_address(&sample.vely,	&disp[5]);
	MPI_Get_address(&sample.velz,	&disp[6]);
	MPI_Get_address(&sample.index,	&disp[7]);

	for (int i = 0; i < 8; ++i) disp[i] -= base;

	MPI_Type_create_struct(8, blocklen, disp, types, &JparticleType);
	MPI_Type_commit(&JparticleType);

	return JparticleType;
}