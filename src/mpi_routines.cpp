#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "global_state.h"
#include "global.h"
#include <mpi.h>
#include <cstddef>
#include "queue.h"
#include "cuda/cuda_defs.h"


template <typename T>
void InitialAssignmentOfTasks(T data, int NumTask);


MPI_Datatype createQueueType();
MPI_Datatype createi_particle_tType();
MPI_Datatype createj_particle_tType();

void initializeMPI(int argc, char *argv[]) {

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
	num_workers = num_processors - 1;

	if (num_processors < 2) {
		std::cerr << "This program requires at least 2 processes.\n";
		MPI_Finalize();
	}

	queue_type_mpi		= createQueueType();
	iparticle_type_mpi	= createi_particle_tType();
	jparticle_type_mpi	= createj_particle_tType();

	// Create a shared memory communicator
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, my_rank, MPI_INFO_NULL, &shared_comm);

	// Get rank and size in the shared communicator
	int shared_rank, shared_size;
	MPI_Comm_rank(shared_comm, &shared_rank);
	MPI_Comm_size(shared_comm, &shared_size);

	if (my_rank == ROOT)
		fprintf(stdout, "my_rank = %d, num_processors = %d : Shared Rank = %d, Shared size = %d\n", my_rank, num_processors, shared_rank, shared_size);

	// Allocate shared memory
	if (shared_rank == 0) {
		MPI_Win_allocate_shared(sizeof(Particle) * MAX_NUM_PARTICLE, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles, &win);
		MPI_Win_allocate_shared(sizeof(GlobalVariable), sizeof(GlobalVariable), MPI_INFO_NULL, shared_comm, &g_state, &win2);
		MPI_Win_allocate_shared(sizeof(int) * MAX_NUM_PARTICLE * MAX_NUM_NEIGHBOR, sizeof(int), MPI_INFO_NULL, shared_comm, &neighbors, &win3);
		MPI_Win_allocate_shared(sizeof(int) * MAX_NUM_PARTICLE * MAX_NUM_NEIGHBOR, sizeof(int), MPI_INFO_NULL, shared_comm, &new_neighbors, &win4);
	} else {
		MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles, &win);
		MPI_Win_allocate_shared(0, sizeof(GlobalVariable), MPI_INFO_NULL, shared_comm, &g_state, &win2);
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, shared_comm, &neighbors, &win3);
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, shared_comm, &new_neighbors, &win4);
	}
	
  // Query shared memory of rank 0
	MPI_Aint size_bytes;
	int disp_unit;

	MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
	MPI_Win_shared_query(win2, 0, &size_bytes, &disp_unit, &g_state);
	MPI_Win_shared_query(win3, 0, &size_bytes, &disp_unit, &neighbors);
	MPI_Win_shared_query(win4, 0, &size_bytes, &disp_unit, &new_neighbors);

	// Allocate SoA particle data via MPI shared memory
	particle_data.allocate_shared(MAX_NUM_PARTICLE, shared_comm);
}

void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG) {
	for (int i=0; i<num_workers; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i],   1, MPI_INT,    i+1, TAG, MPI_COMM_WORLD);
	}
}


void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG) {
	for (int i=0; i<num_workers; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		//MPI_Isend(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);

		MPI_Send(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, MPI_COMM_WORLD);
		MPI_Send(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, MPI_COMM_WORLD);
	}
}


void InitialAssignmentOfTasks(int* data, int NumTask, int TAG) {
	for (int i=0; i<num_workers; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i], 1, MPI_INT, i+1, TAG, MPI_COMM_WORLD);
	}
}

void InitialAssignmentOfTasks(int data, int NumTask, int TAG) {
	for (int i=0; i<num_workers; i++) {
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
	for (int i=0; i<num_workers; i++) {
		if (i >= NumTask) break;
		MPI_Send(&queue, 1, queue_type_mpi, i+1, TAG, MPI_COMM_WORLD);
	}
}



void broadcastFromRoot(int &data) {
	MPI_Bcast(&data, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
}
void broadcastFromRoot(double &data) {
	MPI_Bcast(&data, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
}
void broadcastFromRoot(ull_t &data) {
	MPI_Bcast(&data, 1, MPI_UNSIGNED_LONG_LONG, ROOT, MPI_COMM_WORLD);
}

void ParticleSynchronization() {
		//std::cout << "Particle synchronization." << std::endl;
		MPI_Request requests[num_processors];  // Pointer to the request handle
		MPI_Status statuses[num_processors];    // Pointer to the status object
		int task=100;
		int completed_tasks = 0;
		//int ptcl_id=104;
		InitialAssignmentOfTasks(task, num_workers, TASK_TAG);
		//std::cerr << "before, Rank=" << my_rank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].current_time_irr << std::endl;
		MPI_Win_sync(win);  // TASK_SYNCHRONIZE memory
		// Sync SoA particle data
		particle_data.sync_all();
		MPI_Barrier(shared_comm);
		//std::cerr << "after, Rank=" << my_rank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].current_time_irr << std::endl;
		for (int i=0; i<num_workers; i++) {

			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &requests[i]);
		}
		MPI_Waitall(num_workers, requests, statuses);
		//std::cerr << "Done" << std::endl;
}

void FurtherAssignmentOfTasks() {

}

MPI_Datatype createQueueType() {
    MPI_Datatype queue_type_mpi;
    int blocklen[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT8_T, MPI_INT, MPI_DOUBLE};
	MPI_Aint disp[3], base;

	Queue sample;
    MPI_Get_address(&sample,           &base);
    MPI_Get_address(&sample.task,      &disp[0]);
    MPI_Get_address(&sample.pid,       &disp[1]);
    MPI_Get_address(&sample.next_time, &disp[2]);

	for (int i = 0; i < 3; ++i) disp[i] -= base;

    MPI_Type_create_struct(3, blocklen, disp, types, &queue_type_mpi);
    MPI_Type_commit(&queue_type_mpi);

    return queue_type_mpi;
}

MPI_Datatype createi_particle_tType() {
	MPI_Datatype iparticle_type_mpi;
	int blocklen[8] = {1,1,1,1,1,1,1,1};
#ifdef CUDA_FLOAT
	MPI_Datatype types[8] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
							MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
#else
	MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
							MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
#endif
	MPI_Aint disp[8], base;

	i_particle_t sample;
	MPI_Get_address(&sample,		&base);
	MPI_Get_address(&sample.pos_x,	&disp[0]);
	MPI_Get_address(&sample.pos_y,	&disp[1]);
	MPI_Get_address(&sample.pos_z,	&disp[2]);
	MPI_Get_address(&sample.radius_sq,		&disp[3]);
	MPI_Get_address(&sample.vel_x,	&disp[4]);
	MPI_Get_address(&sample.vel_y,	&disp[5]);
	MPI_Get_address(&sample.vel_z,	&disp[6]);
	MPI_Get_address(&sample.dt_reg,	&disp[7]);

	for (int i = 0; i < 8; ++i) disp[i] -= base;

	MPI_Type_create_struct(8, blocklen, disp, types, &iparticle_type_mpi);
	MPI_Type_commit(&iparticle_type_mpi);

	return iparticle_type_mpi;
}

MPI_Datatype createj_particle_tType() {
	MPI_Datatype jparticle_type_mpi;
	int blocklen[8] = {1,1,1,1,1,1,1,1};
#ifdef CUDA_FLOAT
	MPI_Datatype types[8] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
							MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
#else
	MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
							MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG_LONG};
#endif
	MPI_Aint disp[8], base;

	j_particle_t sample;
	MPI_Get_address(&sample,		&base);
	MPI_Get_address(&sample.pos_x,	&disp[0]);
	MPI_Get_address(&sample.pos_y,	&disp[1]);
	MPI_Get_address(&sample.pos_z,	&disp[2]);
	MPI_Get_address(&sample.mass,	&disp[3]);
	MPI_Get_address(&sample.vel_x,	&disp[4]);
	MPI_Get_address(&sample.vel_y,	&disp[5]);
	MPI_Get_address(&sample.vel_z,	&disp[6]);
	MPI_Get_address(&sample.index,	&disp[7]);

	for (int i = 0; i < 8; ++i) disp[i] -= base;

	MPI_Type_create_struct(8, blocklen, disp, types, &jparticle_type_mpi);
	MPI_Type_commit(&jparticle_type_mpi);

	return jparticle_type_mpi;
}