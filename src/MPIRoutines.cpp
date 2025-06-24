#include <iostream>
#include <cstddef>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "global.h"
#include <mpi.h>
#include "Queue.h"


template <typename T>
void InitialAssignmentOfTasks(T data, int NumTask);


MPI_Datatype createQueueType();
#ifdef MultiNode
MPI_Datatype createUpdateInitAcc1Type();
MPI_Datatype createUpdateInitAcc2Type();
MPI_Datatype createUpdateTimeType();
MPI_Datatype createUpdateTimeCorrType();
MPI_Datatype createUpdateIrrForceType();
MPI_Datatype createUpdateBinaryType();
MPI_Datatype createUpdateFBTermType();
MPI_Datatype createUpdateNewCMType();
MPI_Datatype createUpdateRegCudaType();
MPI_Datatype createUpdateRegCudaUpdateType();
#ifdef SEVN
MPI_Datatype createUpdateSEVN0Type();
MPI_Datatype createUpdateSEVN1Type();
#endif
#endif

void initializeMPI(int argc, char *argv[]) {
	/* MPI Initialization */
	//MPI_Win win;
	//int MyRank;
	//int NumberOfProcessor;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
	MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcessor); // (Query MultiNode) This should be checked; Do I have to explictly set NumberOfWorker & NumberOfProcessor?
	// (Query MultiNode) We need global variable NumberOfNode; How/Where to define it?
	// (Query MultiNode) NumberOfProcessor % NumberOfNode should be 0 !!!
	NumberOfWorker = NumberOfProcessor - 1;

	if (NumberOfProcessor < 2) {
		std::cerr << "This program requires at least 2 processes.\n";
		MPI_Finalize();
	}

	QueueType				= createQueueType();
#ifdef MultiNode
	UpdateInitAcc1Type		= createUpdateInitAcc1Type();
	UpdateInitAcc2Type		= createUpdateInitAcc2Type();
	UpdateTimeType 			= createUpdateTimeType();
	UpdateTimeCorrType 		= createUpdateTimeCorrType();
	UpdateIrrForceType 		= createUpdateIrrForceType();
	UpdateBinaryType 		= createUpdateBinaryType();
	UpdateFBTermType 		= createUpdateFBTermType();
	UpdateNewCMType 		= createUpdateNewCMType();
	UpdateRegCudaType 		= createUpdateRegCudaType();
	UpdateRegCudaUpdateType = createUpdateRegCudaUpdateType();
#ifdef SEVN
	UpdateSEVN0Type 		= createUpdateSEVN0Type();
	UpdateSEVN1Type 		= createUpdateSEVN1Type();
#endif
#endif
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
	int shared_size; // (Query MultiNode) shared_rank is global variable by EW 2025.5.13

	MPI_Comm_rank(shared_comm, &shared_rank);
	MPI_Comm_size(shared_comm, &shared_size);
	fprintf(stderr,"My Rank =%d : Shared Rank = %d, Shared size = %d\n", MyRank, shared_rank, shared_size);

	// Allocate shared memory
	if (shared_rank == 0) {
		//MPI_Win_allocate_shared(sizeof(int), sizeof(int), MPI_INFO_NULL, shared_comm, &shared_mem, &win);
		MPI_Win_allocate_shared(sizeof(Particle) * MaxNumberOfParticle, sizeof(Particle), 
		MPI_INFO_NULL, shared_comm, &particles_original, &win);
		MPI_Win_allocate_shared(sizeof(GlobalVariable), sizeof(GlobalVariable),
		 MPI_INFO_NULL, shared_comm, &global_variable_original, &win2);
		MPI_Win_allocate_shared(sizeof(int)*MaxNumberOfParticle, sizeof(int),
		 MPI_INFO_NULL, shared_comm, &ActiveIndexToOriginalIndex_orginal, &win3);
	} else {
		MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, shared_comm, &particles_original, &win);
		MPI_Win_allocate_shared(0, sizeof(GlobalVariable), MPI_INFO_NULL, shared_comm, &global_variable_original, &win2);
		MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, shared_comm, &ActiveIndexToOriginalIndex_orginal, &win3);
	}
	// Query shared memory of rank 0
	
	MPI_Aint size_bytes;
	int disp_unit;

	MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
	MPI_Win_shared_query(win2, 0, &size_bytes, &disp_unit, &global_variable);
	MPI_Win_shared_query(win3, 0, &size_bytes, &disp_unit, &ActiveIndexToOriginalIndex);

#ifdef MultiNode
	NumberOfNode = 1; // Initialization because it is used in WorkerRoutines.cpp
	int color = (shared_rank == 1) ? 0 : MPI_UNDEFINED;
	MPI_Comm_split(MPI_COMM_WORLD, color, MyRank, &update_comm);
	int update_size = -1; // (Query MultiNode) update_rank is global variable by EW 2025.5.13
	if (shared_rank == 1) {
		MPI_Comm_rank(update_comm, &update_rank);
		MPI_Comm_size(update_comm, &update_size);

		NumberOfNode = update_size;
		fprintf(stderr,"MyRank =%d : Update Rank = %d, Update size = %d\n", MyRank, update_rank, update_size);

		if (MyRank == 1) {
			ranks_update_comm = new int[NumberOfNode];
			MPI_Gather(&MyRank, 1, MPI_INT, ranks_update_comm, 1, MPI_INT, ROOT, update_comm);
			for (int i=0; i<NumberOfNode; i++) {
				fprintf(stderr,"MyRank: %d. Update Rank[%d]: %d\n", MyRank, i, ranks_update_comm[i]);
			}
			assert(MyRank == 1);
			assert(shared_rank == 1);
			assert(update_rank == 0);
		}
		else {
			MPI_Gather(&MyRank, 1, MPI_INT, NULL, 0, MPI_INT, ROOT, update_comm);
		}
	}

	if (MyRank == 1) {
		MPI_Send(&NumberOfNode, 1, MPI_INT, ROOT, 0, MPI_COMM_WORLD);
		MPI_Send(ranks_update_comm, NumberOfNode, MPI_INT, ROOT, 1, MPI_COMM_WORLD);
	}
	else if (MyRank == ROOT) {
		MPI_Recv(&NumberOfNode, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		ranks_update_comm = new int[NumberOfNode];
		MPI_Recv(ranks_update_comm, NumberOfNode, MPI_INT, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		fprintf(stderr, "In ROOT... NumberOfNode: %d\n", NumberOfNode);
		for (int i=0; i<NumberOfNode; i++) {
			fprintf(stderr,"Update Rank[%d]: %d\n", i, ranks_update_comm[i]);
		}
	}
#endif
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
    int block_lengths[3] = {1, 1, 1}; // Number of elements in each field
    MPI_Aint offsets[3];
    MPI_Datatype types[3] = {MPI_INT8_T, MPI_INT, MPI_DOUBLE}; // Match the types in the struct

    // Calculate offsets
    offsets[0] = offsetof(Queue, task);
    offsets[1] = offsetof(Queue, pid);
    offsets[2] = offsetof(Queue, next_time);

    // Create the struct datatype
    MPI_Type_create_struct(3, block_lengths, offsets, types, &QueueType);
    MPI_Type_commit(&QueueType);

    return QueueType;
}

#ifdef MultiNode
MPI_Datatype createUpdateInitAcc1Type() {

	int block_lengths[6] = {1, 1, 3, 3, 3, 3}; // Number of elements in each field
	MPI_Aint offsets[6];
	MPI_Datatype types[6] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateInitAcc1, pid);
	offsets[1] = offsetof(UpdateInitAcc1, numberofneighbors);
	offsets[2] = offsetof(UpdateInitAcc1, airr0);
	offsets[3] = offsetof(UpdateInitAcc1, airr1);
	offsets[4] = offsetof(UpdateInitAcc1, areg0);
	offsets[5] = offsetof(UpdateInitAcc1, areg1);

	// Create the struct datatype
	MPI_Type_create_struct(6, block_lengths, offsets, types, &UpdateInitAcc1Type);
	MPI_Type_commit(&UpdateInitAcc1Type);

	return UpdateInitAcc1Type;
}

MPI_Datatype createUpdateInitAcc2Type() {

	int block_lengths[5] = {1, 3, 3, 3, 3}; // Number of elements in each field
	MPI_Aint offsets[5];
	MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateInitAcc2, pid);
	offsets[1] = offsetof(UpdateInitAcc2, airr2);
	offsets[2] = offsetof(UpdateInitAcc2, airr3);
	offsets[3] = offsetof(UpdateInitAcc2, areg2);
	offsets[4] = offsetof(UpdateInitAcc2, areg3);

	// Create the struct datatype
	MPI_Type_create_struct(5, block_lengths, offsets, types, &UpdateInitAcc2Type);
	MPI_Type_commit(&UpdateInitAcc2Type);

	return UpdateInitAcc2Type;
}

MPI_Datatype createUpdateTimeType() {

	int block_lengths[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Number of elements in each field
	MPI_Aint offsets[11];
	MPI_Datatype types[11] = {MPI_INT, MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG, MPI_INT, 
								MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG, MPI_INT,
								MPI_DOUBLE, MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateTime, pid);
	offsets[1] = offsetof(UpdateTime, timestep_reg);
	offsets[2] = offsetof(UpdateTime, timeblock_reg);
	offsets[3] = offsetof(UpdateTime, timelevel_reg);
	offsets[4] = offsetof(UpdateTime, timestep_irr);
	offsets[5] = offsetof(UpdateTime, timeblock_irr);
	offsets[6] = offsetof(UpdateTime, timelevel_irr);
	offsets[7] = offsetof(UpdateTime, currenttime_irr);
	offsets[8] = offsetof(UpdateTime, currenttime_reg);
	offsets[9] = offsetof(UpdateTime, currentblock_irr);
	offsets[10] = offsetof(UpdateTime, currentblock_reg);

	// Create the struct datatype
	MPI_Type_create_struct(11, block_lengths, offsets, types, &UpdateTimeType);
	MPI_Type_commit(&UpdateTimeType);

	return UpdateTimeType;
}

MPI_Datatype createUpdateTimeCorrType() {

	int block_lengths[4] = {1, 1, 1, 1}; // Number of elements in each field
	MPI_Aint offsets[4];
	MPI_Datatype types[4] = {MPI_INT, MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG, MPI_INT}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateTimeCorr, pid);
	offsets[1] = offsetof(UpdateTimeCorr, timestep_irr);
	offsets[2] = offsetof(UpdateTimeCorr, timeblock_irr);
	offsets[3] = offsetof(UpdateTimeCorr, timelevel_irr);

	// Create the struct datatype
	MPI_Type_create_struct(4, block_lengths, offsets, types, &UpdateTimeCorrType);
	MPI_Type_commit(&UpdateTimeCorrType);

	return UpdateTimeCorrType;
}

MPI_Datatype createUpdateIrrForceType() {

	int block_lengths[12] = {1, 1, 5, 3, 3, 12, 12, 1, 1, 1, 1, 1}; // Number of elements in each field
	MPI_Aint offsets[12];
	MPI_Datatype types[12] = {MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
								MPI_UNSIGNED_LONG_LONG, MPI_INT, MPI_DOUBLE,
								MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateIrrForce, pid);
	offsets[1] = offsetof(UpdateIrrForce, newnumberofneighbors);
	offsets[2] = offsetof(UpdateIrrForce, newneighbors);
	offsets[3] = offsetof(UpdateIrrForce, newposition);
	offsets[4] = offsetof(UpdateIrrForce, newvelocity);
	offsets[5] = offsetof(UpdateIrrForce, airr);
	offsets[6] = offsetof(UpdateIrrForce, atot);
	offsets[7] = offsetof(UpdateIrrForce, newcurrentblock_irr);
	offsets[8] = offsetof(UpdateIrrForce, timelevel_irr);
	offsets[9] = offsetof(UpdateIrrForce, timestep_irr);
	offsets[10] = offsetof(UpdateIrrForce, timeblock_irr);
	offsets[11] = offsetof(UpdateIrrForce, nextblock_irr);

	// Create the struct datatype
	MPI_Type_create_struct(12, block_lengths, offsets, types, &UpdateIrrForceType);
	MPI_Type_commit(&UpdateIrrForceType);

	return UpdateIrrForceType;
}

MPI_Datatype createUpdateBinaryType() {

	int block_lengths[6] = {1, 1, 3, 3, 1, 1}; // Number of elements in each field
	MPI_Aint offsets[6];
	MPI_Datatype types[6] = {MPI_INT, MPI_LONG_LONG_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateBinary, pid);
	offsets[1] = offsetof(UpdateBinary, binary_state);
	offsets[2] = offsetof(UpdateBinary, position);
	offsets[3] = offsetof(UpdateBinary, velocity);
	offsets[4] = offsetof(UpdateBinary, mass);
	offsets[5] = offsetof(UpdateBinary, currenttime_irr);

	// Create the struct datatype
	MPI_Type_create_struct(6, block_lengths, offsets, types, &UpdateBinaryType);
	MPI_Type_commit(&UpdateBinaryType);

	return UpdateBinaryType;
}

MPI_Datatype createUpdateFBTermType() {

	int block_lengrhs[21] = {
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		12, 12, 1, MaxNumNeighbor, 3, 3
	}; // Number of elements in each field
	MPI_Aint offsets[21];
	MPI_Datatype types[21] = {
		MPI_INT, MPI_LONG_LONG_INT, MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE,
		MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, MPI_INT, MPI_DOUBLE,
		MPI_UNSIGNED_LONG_LONG, MPI_INT, MPI_DOUBLE,
		MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
		MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE
	}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateFBTerm, pid);
	offsets[1] = offsetof(UpdateFBTerm, binary_state);
	offsets[2] = offsetof(UpdateFBTerm, currentblock_irr);
	offsets[3] = offsetof(UpdateFBTerm, currenttime_irr);
	offsets[4] = offsetof(UpdateFBTerm, currentblock_reg);
	offsets[5] = offsetof(UpdateFBTerm, currenttime_reg);
	offsets[6] = offsetof(UpdateFBTerm, newcurrentblock_irr);
	offsets[7] = offsetof(UpdateFBTerm, nextblock_irr);
	offsets[8] = offsetof(UpdateFBTerm, timelevel_irr);
	offsets[9] = offsetof(UpdateFBTerm, timestep_irr);
	offsets[10] = offsetof(UpdateFBTerm, timeblock_irr);
	offsets[11] = offsetof(UpdateFBTerm, timelevel_reg);
	offsets[12] = offsetof(UpdateFBTerm, timestep_reg);
	offsets[13] = offsetof(UpdateFBTerm, timeblock_reg);
	offsets[14] = offsetof(UpdateFBTerm, radiusofneighbor);
	offsets[15] = offsetof(UpdateFBTerm, airr);
	offsets[16] = offsetof(UpdateFBTerm, areg);
	offsets[17] = offsetof(UpdateFBTerm, numberofneighbors);
	offsets[18] = offsetof(UpdateFBTerm, neighbors);
	offsets[19] = offsetof(UpdateFBTerm, position);
	offsets[20] = offsetof(UpdateFBTerm, velocity);

	// Create the struct datatype
	MPI_Type_create_struct(21, block_lengrhs, offsets, types, &UpdateFBTermType);
	MPI_Type_commit(&UpdateFBTermType);

	return UpdateFBTermType;
}

MPI_Datatype createUpdateNewCMType() {

	int block_lengths[23] = {
		1, 1, 5, 3, 3, 1, 1, 1, 1,
		1, 1, 1, 1, 
		1, 1, 1, 1, 1, 1, 
		1, MaxNumNeighbor, 12, 12
	}; // Number of elements in each field
	MPI_Aint offsets[23];
	MPI_Datatype types[23] = {
		MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
		MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, 
		MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG, MPI_INT, MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG, MPI_INT,
		MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE
	}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateNewCM, PID);
	offsets[1] = offsetof(UpdateNewCM, numberofmembers);
	offsets[2] = offsetof(UpdateNewCM, members);
	offsets[3] = offsetof(UpdateNewCM, position);
	offsets[4] = offsetof(UpdateNewCM, velocity);
	offsets[5] = offsetof(UpdateNewCM, mass);
	offsets[6] = offsetof(UpdateNewCM, radiusofneighbor);
	offsets[7] = offsetof(UpdateNewCM, currenttime_irr);
	offsets[8] = offsetof(UpdateNewCM, currenttime_reg);
	offsets[9] = offsetof(UpdateNewCM, currentblock_irr);
	offsets[10] = offsetof(UpdateNewCM, currentblock_reg);
	offsets[11] = offsetof(UpdateNewCM, newcurrentblock_irr);
	offsets[12] = offsetof(UpdateNewCM, nextblock_irr);
	offsets[13] = offsetof(UpdateNewCM, timestep_irr);
	offsets[14] = offsetof(UpdateNewCM, timeblock_irr);
	offsets[15] = offsetof(UpdateNewCM, timelevel_irr);
	offsets[16] = offsetof(UpdateNewCM, timestep_reg);
	offsets[17] = offsetof(UpdateNewCM, timeblock_reg);
	offsets[18] = offsetof(UpdateNewCM, timelevel_reg);
	offsets[19] = offsetof(UpdateNewCM, numberofneighbors);
	offsets[20] = offsetof(UpdateNewCM, neighbors);
	offsets[21] = offsetof(UpdateNewCM, airr);
	offsets[22] = offsetof(UpdateNewCM, areg);

	// Create the struct datatype
	MPI_Type_create_struct(23, block_lengths, offsets, types, &UpdateNewCMType);
	MPI_Type_commit(&UpdateNewCMType);

	return UpdateNewCMType;
}

MPI_Datatype createUpdateRegCudaType() {

	int block_lengths[8] = {1, 3, 3, 6, 12, 12, 1, MaxNumNeighbor}; // Number of elements in each field
	MPI_Aint offsets[8];
	MPI_Datatype types[8] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateRegCuda, pid);
	offsets[1] = offsetof(UpdateRegCuda, newposition);
	offsets[2] = offsetof(UpdateRegCuda, newvelocity);
	offsets[3] = offsetof(UpdateRegCuda, airr);
	offsets[4] = offsetof(UpdateRegCuda, areg);
	offsets[5] = offsetof(UpdateRegCuda, atot);
	offsets[6] = offsetof(UpdateRegCuda, newnumberofneighbors);
	offsets[7] = offsetof(UpdateRegCuda, newneighbors);

	// Create the struct datatype
	MPI_Type_create_struct(8, block_lengths, offsets, types, &UpdateRegCudaType);
	MPI_Type_commit(&UpdateRegCudaType);

	return UpdateRegCudaType;
}

MPI_Datatype createUpdateRegCudaUpdateType() {

	int block_lengths[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Number of elements in each field
	MPI_Aint offsets[11];
	MPI_Datatype types[11] = {MPI_INT, MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE, 
		MPI_INT, MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG,
		MPI_INT, MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG,
		MPI_DOUBLE, MPI_UNSIGNED_LONG_LONG
	}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateRegCudaUpdate, pid);
	offsets[1] = offsetof(UpdateRegCudaUpdate, currentblock_reg);
	offsets[2] = offsetof(UpdateRegCudaUpdate, currenttime_reg);
	offsets[3] = offsetof(UpdateRegCudaUpdate, timelevel_reg);
	offsets[4] = offsetof(UpdateRegCudaUpdate, timestep_reg);
	offsets[5] = offsetof(UpdateRegCudaUpdate, timeblock_reg);
	offsets[6] = offsetof(UpdateRegCudaUpdate, timelevel_irr);
	offsets[7] = offsetof(UpdateRegCudaUpdate, timestep_irr);
	offsets[8] = offsetof(UpdateRegCudaUpdate, timeblock_irr);
	offsets[9] = offsetof(UpdateRegCudaUpdate, radiusofneighbor);
	offsets[10] = offsetof(UpdateRegCudaUpdate, nextblock_irr);

	// Create the struct datatype
	MPI_Type_create_struct(11, block_lengths, offsets, types, &UpdateRegCudaUpdateType);
	MPI_Type_commit(&UpdateRegCudaUpdateType);

	return UpdateRegCudaUpdateType;
}
#ifdef SEVN
MPI_Datatype createUpdateSEVN0Type() {

	int block_lengths[4] = {1, 1, 1, 1}; // Number of elements in each field
	MPI_Aint offsets[4];
	MPI_Datatype types[4] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateSEVN0, pid);
	offsets[1] = offsetof(UpdateSEVN0, particletype);
	offsets[2] = offsetof(UpdateSEVN0, mass);
	offsets[3] = offsetof(UpdateSEVN0, radius);

	// Create the struct datatype
	MPI_Type_create_struct(4, block_lengths, offsets, types, &UpdateSEVN0Type);
	MPI_Type_commit(&UpdateSEVN0Type);	

	return UpdateSEVN0Type;
}

MPI_Datatype createUpdateSEVN1Type() {

	int block_lengths[8] = {1, 1, 1, 1, 3, 1, 3, 1}; // Number of elements in each field
	MPI_Aint offsets[8];
	MPI_Datatype types[8] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
		MPI_DOUBLE, MPI_LONG_LONG_INT, MPI_DOUBLE, MPI_CHAR
	}; // Match the types in the struct

	// Calculate offsets
	offsets[0] = offsetof(UpdateSEVN1, pid);
	offsets[1] = offsetof(UpdateSEVN1, dm);
	offsets[2] = offsetof(UpdateSEVN1, mass);
	offsets[3] = offsetof(UpdateSEVN1, radius);
	offsets[4] = offsetof(UpdateSEVN1, velocity);
	offsets[5] = offsetof(UpdateSEVN1, binary_state);
	offsets[6] = offsetof(UpdateSEVN1, a_spin);
	offsets[7] = offsetof(UpdateSEVN1, isactive);

	// Create the struct datatype
	MPI_Type_create_struct(8, block_lengths, offsets, types, &UpdateSEVN1Type);
	MPI_Type_commit(&UpdateSEVN1Type);	

	return UpdateSEVN1Type;
}
#endif
#endif