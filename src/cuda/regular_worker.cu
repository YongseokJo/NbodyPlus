#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <cuda.h>  // CUDA Driver API
#include <cuda_runtime.h>
#include <unordered_set>
#include <mpi.h>
#include <vector>
#include <unordered_set>
#include <cstring>

#include <random>
#include <algorithm>

#include "../def.h"
#include "../global.h"
#include "../Queue.h"
#include "cuda_defs.h"
#include "cuda_kernels.h"
#include "cuda_routines.h"
#include "QueueScheduler.h"


#ifndef NEIGHBOR_COUNT_TAG
#define NEIGHBOR_COUNT_TAG 100
#endif
#ifndef NEIGHBOR_LIST_TAG
#define NEIGHBOR_LIST_TAG 101
#endif

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

// const int ROOT = 0;

// static variables for GPU acceleration
static bool first   = true;
static int j_size = 0;
static int i_size = 0;
static int i_size_loop = 0;
static bool is_open = false;
cudaStream_t stream; //Maximum 4 GPUs

static CUDA_REAL *d_ptcl_i = nullptr; // x,v,m
static CUDA_REAL *d_ptcl_j = nullptr; // x,v,m 
static CUDA_REAL *d_r2 = nullptr; // radius^2
static CUDA_REAL *d_result = nullptr; // result array for reduction
static CUDA_REAL *h_result = nullptr; // host result array for reduction
static CUDA_REAL *d_result_block = nullptr; // result block for reduction
static CUDA_REAL *h_ptcl_j, *h_ptcl_i, *h_r2;
static int *h_indices;
static int *d_target = nullptr; // target indices
static int *d_indices = nullptr; // indices for particles
static int *d_num_neighbor = nullptr; // number of neighbors for each target
static int *h_num_neighbor = nullptr; // host number of neighbors for each target
static int *d_neighbor_block = nullptr; // neighbor block indices
static int *h_neighbor = nullptr; // host neighbor list
static int *d_neighbor = nullptr; // device neighbor list
static int *d_num_neighbor_block = nullptr;
// static int deviceCount = 0;
int deviceCount;

// void RegularWorker(int NumTargetTotal, int Jstart, int Jend, int gpu_id);
// void AllocateDeviceMemory(int N_i, int N_j, int gpu_id);
void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList, int ListSize);
void RegularRoot(int NumTargetTotal, CUDA_REAL Acceleration[], int NumNeighbor[], int *NeighborList);
void SetSize(int N_i, int N_j);
void RegAccelerationWorkThread(int TargetStart, int NumTarget, int NumTargetTotal, int JStart, int Jend, int gpu_id, cudaStream_t stream);


void RegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueScheduler &queue_scheduler){
	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];
    
	int *ACListReceive, *NumNeighborReceive;
	CUDA_REAL *Acceleration;
	int J_start, J_end;

	Particle *ptcl;
	double new_time = NextRegTimeBlock*time_step;  // next regular time

	Acceleration		= new CUDA_REAL[6*ListSize];
	NumNeighborReceive  = new int[ListSize];
	ACListReceive		= new int[ListSize * MaxNumNeighbor];

	std::memset(Acceleration, 0, ListSize * 6 * sizeof(CUDA_REAL));

	Queue queue;
	sendAllParticlesToGPU(new_time, RegularList, IndexList, ListSize);

	// queue_scheduler.initialize(RegCal);
	queue.task = RegCal;
	queue.next_time = -1;
	int* int_lists = new int[3];
	for (int p = 1; p <= deviceCount; p++) {
		int gpu_id = p - 1;
		J_start = (gpu_id * NumberOfParticle + deviceCount - 1) / deviceCount;
		J_end = ((gpu_id + 1) * NumberOfParticle + deviceCount - 1) / deviceCount;
		J_end = std::min(J_end, NumberOfParticle);

		queue.pid = gpu_id;
		MPI_Send(&queue, 1, QueueType, p, QUEUE_TAG, MPI_COMM_WORLD);

		int_lists[0] = ListSize;
		int_lists[1] = J_start;
		int_lists[2] = J_end;
		MPI_Send(int_lists, 3, MPI_INT, p, 1, MPI_COMM_DEVICE);
	}
	delete [] int_lists;
	// Start RegularWorker in the queue_scheduler
	// RegularWorker(ListSize, J_start, J_end, gpu_id);
	SetSize(ListSize, NumberOfParticle);
	RegularRoot(ListSize, Acceleration, NumNeighborReceive, ACListReceive);
	int completed = 0;
	MPI_Status status;
	while (completed < deviceCount) {
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int completed_rank = status.MPI_SOURCE;
		int return_value;
		MPI_Recv(&return_value, 1, MPI_INT, completed_rank, TERMINATE_TAG, MPI_COMM_WORLD, &status);
		completed++;
	}

	for (int i=0; i<ListSize; i++) {
		ptcl = &particles[ActiveIndexToOriginalIndex[IndexList[i]]];
		ptcl->NewNumberOfNeighbor = NumNeighborReceive[i];
		std::memcpy(ptcl->NewNeighbors, &ACListReceive[i * MaxNumNeighbor], NumNeighborReceive[i] * sizeof(int));

		for (int dim=0; dim<Dim; dim++) {
			ptcl->a_irr[dim][0]  = static_cast<double>(Acceleration[_six*i+dim]);
			ptcl->a_irr[dim][1] = static_cast<double>(Acceleration[_six*i+dim+3]);
		}
	}


	delete[] IndexList;
	delete[] Acceleration;
	delete[] NumNeighborReceive;
	delete[] ACListReceive;
}


void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList, int ListSize) {

	// Create a vector of indices from 0 to LastParticleIndex
	std::vector<int> indices(global_variable->LastParticleIndex + 1);
	std::iota(indices.begin(), indices.end(), 0);

	// Shuffle the indices randomly
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);

	// variables for saving variables to send to GPU
	
	CUDA_REAL * Radius2;
	//int size = NumberOfParticle;
	int size=0, j=0, i=0, raw_idx=0;
	int N_target = RegularList.size();
	int N_j_max = (NumberOfParticle + deviceCount - 1) / deviceCount; // As N_j can be differnt among devices

	// allocate memory to the temporary variables
	h_ptcl_j     = new CUDA_REAL[N_j_max*7];
	h_ptcl_i     = new CUDA_REAL[N_target*7];
	Radius2  = new CUDA_REAL[N_target];		
	
	Particle *ptcl;
	double Position[3];
	double Velocity[3];

	int J_start, J_end;
	// copy the data of particles to the arrays to be sent
	
	Queue queue;

	queue.task = RegSend;
	queue.next_time = -1;
	// queue.pid = N_j;


	for (int idx: indices) {
		ptcl       = &particles[idx];
		if (!ptcl->isActive) {
			fprintf(stderr, "Skipping inactive particle (%d)\n", ptcl->PID);
			continue;
		}

		if (RegularList.find(idx) != RegularList.end()) {
			IndexList[j] = size;
			Radius2[j] = (CUDA_REAL)ptcl->RadiusOfNeighbor; // mass weight?

			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position, Velocity);
			h_ptcl_i[j] = (CUDA_REAL) Position[0];
			h_ptcl_i[j + N_target] = (CUDA_REAL) Position[1];
			h_ptcl_i[j + 2 * N_target] = (CUDA_REAL) Position[2];
			h_ptcl_i[j + 3 * N_target] = (CUDA_REAL) Velocity[0];
			h_ptcl_i[j + 4 * N_target] = (CUDA_REAL) Velocity[1];
			h_ptcl_i[j + 5 * N_target] = (CUDA_REAL) Velocity[2];
			j++;
		}

		ActiveIndexToOriginalIndex[size] = idx;
		size++;
	}
	assert(NumberOfParticle ==size); // for debugging by EW 2025.1.25

	
	J_start = 0;
	for (int p = 1; p <= deviceCount; p++) {
		int gpu_id = p - 1; // change this in the future
		// int N_j = std::min(NumberOfParticle / deviceCount, (LastParticleIndex + 1) - J_start);
		int N_j = std::min(N_j_max, NumberOfParticle - J_start);
		j = 0; // reset j for each GPU
		raw_idx = 0;

		while (j < N_j) {
			int idx = indices[J_start + raw_idx];
			ptcl = &particles[idx];

			if (ptcl->isActive) {
				h_ptcl_j[j + 6 * N_j] = (CUDA_REAL)ptcl->Mass;
				ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position, Velocity);
				h_ptcl_j[j] 			= (CUDA_REAL) Position[0];
				h_ptcl_j[j + 	 N_j]	= (CUDA_REAL) Position[1];
				h_ptcl_j[j + 2 * N_j]	= (CUDA_REAL) Position[2];
				h_ptcl_j[j + 3 * N_j]	= (CUDA_REAL) Velocity[0];
				h_ptcl_j[j + 4 * N_j]	= (CUDA_REAL) Velocity[1];
				h_ptcl_j[j + 5 * N_j]	= (CUDA_REAL) Velocity[2];
				j++;
			}
			raw_idx++;
		}


		queue.pid = gpu_id;
		MPI_Send(&queue, 1, QueueType, p, QUEUE_TAG, MPI_COMM_WORLD);
		MPI_Send(&N_target, 1, MPI_INT, p, 1010, MPI_COMM_DEVICE);
		MPI_Send(&N_j, 1, MPI_INT, p, 1011, MPI_COMM_DEVICE);
		MPI_Send(h_ptcl_i, 7 * N_target, MPI_CUDA, p, 1001, MPI_COMM_DEVICE);

		// send the “j”‐particles block
		MPI_Send(h_ptcl_j, 7 * N_j, MPI_CUDA, p, 1002, MPI_COMM_DEVICE);
		MPI_Send(IndexList, N_target, MPI_INT, p, 1003, MPI_COMM_DEVICE);
		MPI_Send(Radius2, N_target, MPI_CUDA, p, 1004, MPI_COMM_DEVICE);

		J_start += N_j; // Update J_start for the next GPU
	}

	assert(J_start == NumberOfParticle);

	int completed = 0;
	MPI_Status status;
	while (completed < deviceCount) {
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int completed_rank = status.MPI_SOURCE;
		int return_value;
		MPI_Recv(&return_value, 1, MPI_INT, completed_rank, TERMINATE_TAG, MPI_COMM_WORLD, &status);
		completed++;
	}
	delete[] h_ptcl_j;
	delete[] h_ptcl_i;
	delete[] Radius2;

}



// To Do list
// 1) Add Global Variable: DeviceCount, list of workers that bound to GPUs
// 2) Write Root Routine: SkeletonRegularWorker

// This should be do the same thing as calculateRegAccelerationOnGPU
//int N_target = RegularList.size();

// Return acc, adot
void RegularRoot(
    int NumTargetTotal,
    CUDA_REAL Acceleration[],
    int NumNeighbor[],
    int *NeighborList
    ){
    
    // Dynamically allocate 2D arrays on the heap to avoid stack overflow.
    // These will store the neighbor data received from each worker.
    int** h_num_neighbor_array = new int*[deviceCount];
    int** NeighborList_array = new int*[deviceCount];
    for (int i = 0; i < deviceCount; ++i) {
        // We only need to allocate for workers (ranks 1 to deviceCount-1),
        // but allocating for all is simpler and safer.
        h_num_neighbor_array[i] = new int[i_size];
        NeighborList_array[i] = new int[i_size * MaxNumNeighbor];
    }


    // The root loops through the particle chunks, receiving and aggregating data for each.
    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += i_size_loop) {
        int NumTarget = std::min(i_size_loop, NumTargetTotal - TargetStart);
        

        // The root process does not compute. It acts as the destination for the reduction.
        MPI_Reduce(MPI_IN_PLACE,
                   Acceleration + TargetStart*_six, // receive buffer
                   _six * NumTarget,
				   MPI_CUDA,
                   MPI_SUM, 
                   0, // root rank
                   MPI_COMM_DEVICE);
        
				
        // Receive neighbor counts and lists from all worker ranks.
        // Workers are assumed to be ranks 1, 2, ..., deviceCount-1.
        for (int p = 1; p <= deviceCount; p++) {
            MPI_Recv(h_num_neighbor_array[p-1],
                     NumTarget,
                     MPI_INT,
                     p, // Receive from worker rank p
                     NEIGHBOR_COUNT_TAG,
                     MPI_COMM_DEVICE,
                     MPI_STATUS_IGNORE);
            MPI_Recv(NeighborList_array[p-1],
                     NumTarget * MaxNumNeighbor,
                     MPI_INT,
                     p, // Receive from worker rank p
                     NEIGHBOR_LIST_TAG,
                     MPI_COMM_DEVICE,
                     MPI_STATUS_IGNORE);
        }

        // --- Aggregate the received neighbor data ---
        for (int j = 0; j < NumTarget; j++) {
            int current_target_idx = TargetStart + j;
            NumNeighbor[current_target_idx] = 0;
            
            // Accumulate neighbor counts from each worker for the current target.
            for (int l = 0; l < deviceCount; l++) {
                NumNeighbor[current_target_idx] += h_num_neighbor_array[l][j];
            }
        }

        for (int k = 0; k < NumTarget; k++) {
            int offset = 0;
            // Merge neighbor lists from all workers for the current target.
            for (int l = 0; l < deviceCount; l++) {
                int count = h_num_neighbor_array[l][k];
                if (offset + count > MaxNumNeighbor) {
                    fprintf(stderr, "ERROR: Sum of neighbors exceeds MaxNumNeighbor for target %d!\n", k);
                    // Handle error, e.g. break or throw
                }
                memcpy(&NeighborList[(TargetStart + k) * MaxNumNeighbor + offset],
                       &NeighborList_array[l][k * MaxNumNeighbor],
                       count * sizeof(int));
                offset += count;
            }
        }
        
    } // end of TargetStart loop

    // --- Free the dynamically allocated memory ---
    for (int i = 0; i < deviceCount; ++i) {
        delete[] h_num_neighbor_array[i];
        delete[] NeighborList_array[i];
    }
    delete[] h_num_neighbor_array;
    delete[] NeighborList_array;
}

void _RegularWorker(int NumTargetTotal, int Jstart, int Jend, int gpu_id){
	// Define h_result, NeighborList, h_num_neighbor_array
	int NumTarget;

    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += i_size_loop) {
		NumTarget = std::min(i_size_loop, NumTargetTotal-TargetStart);

		RegAccelerationWorkThread(TargetStart, NumTarget, NumTargetTotal, Jstart, Jend, gpu_id, stream); //, h_result, NeighborList, h_num_neighbor_array
		cudaStreamSynchronize(stream);
		// Contribute h_result to the global sum on root

		if (MPI_COMM_DEVICE != MPI_COMM_NULL) {
			int local_size, local_rank;
			MPI_Comm_rank(MPI_COMM_DEVICE, &local_rank);
			MPI_Comm_size(MPI_COMM_DEVICE, &local_size);
		}
        MPI_Reduce(h_result,
                   nullptr,
                   _six * NumTarget,
				   MPI_CUDA,
                   MPI_SUM,
                   0,
                   MPI_COMM_DEVICE);

        // Send neighbor counts to root
        MPI_Send(h_num_neighbor,
                 NumTarget,
                 MPI_INT,
                 0,
                 NEIGHBOR_COUNT_TAG,
                 MPI_COMM_DEVICE);

        // Send neighbor lists to root
        MPI_Send(h_neighbor,
                 NumTarget * MaxNumNeighbor,
                 MPI_INT,
                 0,
                 NEIGHBOR_LIST_TAG,
                 MPI_COMM_DEVICE);
        // Ensure GPU work is complete

#define UNUSE
#ifdef UNUSE
        // =================================================================
        // START: CPU VERIFICATION OF GPU NEIGHBOR LIST
        // =================================================================
        // This block recalculates the neighbor list on the CPU and compares
        // it against the list generated by the GPU. This is for debugging.
        // It compares the lists of *shuffled indices*.
        int NumJ = Jend - Jstart;
        bool verification_passed = true;

        for (int i = 0; i < NumTarget; ++i) {
            // This is the target particle we are checking
            int current_target_global_idx = TargetStart + i;

            // Get target particle's data
            CUDA_REAL i_px = h_ptcl_i[current_target_global_idx];
            CUDA_REAL i_py = h_ptcl_i[current_target_global_idx + NumTargetTotal];
            CUDA_REAL i_pz = h_ptcl_i[current_target_global_idx + 2 * NumTargetTotal];
            CUDA_REAL i_r2 = h_r2[current_target_global_idx];
            int target_original_shuffled_idx = h_indices[current_target_global_idx];

            // This vector will store the neighbors found by the CPU
            std::vector<int> cpu_neighbors;
            std::vector<CUDA_REAL> cpu_radii;

            // Iterate over all source particles assigned to this worker
            for (int j_local = 0; j_local < NumJ; ++j_local) {
                int j_global = Jstart + j_local; // The shuffled index of the source particle

                // Get source particle's position
                CUDA_REAL j_px = h_ptcl_j[j_local];
                CUDA_REAL j_py = h_ptcl_j[j_local + NumJ];
                CUDA_REAL j_pz = h_ptcl_j[j_local + 2 * NumJ];

                // Calculate squared distance
                CUDA_REAL dx = j_px - i_px;
                CUDA_REAL dy = j_py - i_py;
                CUDA_REAL dz = j_pz - i_pz;
                CUDA_REAL dist_sq = dx * dx + dy * dy + dz * dz;

                // Avoid self-interaction
                if (target_original_shuffled_idx == j_global) {
                    continue;
                }

                // Check if it's a neighbor (using the correct < condition)
                if (dist_sq < i_r2) {
                    cpu_neighbors.push_back(j_global);
					cpu_radii.push_back(dist_sq / i_r2);
                }
            }

            // Now, compare the CPU-generated list with the GPU-generated list
            int gpu_neighbor_count = h_num_neighbor[i];
            int* gpu_neighbor_list_start = &h_neighbor[i * MaxNumNeighbor];

            if (cpu_neighbors.size() != (size_t)gpu_neighbor_count) {
                fprintf(stderr, "[GPU VERIFY FAIL] Rank %d, Target %d: Neighbor count mismatch! CPU: %zu, GPU: %d\n",
                        gpu_id, target_original_shuffled_idx, cpu_neighbors.size(), gpu_neighbor_count);
                verification_passed = false;
				for (size_t l = 0; l < cpu_neighbors.size(); ++l) {
					fprintf(stderr, "%d, ", cpu_neighbors[l]);
				}
				fprintf(stderr, "\n");
				for (int l = 0; l < gpu_neighbor_count; ++l) {
					fprintf(stderr, "%d, ", gpu_neighbor_list_start[l]);
				}
				fprintf(stderr, "\n");
                continue; // No point in comparing lists if counts differ
            }

            // Sort both lists to perform a consistent comparison
            std::sort(cpu_neighbors.begin(), cpu_neighbors.end());
            std::sort(gpu_neighbor_list_start, gpu_neighbor_list_start + gpu_neighbor_count);

            for (size_t k = 0; k < cpu_neighbors.size(); ++k) {
                if (cpu_neighbors[k] != gpu_neighbor_list_start[k]) {
                    fprintf(stderr, "[GPU VERIFY FAIL] Rank %d, Target %d: Mismatch at neighbor index %zu. (%d, %d) CPU found %d, GPU found %d.\n",
                            gpu_id, target_original_shuffled_idx, k, cpu_neighbors.size(), gpu_neighbor_count, cpu_neighbors[k], gpu_neighbor_list_start[k]);
                    verification_passed = false;		
                    break; // Found a mismatch, no need to check further for this particle
                }
            }
        }

        if (verification_passed) {
            fprintf(stdout, "[GPU VERIFY OK] Rank %d: Neighbor lists for %d targets (starting %d) are correct.\n",
                    gpu_id, NumTarget, TargetStart);
        }
        // =================================================================
        // END: CPU VERIFICATION
        // =================================================================
#endif
	}
}



void RegAccelerationWorkThread(int TargetStart, int NumTarget, int NumTargetTotal, int JStart, int Jend, int gpu_id, cudaStream_t stream){
	// int gpu_id = ...; // Get the GPU ID from myrank. I assume that it starts from 0 to ngpus-1

	// 1) Launch a kernel on each GPU with an offset
	cudaSetDevice(gpu_id);
	dim3 gridDim2(NumTarget, 1);
	dim3 blockDim2(GridDimY, 1);

	int NumJ = Jend - JStart;

	dim3 blockDim(BatchSize, 1, 1);
	dim3 gridDim(
		(NumTarget + BatchSize + blockDim.x - 1) / blockDim.x, 
		GridDimY
	);

	// We send h_ptcl_i to d_ptcl_i for all N_i, however, due to possible memory constraint, we only use NumTarget for each execution
	compute_forces_mpi<<<gridDim, blockDim, 0, stream>>>(
		d_ptcl_i,
		d_r2,
		d_ptcl_j,
		d_result_block,
		d_indices,
		NumTarget,
		NumJ,
		NumTargetTotal,
		d_neighbor_block,
		d_num_neighbor_block,
		TargetStart,   // i_start
		JStart // j_start
	);
	dim3 blockDim3(16, 6);       // 16 threads along X, 6 along Y
	dim3 gridDim3((NumTarget+15)/16, 1);
	reduce_forces_kernel<<<gridDim3, blockDim3>>>(d_result_block, d_result, GridDimY, NumTarget);
	
	gather_neighbor<<<gridDim2, blockDim2, 0, stream>>>(
		d_neighbor_block, 
		d_num_neighbor_block, 
		d_neighbor,
		NumTarget
	);

	gather_numneighbor<<<gridDim2, blockDim2, 0, stream>>>(
		d_num_neighbor_block, 
		d_num_neighbor,
		NumTarget
	);

	toHost(h_result,
			d_result,
			_six * NumTarget,
			stream);

	toHost(h_neighbor,
			d_neighbor,
			NumTarget * MaxNumNeighbor,
			stream);
	
	toHost(h_num_neighbor,
			d_num_neighbor,
			NumTarget,
			stream);

	cudaStreamSynchronize(stream);
}



void _SendToDeviceMPI(
		int N_i,
		int N_j,
		int gpu_id
		){

	cudaError_t cudaStatus;
	cudaSetDevice(gpu_id);

	toDevice(h_ptcl_j, d_ptcl_j, _seven*N_j, stream);
	toDevice(h_ptcl_i, d_ptcl_i, _seven*N_i, stream);
	toDevice(h_r2, d_r2, N_i, stream);
	toDevice(h_indices, d_indices, N_i, stream);

	cudaDeviceSynchronize();
}

void SetSize(int N_i, int N_j) {
	// size of device memory for i and j
	j_size = new_size(N_j);
	i_size = new_size(N_i);
	i_size_loop = i_size; // temporary
}

void _AllocateDeviceMemory(
	int N_i,
	int N_j,
	int gpu_id
	){

	assert(N_i <= nbodymax);
	assert(N_j <= nbodymax);
	cudaError_t cudaStatus;
	cudaSetDevice(gpu_id);

	if ((first) || (new_size(N_j) > j_size) || (new_size(N_i) > i_size)) {
		//varaiable_size should be the number of j, and target size be the number of i
		SetSize(N_i, N_j);
		// size of device memory for i and j

		if (!first) {
			my_free(h_ptcl_i, d_ptcl_i);
			my_free(h_ptcl_j, d_ptcl_j);
			my_free(h_result, d_result);
			my_free(h_num_neighbor , d_num_neighbor);
			cudaFree(d_target);
			my_free(h_r2, d_r2);
			my_free(h_indices, d_indices);
			cudaFree(d_neighbor_block);
			my_free(h_neighbor, d_neighbor);
			cudaFree(d_result_block);
			cudaFree(d_num_neighbor_block);
		}
		else {
			first = false;
		}
		my_allocate(&h_ptcl_i, &d_ptcl_i, _seven*i_size); // x,v,m
		my_allocate(&h_ptcl_j, &d_ptcl_j, _seven*j_size); // x,v,m
		my_allocate(&h_result, &d_result,           _six*i_size);
		my_allocate(&h_num_neighbor, &d_num_neighbor, GridDimY * i_size);
		cudaMalloc((void**)&d_target,  i_size * sizeof(int));
		my_allocate(&h_r2, &d_r2, i_size);
		my_allocate(&h_indices, &d_indices, i_size);
		cudaMalloc((void**)&d_neighbor_block, GridDimY * NNB_per_block * i_size_loop * sizeof(int));
		my_allocate(&h_neighbor, &d_neighbor    , MaxNumNeighbor * i_size);
		cudaMalloc((void**)&d_result_block, _six * GridDimY * i_size_loop * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_num_neighbor_block, GridDimY * i_size_loop * sizeof(CUDA_REAL));

	} //end of if (first) || (new_size(NNB) > j_size)

	cudaDeviceSynchronize();

	MPI_Recv(h_ptcl_i, 7 * N_i, MPI_CUDA, ROOT, 1001, MPI_COMM_DEVICE, MPI_STATUS_IGNORE);
	MPI_Recv(
		h_ptcl_j,                // buffer
		7 * N_j,                 // count (7 fields per j‐particle)
		MPI_CUDA,
		ROOT,
		1002,
		MPI_COMM_DEVICE,
		MPI_STATUS_IGNORE
	);
	MPI_Recv(
		h_indices,               // buffer
		N_i,                // one int per target
		MPI_INT,
		ROOT,
		1003,
		MPI_COMM_DEVICE,
		MPI_STATUS_IGNORE
	);
	MPI_Recv(
		h_r2,
		N_i,
		MPI_CUDA,
		ROOT,
		1004,
		MPI_COMM_DEVICE,
		MPI_STATUS_IGNORE
	);
}

void _InitializeGPU(int gpu_id){
	
	cudaGetDeviceCount(&deviceCount);
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, gpu_id);

	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);

	cudaSetDevice(gpu_id);
	cudaStreamCreate(&stream);
	cudaFree(nullptr);

	// Use CUDA Driver API to get the device associated with the current context
    CUdevice cuDev;
    CUcontext context;
    CUresult resCtx = cuCtxGetCurrent(&context); 
    if ((resCtx == CUDA_SUCCESS) && (context != nullptr)) {
        if (cuCtxGetDevice(&cuDev) == CUDA_SUCCESS) {
            int devId = (int)cuDev;
            // Check if devId is valid
            if (devId < 0 || devId >= deviceCount) {
                std::cerr << "Invalid device ID from context: " << devId << std::endl;
            }
        }
    } else {
        std::cerr << "Failed to get CUDA context on root processor. "
                  << "cuCtxGetCurrent returned: " << resCtx << std::endl;
    }
}







extern "C" {

	void RegularWorker(int NumTargetTotal, int Jstart, int Jend, int gpu_id) {
		_RegularWorker(NumTargetTotal, Jstart, Jend, gpu_id);  // internal function
	}

	void AllocateDeviceMemory(int N_i, int N_j, int gpu_id) {
		_AllocateDeviceMemory(N_i, N_j, gpu_id);
	}

	void SendToDeviceMPI(int N_i, int N_j, int gpu_id) {
		_SendToDeviceMPI(N_i, N_j, gpu_id);
	}
	void InitializeGPU(int gpu_id){
		_InitializeGPU(gpu_id);
	}
}



