#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <cuda.h>  // CUDA Driver API
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "../def.h"
#include "cuda_defs.h"
#include "cuda_kernels.h"
#include "cuda_routines.h"

#ifndef NEIGHBOR_COUNT_TAG
#define NEIGHBOR_COUNT_TAG 100
#endif
#ifndef NEIGHBOR_LIST_TAG
#define NEIGHBOR_LIST_TAG 101
#endif

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

// static variables for GPU acceleration
static int deviceCount;
static bool first   = true;
static int j_size = 0;
static int i_size = 0;
static bool is_open = false;
cudaStream_t stream; //Maximum 4 GPUs

static CUDA_REAL *d_ptcl_i = nullptr; // x,v,m
static CUDA_REAL *d_ptcl_j = nullptr; // x,v,m 
static CUDA_REAL *d_r2 = nullptr; // radius^2
static CUDA_REAL *d_diff = nullptr; // difference array for reduction
static CUDA_REAL *d_result = nullptr; // result array for reduction
static CUDA_REAL *h_result = nullptr; // host result array for reduction
static CUDA_REAL *d_result_block = nullptr; // result block for reduction
static int *d_target = nullptr; // target indices
static int *d_indices = nullptr; // indices for particles
static int *d_num_neighbor = nullptr; // number of neighbors for each target
static int *h_num_neighbor = nullptr; // host number of neighbors for each target
static int *d_neighbor_block = nullptr; // neighbor block indices
static int *h_neighbor = nullptr; // host neighbor list
static int *d_neighbor = nullptr; // device neighbor list


// To Do list
// 1) Add Global Variable: DeviceCount, list of workers that bound to GPUs
// 2) Write Root Routine: SkeletonRegularWorker

// This should be do the same thing as calculateRegAccelerationOnGPU
//int N_target = RegularList.size();

void calculateRegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueScheduler &queue_scheduler){
	int ListSize = RegularList.size();
	int *IndexList = new int[ListSize];

	int *ACListReceive;
	int *NumNeighborReceive;

	Particle *ptcl;
	double new_time = NextRegTimeBlock*time_step;  // next regular time

	Acceleration		= new CUDA_REAL[6*ListSize];
	NumNeighborReceive  = new int[ListSize];
	ACListReceive = new int[ListSize * MaxNumNeighbor];

	std::memset(Acceleration, 0, ListSize * 6 * sizeof(CUDA_REAL));

	// Run this in worker's routine in queue_scheduler
	// two workers launch sendAllParitclesToGPU using MPI
	// replace this with the queue_scheduler

	/* // This is the original code that uses MPI to send particles to GPU
	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	assert(mpi_size > deviceCount);  // at least ranks 0 (root), 1, 2 (workers)
	if (mpi_rank > 0 && mpi_rank <= deviceCount) {
		int gpu_id = mpi_rank - 1; // temporarily
		int N_start = (gpu_id * ListSize) / deviceCount;
		int N_end = ((gpu_id + 1) * ListSize) / deviceCount;
		sendAllParticlesToGPU(new_time, RegularList, IndexList, N_start, N_end, gpu_id);
		// But this code, only processeses bound to GPUs participate in to the preprocessing
	}
	// synchronize all ranks before proceeding
	MPI_Barrier(MPI_COMM_WORLD);
	*/

	int gpu_id = 0;
	int N_start = (gpu_id * ListSize) / deviceCount;
	int N_end = ((gpu_id + 1) * ListSize) / deviceCount;

	queue_scheduler.initialize(RegSend);
	// Start sendAllParticleToGPU in the queue_scheduler
	sendAllParticlesToGPU(new_time, RegularList, IndexList, N_start, N_end, gpu_id);

	queue_scheduler.initialize(RegCal);
	// Start RegularWorker in the queue_scheduler
	RegularWorker(ListSize, gpu_id);
	// But this code, only processeses bound to GPUs participate in to the preprocessing

	RegularRoot(ListSize, Acceleration, NumNeighborReceive, ACListReceive);

	for (int i=0; i<ListSize; i++) {
		ptcl = &particles[ActiveIndexToOriginalIndex[IndexList[i]]];
		ptcl->NewNumberOfNeighbor = NumNeighborReceive[i];
		std::memcpy(ptcl->NewNeighbors, &ACListReceive[i * MaxNumNeighbor], NumNeighborReceive[i] * sizeof(int));

		for (int dim=0; dim<Dim; dim++) {
			ptcl->a_irr[dim][0]  = static_cast<double>(Acceleration[_six*i+dim]);
			ptcl->a_irr[dim][1] = static_cast<double>(Acceleration[_six*i+dim+3]);
		}
	}

	queue_scheduler.initialize(RegCuda);
	queue_scheduler.takeQueueRegularList(RegularList);
	do
	{
		queue_scheduler.assignQueueAutoRegularList();
		queue_scheduler.runQueueAuto();
		queue_scheduler.waitQueue(0); // blocking wait
	} while (queue_scheduler.isComplete());

	delete[] IndexList;
	delete[] Acceleration;
	delete[] NumNeighborReceive;
	delete[] ACListReceive;
}

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
    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += i_size) {
        int NumTarget = std::min(i_size, NumTargetTotal - TargetStart);
        
        // The root process does not compute. It acts as the destination for the reduction.
        MPI_Reduce(nullptr, // sendbuf is ignored on the root if not using MPI_IN_PLACE
                   Acceleration + TargetStart*_six, // receive buffer
                   _six * NumTarget,
                   MPI_FLOAT, 
                   MPI_SUM, 
                   0, // root rank
                   MPI_COMM_WORLD);
        
        // Receive neighbor counts and lists from all worker ranks.
        // Workers are assumed to be ranks 1, 2, ..., deviceCount-1.
        for (int p = 1; p < deviceCount; p++) {
            MPI_Recv(h_num_neighbor_array[p],
                     NumTarget,
                     MPI_INT,
                     p, // Receive from worker rank p
                     NEIGHBOR_COUNT_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Recv(NeighborList_array[p],
                     NumTarget * MaxNumNeighbor,
                     MPI_INT,
                     p, // Receive from worker rank p
                     NEIGHBOR_LIST_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        }

        // --- Aggregate the received neighbor data ---
        for (int j = 0; j < NumTarget; j++) {
            int current_target_idx = TargetStart + j;
            NumNeighbor[current_target_idx] = 0;
            
            // Accumulate neighbor counts from each worker for the current target.
            for (int p = 1; p < deviceCount; p++) {
                NumNeighbor[current_target_idx] += h_num_neighbor_array[p][j];
            }
        }

        for (int k = 0; k < NumTarget; k++) {
            int offset = 0;
            // Merge neighbor lists from all workers for the current target.
            for (int p = 1; p < deviceCount; p++) {
                int count = h_num_neighbor_array[p][k];
                if (offset + count > MaxNumNeighbor) {
                    fprintf(stderr, "ERROR: Sum of neighbors exceeds MaxNumNeighbor for target %d!\n", k);
                    // Handle error, e.g. break or throw
                }
                memcpy(&NeighborList[(TargetStart + k) * MaxNumNeighbor + offset],
                       &NeighborList_array[p][k * MaxNumNeighbor],
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

void RegularWorker(int NumTargetTotal, int gpu_id){
	// Define h_result, NeighborList, h_num_neighbor_array
	int NumTarget;

    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += i_size) {
		NumTarget = std::min(i_size, NumTargetTotal-TargetStart);
		
		RegAccelerationWorkThread(TargetStart, NumTarget, gpu_id, stream); //, h_result, NeighborList, h_num_neighbor_array
		cudaStreamSynchronize(stream);
		// Contribute h_result to the global sum on root
        MPI_Reduce(h_result,
                   nullptr,
                   _six * NumTarget,
                   MPI_FLOAT,
                   MPI_SUM,
                   0,
                   MPI_COMM_WORLD);

        // Send neighbor counts to root
        MPI_Send(h_num_neighbor,
                 NumTarget,
                 MPI_INT,
                 0,
                 NEIGHBOR_COUNT_TAG,
                 MPI_COMM_WORLD);

        // Send neighbor lists to root
        MPI_Send(NeighborList,
                 NumTarget * MaxNumNeighbor,
                 MPI_INT,
                 0,
                 NEIGHBOR_LIST_TAG,
                 MPI_COMM_WORLD);

        // Ensure GPU work is complete
	}
}



void RegAccelerationWorkThread(int TargetStart, int NumTarget, int gpu_id, int JStart, cudaStream_t stream){
	// int gpu_id = ...; // Get the GPU ID from myrank. I assume that it starts from 0 to ngpus-1
	int chunkPerGpu = (NNB + deviceCount - 1) / deviceCount;
	int NumTarget;

	// 1) Launch a kernel on each GPU with an offset
	cudaSetDevice(gpu_id);
	dim3 gridDim2(NumTarget, 1);
	dim3 blockDim2(GridDimY, 1);

	int deviceJStart = gpu_id * chunkPerGpu;
	int deviceNumJ = std::min(chunkPerGpu, NNB - deviceJStart);
	if (deviceNumJ <= 0) return; 

	toDevice(h_target_list + TargetStart,
		d_target,
		NumTarget,
		stream);

	dim3 blockDim(BatchSize, 1, 1);
	dim3 gridDim(
		(NumTarget + BatchSize + blockDim.x - 1) / blockDim.x, 
		GridDimY
	);

	// Launch compute_forces on GPU i -> change i into gpu_id
	compute_forces_mpi<<<gridDim, blockDim, 0, stream>>>(
		d_ptcl_i,
		d_r2,
		d_ptcl_j,
		d_result_block,
		d_indices,
		NumTarget, //NNB
		d_neighbor_block,
		d_num_neighbor_block,
		TargetStart,   // i_start
		JStart, // j_start
		NNB
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

	toHost(NeighborList,
			d_neighbor,
			NumTarget * MaxNumNeighbor,
			stream);
	
	toHost(h_num_neighbor,
			d_num_neighbor,
			NumTarget,
			stream);

	cudaStreamSynchronize(stream);

}



void SendToDeviceMPI(
		int _NNB,
		CUDA_REAL h_ptcl_j[],
		CUDA_REAL h_ptcl_i[],
		CUDA_REAL r2_i[],
		int i_indices[],
		cudaStream_t stream,
		int gpu_id
		){

	cudaError_t cudaStatus;
	cudaSetDevice(gpu_id);

	toDevice(h_ptcl_j, d_ptcl_j, _seven*j_size, stream);
	toDevice(h_ptcl_i, d_ptcl_i, _seven*i_size, stream);
	toDevice(r2_i, d_r2, i_size, stream);
	toDevice(i_indices, d_indices, i_size, stream);

	cudaDeviceSynchronize();
}


void AllocateDeviceMemory(
	int N_i,
	int N_j,
	cudaStream_t stream,
	int gpu_id
	){

	assert(N_i <= nbodymax);
	assert(N_j <= nbodymax);
	cudaError_t cudaStatus;
	cudaSetDevice(gpu_id);

	if ((first) || (new_size(N_j) > j_size) || (new_size(N_i) > i_size)) {
		//varaiable_size should be the number of j, and target size be the number of i
		
		j_size = new_size(N_j); // size of device memory for j
		i_size = new_size(N_i);
		fprintf(stderr, "background_size=%d, target_size=%d\n", j_size, i_size);

		if (!first) {
			cudaFree(d_ptcl_i);
			cudaFree(d_ptcl_j);
			my_free(h_result, d_result);
			my_free(h_num_neighbor , d_num_neighbor);
			cudaFree(d_target);
			cudaFree(d_r2);
			cudaFree(d_diff);
			cudaFree(d_indices);
			cudaFree(d_neighbor_block);
			my_free(h_neighbor, d_neighbor);
			cudaFree(d_result_block);
		}
		else {
			first = false;
		}

		cudaMalloc((void**)&d_ptcl_i, _seven*i_size*sizeof(CUDA_REAL) ); // x,v,m
		cudaMalloc((void**)&d_ptcl_j, _seven*j_size*sizeof(CUDA_REAL) ); // x,v,m
		my_allocate(&h_result, &d_result,           _six*i_size);
		my_allocate(&h_num_neighbor, &d_num_neighbor, GridDimY * i_size);
		cudaMalloc((void**)&d_target,  i_size * sizeof(int));
		cudaMalloc((void**)&d_r2, i_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_diff, _six * GridDimY * i_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_indices, i_size * sizeof(int));
		cudaMalloc((void**)&d_neighbor_block, GridDimY * NNB_per_block * i_size * sizeof(int));
		my_allocate(&h_neighbor, &d_neighbor    , MaxNumNeighbor * i_size);
		cudaMalloc((void**)&d_result_block, _six * GridDimY * i_size * sizeof(CUDA_REAL));
	} //end of if (first) || (new_size(NNB) > j_size)

	cudaDeviceSynchronize();
}

void _InitializeDevice(int gpu_id){

	cudaGetDeviceCount(&deviceCount);

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");

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
            std::cout << "[Rank " << MyRank << "] Current device from driver context = " << devId << std::endl;
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

void _OpenDevice(const int irank){
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
	icall = ini = isend = 0;

	_InitializeDevice(irank);
	if(is_open){
		fprintf(stderr, "gpunb: it is already open\n");
		return;
	}
	is_open = true;
}

void _CloseDevice() {
	if(!is_open) {
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;
	cudaError_t error;

	error = cudaGetLastError();
	if (error != cudaSuccess) {
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		// Handle error
	}
}


void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList, int N_start, int N_end, int gpu_id) {

		// Create a vector of indices from 0 to LastParticleIndex
		std::vector<int> indices(global_variable->LastParticleIndex + 1);
		std::iota(indices.begin(), indices.end(), 0);
	
		// Shuffle the indices randomly
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(indices.begin(), indices.end(), g);
	
		// variables for saving variables to send to GPU
		CUDA_REAL *h_ptcl_j, *h_ptcl_i;
		CUDA_REAL * Radius2;
		//int size = NumberOfParticle;
		int size=0, j=0, i=0;
		int N_target = RegularList.size();

		// allocate memory to the temporary variables
		int N_j = N_end - N_start;
		h_ptcl_j     = new CUDA_REAL[N_j*7];
		h_ptcl_i     = new CUDA_REAL[N_target*7];
		Radius2  = new CUDA_REAL[N_target];		
		
		Particle *ptcl;
		double Position[3];
		double Velocity[3];

		// copy the data of particles to the arrays to be sent
		for (int i = N_start; i < N_end; ++i) {
			int idx = indices[i];
			ptcl = &particles[idx];
	
			if (!ptcl->isActive) {
				// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
				continue;
			}
	
			// h_ptcl_j[size + k * NumberOfParticle] = PREDICTED POSITION AND VELCOITY;
			h_ptcl_j[i + 6 * N_j] = (CUDA_REAL)ptcl->Mass;
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position, Velocity);
			h_ptcl_j[i] = (CUDA_REAL) Position[0];
			h_ptcl_j[i + N_j] = (CUDA_REAL) Position[1];
			h_ptcl_j[i + 2 * N_j] = (CUDA_REAL) Position[2];
			h_ptcl_j[i + 3 * N_j] = (CUDA_REAL) Velocity[0];
			h_ptcl_j[i + 4 * N_j] = (CUDA_REAL) Velocity[1];
			h_ptcl_j[i + 5 * N_j] = (CUDA_REAL) Velocity[2];

		}
	

		for (int idx: indices) {
			ptcl       = &particles[idx];
			size++;
			if (!ptcl->isActive) {
				// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
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


		assert(NumberOfParticle == size); // for debugging by EW 2025.1.25
		AllocateDeviceMemory(NumberOfParticle, N_j, stream, gpu_id); // this part cannot be parallelized
		SendToDeviceMPI(&size, h_ptcl_j, h_ptcl_i, Radius2, IndexList, stream, gpu_id); //this part should be parallelized

		delete[] h_ptcl_j;
		delete[] h_ptcl_i;
		delete[] Radius2;
	}
	