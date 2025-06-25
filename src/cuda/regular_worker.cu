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

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif


// To Do list
// 1) Add Global Variable: DeviceCount, list of workers that bound to GPUs
// 2) Write Root Routine: SkeletonRegularWorker

// Return acc, adot
void SkeletonRegularRoot(){

	// Define h_result, NeighborList_array, NeighborList, h_num_neighbor_array


    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += target_size) {
        // How many targets remain in this chunk
        // int bigChunk = std::min(target_size, NumTargetTotal - TargetStart);
		NumTarget = std::min(target_size, NumTargetTotal-TargetStart);

		MPI_Reduce();
		MPI_Gather(NeighborList_array);
		MPI_Gather(NeighborList_array, ...);

		for (int k = 0; k < NumTarget; k++) {
			int offset = 0;
			for (int i = 0; i < deviceCount; i++) {
				int count = h_num_neighbor_array[i][k];
				if (offset + count > MaxNumNeighbor) {
					fprintf(stderr, "ERROR: Sum of neighbors exceeds MaxNumNeighbor for target %d!\n", k);
					// Handle error, e.g. break or throw
				}
				memcpy(&NeighborList[k * MaxNumNeighbor + offset],
					&NeighborList_array[i][k * MaxNumNeighbor],
					count * sizeof(int));
				offset += count; 
			}
		}

		for (int i=0; i<NumTarget; i++) {
			acc[i+TargetStart][0]  = h_result[_six*i];
			acc[i+TargetStart][1]  = h_result[_six*i+1];
			acc[i+TargetStart][2]  = h_result[_six*i+2];
			adot[i+TargetStart][0] = h_result[_six*i+3];
			adot[i+TargetStart][1] = h_result[_six*i+4];
			adot[i+TargetStart][2] = h_result[_six*i+5];
		
	} // end of TargetStart loop
}

void SkeletonRegularWorker(){
	// Define h_result, NeighborList, h_num_neighbor_array


	sendAllParticlesToGPU();
	// This function contains SendToDeviceMPI
	// But this code, only processeses bound to GPUs participate in to the preprocessing


    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += target_size) {
        // How many targets remain in this chunk
        // int bigChunk = std::min(target_size, NumTargetTotal - TargetStart);
		NumTarget = std::min(target_size, NumTargetTotal-TargetStart);

		RegAccelerationWorkThread(TargetStart, NumTarget, gpu_id, stream, h_result, NeighborList, h_num_neighbor_array);
		MPI_Reduce(h_result,      // send buffer
			nullptr,       // receive buffer (only on root)
			_six * NumTarget, // number of elements
			MPI_FLOAT,        // data type
			MPI_SUM,          // operation: sum
			0,                // root rank
			MPI_COMM_WORLD);  // communicator
		
			MPI_Send(h_num_neighbor, // We couldn't use MPI_Reduce here because we need to use this data for gathering NeighborList
			nullptr,
			_six * NumTarget,
			MPI_FLOAT,
			MPI_SUM,
			0,
			MPI_COMM_WORLD);
		
		MPI_Send(NeighborList, ...); // This process would not be trivial
		cudaStreamSynchronize(); // or deviceSynchronize();
	}
}



void RegAccelerationWorkThread(int TargetStart, int NumTarget, int gpu_id, cudaStream_t stream){
	// int gpu_id = ...; // Get the GPU ID from myrank. I assume that it starts from 0 to ngpus-1
	int chunkPerGpu = (NNB + deviceCount - 1) / deviceCount;
	int NumTarget;

	// 1) Launch a kernel on each GPU with an offset
	cudaSetDevice(gpu_id);
	dim3 gridDim2(NumTarget, 1);
	dim3 blockDim2(GridDimY, 1);

	int deviceJStart = gpu_id * chunkPerGpu;
	int deviceNumJ = std::min(chunkPerGpu, NNB - deviceJStart);
	if (deviceNumJ <= 0) break;

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
	compute_forces_mpi<<<gridDim, blockDim, 0, streams[i]>>>(
		d_ptcl_i,
		d_r2,
		d_ptcl_j,
		d_acc,
		d_target_indices,
		NumI, //NNB
		NumJ,
		d_neighbor_block,
		d_num_neighbor_block,
		IStart,   // i_start
		JStart, // j_start
		NNB
	);
	dim3 blockDim3(16, 6);       // 16 threads along X, 6 along Y
	dim3 gridDim3((NumTarget+15)/16, 1);
	reduce_forces_kernel<<<gridDim3, blockDim3>>>(d_diff_array[i], d_result_array[i], GridDimY, NumTarget);

	gather_neighbor<<<gridDim2, blockDim2, 0, streams[i]>>>(
		d_neighbor_block, 
		d_num_neighbor_block, 
		d_neighbor,
		NumTarget
	);

	gather_numneighbor<<<gridDim2, blockDim2, 0, streams[i]>>>(
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
		cudaStream_t stream,
		int gpu_id
		){

	NNB            = _NNB;
	isend++;
	assert(NNB <= nbodymax);
	cudaError_t cudaStatus;

	cudaSetDevice(gpu_id);

	if ((first) || (new_size(NNB) / deviceCount + 1 > j_size )) {
		//varaiable_size should be the number of j, and target size be the number of i
		
		j_size = new_size(NNB) / deviceCount + 1; // +1 to avoid zero size
		i_size = new_size(NNB);
		fprintf(stderr, "background_size=%d, target_size=%d\n", j_size, i_size);

		if (!first) {
			cudaFree(d_ptcl_i);
			cudaFree(d_ptcl_j);
			my_free(h_result, d_result);
			my_free(h_num_neighbor , d_num_neighbor);
			cudaFree(d_target);
			cudaFree(d_r2);
			cudaFree(d_diff);
			my_free(h_neighbor, d_neighbor_block);
		}
		else {
			first = false;
		}

		cudaMalloc((void**)&d_ptcl_j, _seven*j_size*sizeof(CUDA_REAL) ); // x,v,m
		cudaMalloc((void**)&d_ptcl_i, _seven*i_size*sizeof(CUDA_REAL) ); // x,v,m
		my_allocate(&h_result, &d_result,           _six*i_size);
		cudaMalloc((void**)&d_r2, i_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_target,  i_size * sizeof(int));
		cudaMalloc((void**)&d_diff, _six * GridDimY * i_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_neighbor_block, GridDimY * NNB_per_block * i_size * sizeof(int));
		my_allocate(&h_neighbor, &d_neighbor    , MaxNumNeighbor * i_size);
		my_allocate(&h_num_neighbor, &d_num_neighbor, GridDimY * i_size);
	} //end of if (first) || (new_size(NNB) > j_size)

	
	toDevice(h_ptcl, d_ptcl, _seven*NNB, stream);
	cudaDeviceSynchronize();
	toDevice(r2    , d_r2  ,        NNB, stream);
	}
}



void _InitializeDevice(int irank){

	if (MyRank == ROOT) {
	std::cout << "Initializing CUDA ..." << std::endl;
	}
	// Select CUDA device (optional)
	cudaGetDeviceCount(&deviceCount);

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");

	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);

	if (MyRank == ROOT) {
	fprintf(stderr, "# GPU initialization - rank: %d; HOST %s; NGPU %d; device: %d %s\n", irank, hostname,numGPU, devid, prop.name);
	}

	for (int deviceNum = 0; deviceNum < deviceCount; deviceNum++) {
		cudaSetDevice(deviceNum);
		cudaStreamCreate(&streams[deviceNum]);
		initializeCudaAndCublas(&cublasHandles[deviceNum]);	
	    cublasSetStream(cublasHandles[deviceNum], streams[deviceNum]);
        // Force runtime to initialize driver context for this device
        cudaFree(nullptr);
		}

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

    if (MyRank == ROOT) {
        std::cout << "There are " << deviceCount << " GPUs." << std::endl;
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


void sendAllParticlesToGPU(double new_time, std::unordered_set<int> RegularList, int *IndexList) {

	#ifdef MultiNode
	#ifdef PERFORMANCETRACE
		std::chrono::high_resolution_clock::time_point start_point_update;
		std::chrono::high_resolution_clock::time_point end_point_update;
	#endif
	#endif
		
	
		// Create a vector of indices from 0 to LastParticleIndex
		std::vector<int> indices(global_variable->LastParticleIndex + 1);
		std::iota(indices.begin(), indices.end(), 0);
	
		// Shuffle the indices randomly
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(indices.begin(), indices.end(), g);
	
		// variables for saving variables to send to GPU
		CUDA_REAL * h_ptcl_j, h_ptcl_i;
		CUDA_REAL * Radius2;
		//int size = NumberOfParticle;
		int size=0, j=0;
		int N_target = RegularList.size();

		// allocate memory to the temporary variables
		h_ptcl_j     = new CUDA_REAL[NumberOfParticle*7];
		h_ptcl_i     = new CUDA_REAL[N_target*7];
		Radius2  = new CUDA_REAL[N_target];		
		
		Particle *ptcl;
		double Position[3];
		double Velocity[3];

		// copy the data of particles to the arrays to be sent
		for (int idx: indices) {
			ptcl       = &particles[idx];
	
			if (!ptcl->isActive) {
				// fprintf(stdout, "Skipping inactive particle (%d)\n", ptcl->PID);
				continue;
			}
	
			// h_ptcl_j[size + k * NumberOfParticle] = PREDICTED POSITION AND VELCOITY;
			h_ptcl_j[size + 6 * NumberOfParticle] = (CUDA_REAL)ptcl->Mass;
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Position, Velocity);
			h_ptcl_j[size] = (CUDA_REAL) Position[0];
			h_ptcl_j[size + NumberOfParticle] = (CUDA_REAL) Position[1];
			h_ptcl_j[size + 2 * NumberOfParticle] = (CUDA_REAL) Position[2];
			h_ptcl_j[size + 3 * NumberOfParticle] = (CUDA_REAL) Velocity[0];
			h_ptcl_j[size + 4 * NumberOfParticle] = (CUDA_REAL) Velocity[1];
			h_ptcl_j[size + 5 * NumberOfParticle] = (CUDA_REAL) Velocity[2];

			if (RegularList.find(idx) != RegularList.end()) {
				IndexList[j] = size;
				Radius2[j] = (CUDA_REAL)ptcl->RadiusOfNeighbor; // mass weight?
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
		SendToDeviceMPI(&size, h_ptcl_j, h_ptcl_i, Radius2, stream, gpu_id);

		delete[] h_ptcl_j;
		delete[] h_ptcl_i;
		delete[] Radius2;
	
	#ifdef MultiNode
	#ifdef PERFORMANCETRACE
		start_point_update = std::chrono::high_resolution_clock::now();
	#endif
		Queue queue = {UpdateActiveIndexToOriginalIndex, size, -1.0};
		for (int i = 1; i < NumberOfNode; i++) {
			MPI_Send(&queue, 1, QueueType, ranks_update_comm[i], QUEUE_TAG, MPI_COMM_WORLD);
			MPI_Send(ActiveIndexToOriginalIndex, size, MPI_INT, ranks_update_comm[i], 1, MPI_COMM_WORLD);
		}
		int completed = 0;
		MPI_Status status;
		while (completed < NumberOfNode - 1) { // Root node has already completed its job by EW 2025.1.20
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			int completed_rank = status.MPI_SOURCE;
			int return_value;
			MPI_Recv(&return_value, 1, MPI_INT, completed_rank, TERMINATE_TAG, MPI_COMM_WORLD, &status);
			completed++;
		}
	#ifdef PERFORMANCETRACE
		end_point_update = std::chrono::high_resolution_clock::now();
		performance.Update +=
			std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_update - start_point_update).count();
	#endif
	#endif
	}
	