#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <vector>
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



extern int MyRank;
const int ROOT = 0;

static int NNB;
static CUDA_REAL time_send, time_grav, time_out, time_nb;
static long long numInter;
static int icall,ini,isend;
// static int nbodymax;
static int deviceCount;
static int target_size_per_gpu;

static int devid, numGPU;
static bool is_open = false;
static bool devinit = false;
static bool first   = true;
static int J_capacity;
static int I_capacity;

#ifdef unuse // will be depricated soon
extern CUDA_REAL *h_ptcl, *d_ptcl; //, *background;
extern CUDA_REAL *h_result, *d_result;
extern CUDA_REAL *d_r2, *d_diff; //, *d_magnitudes
extern int *d_target;
// extern double3 *d_acc, *d_adot;

CUDA_REAL *h_ptcl=nullptr, *d_ptcl=nullptr;; //, *background;
CUDA_REAL *h_result=nullptr, *d_result=nullptr;
CUDA_REAL *d_r2=nullptr, *d_diff=nullptr; // ,*d_magnitudes=nullptr, 
extern int *h_neighbor, *d_neighbor;
extern int *h_num_neighbor, *d_num_neighbor, *d_neighbor_block;
int *h_neighbor=nullptr, *d_neighbor=nullptr;
int *h_num_neighbor=nullptr, *d_num_neighbor=nullptr;
int *d_neighbor_block=nullptr;
int *d_target=nullptr;
#endif

struct GPU {
    int id = -1;
    cudaStream_t stream = 0;

    // device
    // CUDA_REAL* d_diff = nullptr;
    int*       d_neighbor_block = nullptr;
    int*       d_neighbor  = nullptr;
    int*       d_neighbor_count = nullptr;
    int*       d_neighbor_count_block = nullptr;
    CUDA_REAL* d_result= nullptr;
	Jparticle* dJ = nullptr;
	Iparticle* dI = nullptr;
	//CUDA_INT* dIdx = nullptr;

    // host (should be pinned, but currently not)
    CUDA_REAL* h_result = nullptr;
    int*       h_neighbor_count = nullptr;
    int*       h_neighbor  = nullptr;
};
static std::vector<GPU> gpu;

// double3 *d_adot=nullptr, *d_acc=nullptr;

//#define debuggig_verification
#ifdef debuggig_verification
extern CUDA_REAL *h_r2;
CUDA_REAL *h_r2=nullptr; //only for verification
#endif

#ifdef MultiGPU
cudaStream_t streams[4]; //Maximum 4 GPUs
// cublasHandle_t cublasHandles[4];

#endif
extern CUDA_REAL *h_diff, *h_magnitudes;
CUDA_REAL *h_diff, *h_magnitudes;


/*************************************************************************
 *	 Computing Acceleration
 *************************************************************************/
void GetAcceleration(
    int NumTargetTotal,
    int h_target_list[],
    CUDA_REAL acc[][3],
    CUDA_REAL adot[][3],
    int NumNeighbor[],
    int *NeighborList
) {
    assert(is_open);
	assert((NumTargetTotal > 0) && (NumTargetTotal <= NNB));

    // -----------------------------------------------------------------
    // MULTI-GPU PATH (illustration)
    // -----------------------------------------------------------------
    // Instead of a single handle, each GPU has cublasHandles[i].
    // Also, each GPU has its own d_ptcl_array[i], d_diff_array[i], etc.
	cudaGetDeviceCount(&deviceCount);
	// cublasHandle_t cublasHandles[4];
#ifdef DEBUG
	fprintf(stderr, "Number of GPUs (GetAcceleration): %d\n", deviceCount);
#endif



    // Let’s define chunk = how many targets each GPU will handle at a time:
    // int chunkPerGpu = (variable_size + deviceCount - 1) / deviceCount; 
    int chunkPerGpu = (NNB + deviceCount - 1) / deviceCount; 

    // or any chunk size you prefer
	int NumTarget;

    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += target_size) {
        // How many targets remain in this chunk
        // int bigChunk = std::min(target_size, NumTargetTotal - TargetStart);
		NumTarget = std::min(I_capacity, NumTargetTotal-TargetStart);
#ifdef DEBUG
		fprintf(stderr, "TargetStart = %d, NumTargetTotal = %d, NumTarget = %d\n", TargetStart, NumTargetTotal, NumTarget);
#endif

        // 1) Launch a kernel on each GPU with an offset
        for (int i = 0; i < deviceCount; i++) {
            cudaSetDevice(i);
			dim3 gridDim2(NumTarget, 1);
            dim3 blockDim2(GridDimY, 1);

            int deviceJStart = i * chunkPerGpu;
            // int deviceNumJ   = std::min(chunkPerGpu, variable_size - deviceJStart);
			int deviceNumJ   = std::min(chunkPerGpu, NNB - deviceJStart);

            if (deviceNumJ <= 0) break;  // No more work
#ifdef DEBUG
			fprintf(stderr, "%d, deviceJStart = %d, deviceNumJ = %d\n", i, deviceJStart, deviceNumJ);
#endif

			
            // Prepare kernel dimensions
            dim3 blockDim(BatchSize, 1, 1);
            dim3 gridDim(
                (NumTarget + BatchSize + blockDim.x - 1) / blockDim.x, 
                GridDimY
            );

            // Launch compute_forces on GPU i
			/*
            compute_forces<<<gridDim, blockDim, 0, streams[i]>>>(
                d_ptcl_array[i],
                d_r2_array[i],
                d_diff_array[i],
                NumTarget,
                deviceNumJ, //NNB
                d_target_array[i],
                d_neighbor_block_array[i],
                d_num_neighbor_block_array[i],
                TargetStart,   // i_start
				deviceJStart, // j_start
				NNB
            );
			*/
            compute_forces<<<gridDim, blockDim, 0, streams[i]>>>(
                gpu[i].dI,
                gpu[i].dJ,
                gpu[i].d_result_block,
				gpu[i].d_neighbor_block,
				gpu[i].d_neighbor_count_block,
                NumTarget,
                deviceNumJ, //NNB
                TargetStart,   // i_start
            );
			//cudaStreamSynchronize(streams[i]);

			dim3 blockDim3(16, 6);       // 16 threads along X, 6 along Y
			dim3 gridDim3((NumTarget+15)/16, 1);
			reduce_forces_kernel<<<gridDim3, blockDim3, 0, streams[i]>>>(d_diff_array[i], d_result_array[i], GridDimY, NumTarget);

            gather_neighbor<<<gridDim2, blockDim2, 0, streams[i]>>>(
                d_neighbor_block_array[i], 
                d_num_neighbor_block_array[i], 
                d_neighbor_array[i],
                NumTarget
            );


            gather_numneighbor<<<gridDim2, blockDim2, 0, streams[i]>>>(
                d_num_neighbor_block_array[i], 
                d_num_neighbor_array[i], 
                NumTarget
            );
			//cudaStreamSynchronize(streams[i]);

            // Copy back partial results
            toHost(h_result_array[i],
                   d_result_array[i],
                   _six * NumTarget,
                   streams[i]);


            toHost(NeighborList_array[i],
                   d_neighbor_array[i],
                   NumTarget * MaxNumNeighbor,
                   streams[i]);
			
            toHost(h_num_neighbor_array[i],
                   d_num_neighbor_array[i],
                   NumTarget,
                   streams[i]);
        }

        // 2) Synchronize all devices, then copy results to host
        for (int i = 0; i < deviceCount; i++) {
            cudaSetDevice(i);
			cudaStreamSynchronize(streams[i]);
        }
		
		for (int j = 0; j < NumTarget; j++) {
			int target_idx = TargetStart + j;  // Precompute base index
			int result_idx = _six * target_idx;  // Precompute h_result index
			NumNeighbor[j] = 0;

			// Initialize h_result for this target
			for (int k = 0; k < _six; k++) {
				h_result[result_idx + k] = 0.0;
			}

			// Accumulate results across devices
			for (int i = 0; i < deviceCount; i++) {
				for (int k = 0; k < _six; k++) {
					h_result[result_idx + k] += h_result_array[i][j * _six + k];
				}
				NumNeighbor[j] += h_num_neighbor_array[i][j];
			}
		}

#ifdef NSIGHT
		nvtxRangePushA("NeighborList_array to NeighborList");
#endif

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
#ifdef NSIGHT
		nvtxRangePop();
#endif


#ifdef NSIGHT
		nvtxRangePushA("h_result to acc and adot");
#endif
		/* // modified code by EW 2025.5.24 // not tested yet!!!
		memcpy(acc 	+ TargetStart	, h_result					, NumTarget * 3 * sizeof(CUDA_REAL));
		memcpy(adot + TargetStart	, h_result + NumTarget * 3	, NumTarget * 3 * sizeof(CUDA_REAL));
		*/
		// /* // original code
		for (int i=0; i<NumTarget; i++) {
			acc[i+TargetStart][0]  = h_result[_six*i];
			acc[i+TargetStart][1]  = h_result[_six*i+1];
			acc[i+TargetStart][2]  = h_result[_six*i+2];
			adot[i+TargetStart][0] = h_result[_six*i+3];
			adot[i+TargetStart][1] = h_result[_six*i+4];
			adot[i+TargetStart][2] = h_result[_six*i+5];
			

#ifdef debuggig_verification
			cudaSetDevice(0);
			toHost(h_r2, d_r2_array[0], NNB); // only for verification

			fprintf(stderr, "%d (%d) neighbors of %d = ", i, h_target_list[i], NumNeighbor[i]);
			for (int j=0;j<NumNeighbor[i];j++) {
				fprintf(stderr, "%d, ", NeighborList[i * MaxNumNeighbor + j]);
			}
			fprintf(stderr, "\n");

			// verification
			fprintf(stderr, "%d (%d) neighbors of %d (veri)= ", i, h_target_list[i], NumNeighbor[i]);
			double ix = h_ptcl[h_target_list[i]];
			double iy = h_ptcl[h_target_list[i] + NNB * 1];
			double iz = h_ptcl[h_target_list[i] + NNB * 2];
			double i_r2 = h_r2[h_target_list[i]];
			
			fprintf(stderr, "h_r2 = %e \n", i_r2);

			for (int j=0; j<NNB; j++) {
				double dx = ix - h_ptcl[j];
				double dy = iy - h_ptcl[j + NNB * 1];
				double dz = iz - h_ptcl[j + NNB * 2];
				double r2_temp = dx*dx + dy*dy + dz*dz;
				if (r2_temp < i_r2) {
					fprintf(stderr, "%d, (%e)", j, r2_temp);
				}
			}
			fprintf(stderr, "\n");
			exit(1);
#endif
		}
		// */
#ifdef NSIGHT
		nvtxRangePop();
#endif

    } // end of TargetStart loop
	/*
	for (int i = 0; i < deviceCount; i++){
		cudaSetDevice(i);
		cublasDestroy(cublasHandles[i]);
	}
	*/

}




/*************************************************************************
 *	 Communication with HOST
 *************************************************************************/


void _ReceiveFromHost(
		std::vector<Jparticle>& hJ,
		std::vector<Iparticle>& hI
		){

	//variable_size stands for the (maximum) number of j (background particles)
	//target_size stands for the (maximum) number of i (target particles)
	//in this version of the codes, the two are defined as the same, with higher than NNB

    const size_t nI = hI.size();
    const size_t nJ = hJ.size();
	//nbodymax       = 100000000;
	NNB            = nJ;
	isend++;
	assert(NNB <= nbodymax);
	cudaError_t cudaStatus;


	//printf("CUDA: receive starts\n");
	cudaGetDeviceCount(&deviceCount);

	

	if ((first) || (new_size(NNB) > J_capacity )) {
		//varaiable_size should be the number of j, and target size be the number of i
		J_capacity = new_size(NNB);
		// target_size = ((NNB > nbodymax/NNB) ? int(pow(2,ceil(log(nbodymax/NNB)/log(2.0)))) : NNB);
		I_capacity = J_capacity;
		target_size_per_gpu = I_capacity / deviceCount + 1;
		fprintf(stderr, "variable_size=%d, target_size=%d\n", J_capacity, I_capacity);

		if (!first) {
			for (int i = 0; i < deviceCount; i++){
				cudaSetDevice(i);
				my_free(gpu[i].h_result, gpu[i].d_result);
				my_free(gpu[i].h_neighbor_count, gpu[i].d_neighbor_count);
				my_free(gpu[i].h_neighbor, gpu[i].d_neighbor);
				my_free_d(gpu[i].dJ);
				my_free_d(gpu[i].dI);
				my_free_d(gpu[i].d_neighbor_count_block);
				my_free_d(gpu[i].d_neighbor_block);
			}
		}
		else {
			first = false;
		}

		for (int i = 0; i < deviceCount; i++){
			cudaSetDevice(i);
			my_allocate(&gpu[i].h_result, &gpu[i].d_result, _six*I_capacity); // x,v,m
			my_allocate(&gpu[i].h_neighbor_count, &gpu[i].d_neighbor_count, I_capacity);
			my_allocate(&gpu[i].h_neighbor, &gpu[i].d_neighbor, I_capacity * MaxNumNeighbor);
			my_allocate_d(&gpu[i].d_neighbor_block, GridDimY * NNB_per_block * I_capacity);
			my_allocate_d(&gpu[i].d_neighbor_count_block, GridDimY * I_capacity);
			my_allocate_d(&gpu[i].dJ, J_capacity);
			my_allocate_d(&gpu[i].dI, I_capacity);
		}
		/*
		my_allocate(&h_result       , d_result_array      ,           _six*J_capacity, deviceCount, 0);
		my_allocate(&h_num_neighbor , d_num_neighbor_array, J_capacity, deviceCount, 0);
		my_allocate_d(d_r2_array,        J_capacity, deviceCount, 0);
		my_allocate_d(d_target_array,        J_capacity, deviceCount, 0);
		my_allocate_d(d_diff_array      , _six * GridDimY * I_capacity, deviceCount, 0);
		// C * m / N_device
		// C * (m / N_device)
		my_allocate_d(d_num_neighbor_block_array, GridDimY * I_capacity, deviceCount, 0);
		my_allocate_d(d_neighbor_block_array, GridDimY * NNB_per_block * I_capacity, deviceCount, 0);
		my_allocate_d(d_neighbor_array, MaxNumNeighbor * I_capacity, deviceCount, 0);
		for (int i = 0; i < deviceCount; i++) {
			cudaSetDevice(i);
			cudaMallocHost(&h_result_array[i], _six*J_capacity * sizeof(CUDA_REAL));
			cudaMallocHost(&h_num_neighbor_array[i], J_capacity * sizeof(int));
			cudaMallocHost(&NeighborList_array[i], J_capacity * MaxNumNeighbor * sizeof(int));
		*/
		}
#ifdef debuggig_verification
		cudaMallocHost((void**)&h_r2        ,        variable_size * sizeof(CUDA_REAL)); // only for verification
#endif
	} //end of if (first) || (new_size(NNB) > variable_size)
#ifdef DEBUG
	size_t freeMem, totalMem;
	for (int i = 0; i < deviceCount; i++) {
		cudaSetDevice(i);
		cudaMemGetInfo(&freeMem, &totalMem);
		std::cout << "Device " << i 
				<< " Free memory: " << freeMem / (1024.0 * 1024.0) << " MB, "
				<< "Total memory: " << totalMem / (1024.0 * 1024.0) << " MB" 
				<< std::endl;
	}
#endif
/*
template <typename T>
void toDevice(T *host, T *device, const int size, cudaStream_t &stream) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpyAsync(device, host, size * sizeof(T), cudaMemcpyHostToDevice, stream);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyHostToDevice failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}
*/

    // -------------------- H2D copies --------------------
    // J: partition across devices
    for (int i = 0; i < deviceCount; ++i) {
        cudaSetDevice(i);
        auto [aJ, bJ] = block_range(nJ, i, deviceCount);
        size_t numJ = (bJ > aJ) ? (bJ - aJ) : 0;
        // copy the J slice for this device
        toDevice(hJ.data() + aJ, gpu[i].dJ, numJ, gpu[i].stream);
        // I: replicate to every device (common pattern). If you want to partition I, do block_range on nI instead.
        toDevice(hI.data(), gpu[i].dI, nI, gpu[i].stream);
    }

}


#ifdef MultiGPU
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

	gpu.resize(deviceCount);

	for (int deviceNum = 0; deviceNum < deviceCount; deviceNum++) {
		cudaSetDevice(deviceNum);

		gpu[i].id = i;
		cudaStreamCreate(&gpu[deviceNum].stream);

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

	// Initialize CUDA context
	/*
	cudaError_t cudaStatus = cudaFree(0);
	if (cudaStatus != cudaSuccess) {
		std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return;
	}
	*/

	// CUDA is now initialized and ready to be used
	//std::cout << "CUDA initialized successfully!" << std::endl;

	/*
	if(devinit) return;

	cudaGetDeviceCount(&numGPU);
	assert(numGPU > 0);
	char *gpu_list = getenv("GPU_LIST");
	if(gpu_list)
	{
		numGPU = 0;
		char *p = strtok(gpu_list, " ");
		if (p) {
			devid = atoi(p);
			numGPU++;
		}
		assert(numGPU > 0);
	}else{
		devid=irank%numGPU;
	}
	cudaSetDevice(devid);

#ifdef PROFILE
	//  if(!irank)fprintf(stderr, "***********************\n");
	//  if(!irank)fprintf(stderr, "Initializing NBODY6/GPU library\n");
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");
	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);
	fprintf(stderr, "# GPU initialization - rank: %d; HOST %s; NGPU %d; device: %d %s\n", irank, hostname,numGPU, devid, prop.name);
	//  if(!irank)fprintf(stderr, "***********************\n");
#endif
	devinit = true;
	*/
}
#else //the regacy
void _InitializeDevice(int irank){

	if (MyRank == ROOT) {
	std::cout << "Initializing CUDA ..." << std::endl;
	}
	// Select CUDA device (optional)
	int deviceNum = 0; // Choose GPU device 0
	int deviceCount;
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


	cudaSetDevice(deviceNum);

	// Use CUDA Driver API to get the device associated with the current context
	CUdevice device;
	CUcontext context;
	cuCtxGetCurrent(&context); // Get current CUDA context

	if (context != nullptr) {
		cuCtxGetDevice(&device); // Get the device associated with the current context
		int deviceId=1;
		//cuDeviceGetAttribute(&deviceId, CU_DEVICE_ATTRIBUTE_DEVICE_PARTITIONABLE, device);
		//cuDeviceGetAttribute(&deviceId, CU_DEVICE_ATTRIBUTE_DEVICE_PARTITIONABLE, device);
		//cuDeviceGetAttribute(&deviceId, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY, device);
		//std::cout << "Root processor's current device is: " << deviceId << std::endl;
	} else {
		std::cerr << "Failed to get CUDA context on root processor." << std::endl;
	}

	cudaStreamCreate(&stream);

	if (MyRank == ROOT) {
	std::cout << "There are " << deviceCount << " GPUs." << std::endl;
	}
	if (device < 0 || device >= deviceCount) {
		    // Handle invalid device index
	}

	// Initialize CUDA context
	/*
	cudaError_t cudaStatus = cudaFree(0);
	if (cudaStatus != cudaSuccess) {
		std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return;
	}
	*/

	// CUDA is now initialized and ready to be used
	//std::cout << "CUDA initialized successfully!" << std::endl;

	/*
	if(devinit) return;

	cudaGetDeviceCount(&numGPU);
	assert(numGPU > 0);
	char *gpu_list = getenv("GPU_LIST");
	if(gpu_list)
	{
		numGPU = 0;
		char *p = strtok(gpu_list, " ");
		if (p) {
			devid = atoi(p);
			numGPU++;
		}
		assert(numGPU > 0);
	}else{
		devid=irank%numGPU;
	}
	cudaSetDevice(devid);

#ifdef PROFILE
	//  if(!irank)fprintf(stderr, "***********************\n");
	//  if(!irank)fprintf(stderr, "Initializing NBODY6/GPU library\n");
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");
	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);
	fprintf(stderr, "# GPU initialization - rank: %d; HOST %s; NGPU %d; device: %d %s\n", irank, hostname,numGPU, devid, prop.name);
	//  if(!irank)fprintf(stderr, "***********************\n");
#endif
	devinit = true;
	*/
}
#endif



void _OpenDevice(const int irank){
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
	icall = ini = isend = 0;

	//select GPU========================================//
	_InitializeDevice(irank);

	if(is_open){
		fprintf(stderr, "gpunb: it is already open\n");
		return;
	}
	is_open = true;


#ifdef PROFILE
	//	fprintf(stderr, "RANK: %d ******************\n",irank);
	//	fprintf(stderr, "Opened NBODY6/GPU library\n");
	fprintf(stderr, "# Open GPU regular force - rank: %d\n", irank);
	//fprintf(stderr, "***********************\n");
#endif
}



void _CloseDevice() {
	if(!is_open) {
		fprintf(stderr, "gpunb: it is already close\n");
		return;
	}
	is_open = false;


	cudaError_t error;

	printf("CUDA: ?!! ...\n");
	//my_free(&h_result    , &d_result);
	fprintf(stderr, "result ...\n");
	//my_free(&h_target    , &d_target);
	fprintf(stderr, "target ...\n");
	//my_free(&h_neighbor  , &d_neighbor);
	fprintf(stderr, "neighbor ...\n");
	//my_free(&h_background, &d_background);

	error = cudaGetLastError();
	if (error != cudaSuccess) {
		printf("CUDA error: %s\n", cudaGetErrorString(error));
		// Handle error
	}

#ifdef PROFILE
	fprintf(stderr, "Closed NBODY6/GPU library\n");
	fprintf(stderr, "rank: %d***************\n",devid);
	fprintf(stderr, "time send : %f sec\n", time_send);
	fprintf(stderr, "time grav : %f sec\n", time_grav);
	fprintf(stderr, "time nb   : %f sec\n", time_nb);
	fprintf(stderr, "time out  : %f sec\n", time_out);
	fprintf(stderr, "%f Gflops (gravity part only)\n", 60.e-9 * numInter / time_grav);
	fprintf(stderr, "***********************\n");
#endif
}



void _ProfileDevice(int irank) {
#ifdef PROFILE
	if(icall) {
		fprintf(stderr,"[R.%d-D.%d GPU Reg.F ] Nsend %d  Ngrav %d  <Ni> %d   send(s) %f grav(s) %f  nb(s) %f  out(s) %f  Perf.(Gflops) %f\n",irank,devid,isend,icall,ini/isend,time_send,time_grav,time_nb,time_out,60.e-9*numInter/time_grav);
	}
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
	icall = ini = isend= 0;
#else
	return;
#endif
}


#define mexPrintf printf

inline void gpuMemReport(size_t * avail, size_t * total, 
		        const char * title = 0, const size_t * free = 0, const bool sense = true) 
{
	char tstring[32] = { '\0' };
	cudaMemGetInfo(avail, total);  

	if (free) {
		if (title) {
			strncpy(tstring, title, 31);
		}
		mexPrintf("%s Memory avaliable: Free: %zu, Total: %zu, %s: %zu\n",
				tstring, *avail, *total, (sense) ? "Allocated\0" : "Freed\0", 
				(sense) ? (*free - *avail) : (*avail - *free));
	} else {
		mexPrintf("Memory avaliable: Free: %zu, Total: %zu\n", *avail, *total);  
	}
}



extern "C" {
	void InitializeDevice(int *irank){
		_InitializeDevice(*irank);
	}
	void OpenDevice(const int *irank){
		_OpenDevice(*irank);
	}
	void CloseDevice(){
		_CloseDevice();
	}
	void SendToDevice(std::vector<Jparticle>& hJ, std::vector<Iparticle>& hI) { {
		_ReceiveFromHost(hJ, hI);
	}
	void ProfileDevice(int *irank){
		_ProfileDevice(*irank);
	}
	void CalculateAccelerationOnDevice(int *NumTargetTotal, int *h_target_list, CUDA_REAL acc[][3], CUDA_REAL adot[][3], int NumNeighbor[], int *NeighborList) {
		GetAcceleration(*NumTargetTotal, h_target_list, acc, adot, NumNeighbor, NeighborList);
	}
}

