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
#include "../particle.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif



extern int MyRank;
const int ROOT = 0;
extern Particle *particles;
extern int* NewNeighbors;

static int NNB;
static CUDA_REAL time_send, time_grav, time_out, time_nb;
static long long numInter;
static int icall,ini,isend;
// static int nbodymax;
static int deviceCount = 0;

static int devid;
static bool is_open = false;
static bool first   = true;
static int J_capacity;
static int I_capacity;

struct GPU {
    int id = -1;
    cudaStream_t stream = 0;

    // device
    // CUDA_REAL* d_diff = nullptr;
    int*       d_neighbor_block = nullptr;
    int*       d_neighbor  = nullptr;
    int*       d_neighbor_count = nullptr;
    int*       d_neighbor_count_block = nullptr;
    CUDA_REAL* d_result = nullptr;
    CUDA_REAL* d_result_block = nullptr;
	Jparticle* dJ = nullptr;
	Iparticle* dI = nullptr;
	//CUDA_INT* dIdx = nullptr;

    // host (should be pinned, but currently not)
    CUDA_REAL* h_result = nullptr;
    int*       h_neighbor_count = nullptr;
    int*       h_neighbor  = nullptr;
	int        J_start = 0;
	int        J_count = 0;
};
static std::vector<GPU> gpu;


//#define debuggig_verification
#ifdef debuggig_verification
extern CUDA_REAL *h_r2;
CUDA_REAL *h_r2=nullptr; //only for verification
#endif


void block_range(size_t N, int i, int P, size_t &a, size_t &b) {
    a = (N * i) / P;        // floor
    b = (N * (i + 1)) / P;  // floor
    // [a, b) half-open interval
}


/*************************************************************************
 *	 Computing Acceleration
 *************************************************************************/

void GetAcceleration(
    int NumTargetTotal,
	std::vector<int>& RegularList
) {
    assert(is_open);
	assert((NumTargetTotal > 0) && (NumTargetTotal <= NNB));

    // -----------------------------------------------------------------
    // MULTI-GPU PATH (illustration)
    // -----------------------------------------------------------------
    // Instead of a single handle, each GPU has cublasHandles[i].
    // Also, each GPU has its own d_ptcl_array[i], d_diff_array[i], etc.

	// cublasHandle_t cublasHandles[4];
#ifdef DEBUG
	fprintf(stderr, "Number of GPUs (GetAcceleration): %d\n", deviceCount);
#endif



    // Let’s define chunk = how many targets each GPU will handle at a time:
    // or any chunk size you prefer
	int NumTarget;

    for (int TargetStart = 0; TargetStart < NumTargetTotal; TargetStart += I_capacity) {
        // How many targets remain in this chunk
        // int bigChunk = std::min(target_size, NumTargetTotal - TargetStart);
		NumTarget = std::min(I_capacity, NumTargetTotal-TargetStart);
#ifdef DEBUG
		fprintf(stderr, "TargetStart = %d, NumTargetTotal = %d, NumTarget = %d\n", TargetStart, NumTargetTotal, NumTarget);
#endif

        // 1) Launch a kernel on each GPU with an offset
        for (int i = 0; i < deviceCount; i++) {
            cudaSetDevice(i);

			// fprintf(stderr, "GPU %d: J_start = %d, J_count = %d\n", i, gpu[i].J_start, gpu[i].J_count);

			dim3 gridDim2(NumTarget, 1);
            dim3 blockDim2(GridDimY, 1);

            if (gpu[i].J_count <= 0) break;  // No more work

			
            // Prepare kernel dimensions
            dim3 blockDim(BatchSize, 1, 1);
            dim3 gridDim(
            	(NumTarget + BatchSize + blockDim.x - 1) / blockDim.x, 
            	GridDimY
            );
            compute_forces<<<gridDim, blockDim, 0, gpu[i].stream>>>(
                gpu[i].dI,
                gpu[i].dJ,
                gpu[i].d_result_block,
				gpu[i].d_neighbor_block,
				gpu[i].d_neighbor_count_block,
                NumTarget,
                gpu[i].J_count, //NNB
                TargetStart   // i_start
            );
			//cudaStreamSynchronize(gpu[i].stream);

			dim3 blockDim3(16, 6);       // 16 threads along X, 6 along Y
			dim3 gridDim3((NumTarget+15)/16, 1);
			reduce_forces_kernel<<<gridDim3, blockDim3, 0, gpu[i].stream>>>(
				gpu[i].d_result_block, gpu[i].d_result, GridDimY, NumTarget);

            gather_neighbor<<<gridDim2, blockDim2, 0, gpu[i].stream>>>(
                gpu[i].d_neighbor_block, 
                gpu[i].d_neighbor_count_block, 
                gpu[i].d_neighbor,
                NumTarget
            );


            gather_numneighbor<<<gridDim2, blockDim2, 0, gpu[i].stream>>>(
                gpu[i].d_neighbor_count_block, 
                gpu[i].d_neighbor_count, 
                NumTarget
            );
			//cudaStreamSynchronize(streams[i]);

            // Copy back partial results
            toHost(gpu[i].h_result,
                   gpu[i].d_result,
                   _six * NumTarget,
                   gpu[i].stream);


            toHost(gpu[i].h_neighbor,
                   gpu[i].d_neighbor,
                   NumTarget * MaxNumNeighbor,
                   gpu[i].stream);
			
            toHost(gpu[i].h_neighbor_count,
                   gpu[i].d_neighbor_count,
                   NumTarget,
                   gpu[i].stream);
        }

        // 2) Synchronize all devices, then copy results to host
        for (int i = 0; i < deviceCount; i++) {
            cudaSetDevice(i);
			cudaStreamSynchronize(gpu[i].stream);
        }
		
		Particle* ptcl;
		for (int j = 0; j < NumTarget; j++) {

			int target_idx = TargetStart + j;  // Precompute base index
			ptcl = &particles[RegularList[j]];

			ptcl->NewNumberOfNeighbor = 0;
			for (int dim=0; dim < Dim; dim++) {
				ptcl->a_irr[dim][0] = 0.0;
				ptcl->a_irr[dim][1] = 0.0;
			}

			for (int i = 0; i < deviceCount; i++) {
				memcpy(NewNeighbors + ptcl->NeighborsOffset + ptcl->NewNumberOfNeighbor, &gpu[i].h_neighbor[j * MaxNumNeighbor], gpu[i].h_neighbor_count[j] * sizeof(int));
				ptcl->NewNumberOfNeighbor += gpu[i].h_neighbor_count[j];
				for (int k = 0; k < 3; k++) {
					ptcl->a_irr[k][0] += static_cast<double>(gpu[i].h_result[j * _six + k]);
					ptcl->a_irr[k][1] += static_cast<double>(gpu[i].h_result[j * _six + k + 3]);
				}
			}
		}

#ifdef NSIGHT
		nvtxRangePushA("NeighborList_array to NeighborList");
#endif

#ifdef NSIGHT
		nvtxRangePop();
#endif

		// Debugging: print out the first few accelerations
		/*
		for (int i = 0; i < std::min(5, NumTargetTotal); i++) {
			printf("Target %d: acc = (%e, %e, %e), adot = (%e, %e, %e), NumNeighbor = %d\n",
				i, acc[i][0], acc[i][1], acc[i][2],
				adot[i][0], adot[i][1], adot[i][2],
				NumNeighbor[i]);
			for (int j = 0; j < NumNeighbor[i]; j++) {
				printf("  Neighbor %d: %d\n", j, NeighborList[i * MaxNumNeighbor + j]);
			}
		}
		*/

		// 3) Now acc and adot arrays on host have the accumulated results for this chunk
		// You can process them as needed before the next chunk

#ifdef NSIGHT
		nvtxRangePushA("h_result to acc and adot");
#endif

#ifdef NSIGHT
		nvtxRangePop();
#endif

    } // end of TargetStart loop
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
	

	if ((first) || (new_size(NNB) > J_capacity )) {
		//varaiable_size should be the number of j, and target size be the number of i
		J_capacity = new_size(NNB);
		// target_size = ((NNB > nbodymax/NNB) ? int(pow(2,ceil(log(nbodymax/NNB)/log(2.0)))) : NNB);
		I_capacity = J_capacity;

		fprintf(stderr, "variable_size=%d, target_size=%d\n", J_capacity, I_capacity);

		if (!first) {
			for (int i = 0; i < deviceCount; i++) {
				cudaSetDevice(i);
				my_free(gpu[i].h_result, gpu[i].d_result);
				my_free(gpu[i].h_neighbor_count, gpu[i].d_neighbor_count);
				my_free(gpu[i].h_neighbor, gpu[i].d_neighbor);
				my_free_d(gpu[i].dJ);
				my_free_d(gpu[i].dI);
				my_free_d(gpu[i].d_neighbor_count_block);
				my_free_d(gpu[i].d_neighbor_block);
				my_free_d(gpu[i].d_result_block);
			}
		}
		else {
			first = false;
		}
		for (int i = 0; i < deviceCount; i++) {
			cudaSetDevice(i);
			my_allocate(&gpu[i].h_result, &gpu[i].d_result, _six*I_capacity); // x,v,m
			my_allocate(&gpu[i].h_neighbor_count, &gpu[i].d_neighbor_count, I_capacity);
			my_allocate(&gpu[i].h_neighbor, &gpu[i].d_neighbor, I_capacity * MaxNumNeighbor);
			my_allocate_d(&gpu[i].d_neighbor_block, GridDimY * NNB_per_block * I_capacity);
			my_allocate_d(&gpu[i].d_neighbor_count_block, GridDimY * I_capacity);
			my_allocate_d(&gpu[i].dJ, J_capacity);
			my_allocate_d(&gpu[i].dI, I_capacity);
			my_allocate_d(&gpu[i].d_result_block, _six * GridDimY * I_capacity);
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

    // -------------------- H2D copies --------------------
    // J: partition across devices
	size_t aJ, bJ;
	int chunkPerGpu = (NNB + deviceCount - 1) / deviceCount;

    for (int i = 0; i < deviceCount; ++i) {
        cudaSetDevice(i);
		block_range(nJ, i, deviceCount, aJ, bJ);
		size_t numJ = (bJ > aJ) ? (bJ - aJ) : 0; // I think this is redundant! by EW 2025.8.24
        // copy the J slice for this device
        toDevice(hJ.data() + aJ, gpu[i].dJ, numJ, gpu[i].stream);
        // I: replicate to every device (common pattern). If you want to partition I, do block_range on nI instead.
        toDevice(hI.data(), gpu[i].dI, nI, gpu[i].stream);
		gpu[i].J_start = aJ;
		gpu[i].J_count = numJ;
    }
}


void _InitializeDevice(){

	if (MyRank == ROOT) {
		std::cout << "Initializing CUDA ..." << std::endl;
	}
	// Select CUDA device (optional)
	cudaGetDeviceCount(&deviceCount);
	if (MyRank == ROOT) {
        std::cout << "There are " << deviceCount << " GPUs." << std::endl;
    }
	gpu.resize(deviceCount);

	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);

	for (int deviceNum = 0; deviceNum < deviceCount; deviceNum++) {
		cudaSetDevice(deviceNum);

		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, deviceNum);
		fprintf(stdout, "# GPU initialization - MyRank: %d; HOST: %s; NGPU: %d; device: %d %s\n", MyRank, hostname, deviceCount, deviceNum, prop.name);

		gpu[deviceNum].id = deviceNum;
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
            devid = (int)cuDev;
            std::cout << "[Rank " << MyRank << "] Current device from driver context = " << devid << std::endl;
            // Check if devId is valid
            if (devid < 0 || devid >= deviceCount) {
                std::cerr << "Invalid device ID from context: " << devid << std::endl;
            }
        }
    } else {
        std::cerr << "Failed to get CUDA context on root processor. "
                  << "cuCtxGetCurrent returned: " << resCtx << std::endl;
    }
}



void _OpenDevice(){
	time_send = time_grav = time_nb = time_out = 0.0;
	numInter = 0;
	icall = ini = isend = 0;

	//select GPU========================================//
	_InitializeDevice();

	if(is_open){
		fprintf(stderr, "gpunb: it is already open\n");
		return;
	}
	is_open = true;


#ifdef PROFILE
	//	fprintf(stderr, "RANK: %d ******************\n",MyRank);
	//	fprintf(stderr, "Opened NBODY6/GPU library\n");
	fprintf(stderr, "# Open GPU regular force - rank: %d\n", MyRank);
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



void _ProfileDevice() {
#ifdef PROFILE
	if(icall) {
		fprintf(stderr,"[R.%d-D.%d GPU Reg.F ] Nsend %d  Ngrav %d  <Ni> %d   send(s) %f grav(s) %f  nb(s) %f  out(s) %f  Perf.(Gflops) %f\n",MyRank,devid,isend,icall,ini/isend,time_send,time_grav,time_nb,time_out,60.e-9*numInter/time_grav);
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
	void InitializeDevice(){
		_InitializeDevice();
	}
	void OpenDevice(){
		_OpenDevice();
	}
	void CloseDevice(){
		_CloseDevice();
	}
	void SendToDevice(std::vector<Jparticle>& hJ, std::vector<Iparticle>& hI){
		_ReceiveFromHost(hJ, hI);
	}
	void ProfileDevice(){
		_ProfileDevice();
	}
	void CalculateAccelerationOnDevice(int *NumTargetTotal, std::vector<int>& RegularList) {
		GetAcceleration(*NumTargetTotal, RegularList);
	}
}

