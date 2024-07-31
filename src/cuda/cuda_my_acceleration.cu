#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include "../defs.h"
#include "cuda_defs.h"
#include "cuda_kernels.h"
#include "cuda_routines.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif



static int NNB;
static CUDA_REAL time_send, time_grav, time_out, time_nb;
static long long numInter;
static int icall,ini,isend;
static int nbodymax;


static int devid, numGPU;
static bool is_open = false;
static bool devinit = false;
static bool first   = true;
static int variable_size;
static int target_size;

extern CUDA_REAL *h_ptcl, *d_ptcl; //, *background;
extern CUDA_REAL *h_result, *d_result;
extern CUDA_REAL *d_diff, *d_magnitudes, *d_r2;
extern int *d_target;

CUDA_REAL *h_ptcl=nullptr, *d_ptcl=nullptr;; //, *background;
CUDA_REAL *h_result=nullptr, *d_result=nullptr;
CUDA_REAL *d_diff=nullptr,*d_magnitudes=nullptr, *d_r2=nullptr;
int *d_target=nullptr;

#define TEST_CUBLAS
#ifndef TEST_CUBLAS
extern int *h_neighbor, *d_neighbor, *h_num_neighbor, *d_num_neighbor;
int *h_neighbor=nullptr, *d_neighbor=nullptr, *d_num_neighbor=nullptr, *h_num_neighbor=nullptr;

#else
extern bool *h_neighbor, *d_neighbor;
extern int *h_num_neighbor;
bool *h_neighbor=nullptr, *d_neighbor=nullptr;
int *h_num_neighbor=nullptr; // added by wispedia
#endif

extern cudaStream_t stream;
cudaStream_t stream;

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
		int **NeighborList
		) {

	assert(is_open);
	assert((NumTargetTotal > 0) && (NumTargetTotal <= NNB));

	int minGridSize, blockSize, gridSize;
	int sharedMemSize;

	//cudaStreamCreate(&stream);

	cublasHandle_t handle;
	initializeCudaAndCublas(&handle);
	
	/*
	for(int i=0; i<NumTarget; i++) {
		d_result[i].clear();
		d_neighbor[i].clear();
		d_dist = 0.;
	}
	*/
	/*
	fprintf(stderr,"\ntargets=");
	for(int i=0; i<NumTarget; i++) {
		fprintf(stderr,"%d, ", h_target_list[i]);
	}
	fprintf(stderr,"\n");
	*/
	int total_data_num;
	int NumTarget;

	for (int TargetStart=0; TargetStart < NumTargetTotal; TargetStart+=target_size){
		NumTarget = std::min(target_size, NumTargetTotal-TargetStart);
		fprintf(stdout, "TargetStart=%d, NumTargetTotal=%d, NumTarget=%d\n", TargetStart, NumTargetTotal, NumTarget);

		toDevice(h_target_list + TargetStart, d_target, NumTarget, stream);

		// Compute pairwise differences for the subset
		//blockSize = variable_size;
		//gridSize = NumTarget;
		total_data_num = new_size(NNB*NumTarget);
		/******* Initialize *********/
		checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
					initialize, 0, 0));	
		gridSize = (total_data_num + blockSize - 1) / blockSize;

		initialize<<<gridSize, blockSize, 0, stream>>>\
			(d_result, d_diff, d_magnitudes, NNB, NumTarget, d_target);
		cudaDeviceSynchronize();


		/******* Differencese *********/
		checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
					compute_pairwise_diff_subset, 0, 0));	
		gridSize = (total_data_num + blockSize - 1) / blockSize;

		compute_pairwise_diff_subset<<<gridSize, blockSize, 0, stream>>>\
			(d_ptcl, d_diff, NNB, NumTarget, d_target, TargetStart);
		cudaDeviceSynchronize();

		/******* Magnitudes *********/
		checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
					compute_magnitudes_subset, 0, 0));	
		gridSize = (total_data_num + blockSize - 1) / blockSize;

		compute_magnitudes_subset<<<gridSize, blockSize, 0, stream>>>\
			(d_r2, d_diff, d_magnitudes, NNB, NumTarget, d_target, d_neighbor, TargetStart);
		cudaDeviceSynchronize();

		/******* Force *********/
		checkCudaError(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,
					compute_forces_subset, 0, 0));
		gridSize = (total_data_num + blockSize - 1) / blockSize;

		compute_forces_subset<<<gridSize, blockSize, 0, stream>>>\
			(d_ptcl, d_diff, d_magnitudes, NNB, NumTarget, d_target);



		/******* Neighborhood (new) *********/
		// reduce_neighbors(handle, d_neighbor, d_num_neighbor, d_magnitudes, NNB, NumTarget, h_target_list);
		cudaDeviceSynchronize();

		#ifdef NSIGHT
		nvtxRangePushA("Reduction");
		#endif
		/******* Reduction *********/
		reduce_forces_cublas(handle, d_diff, d_result, NNB, NumTarget); //test by wispedia
		//reduce_forces_thrust(d_diff, d_result, NNB, NumTarget);
		cudaDeviceSynchronize();

		#ifdef NSIGHT
		nvtxRangePop();
		#endif


		//print_forces_subset<<<gridSize, blockSize>>>\
			(d_result, NumTarget);	


		cudaStreamSynchronize(stream); // Wait for all operations to finish
		toHost(h_result + _six * TargetStart, d_result, _six * NumTarget);

		#ifdef NSIGHT
		nvtxRangePushA("Neighbor in CPU");
		#endif


		toHost(h_neighbor, d_neighbor, NNB * NumTarget);
		for (int i=0;i<NumTarget;i++) {
			int k = 0;
			int* targetNeighborList = NeighborList[i + TargetStart]; // Cache the row pointer
			int target = h_target_list[i + TargetStart]; // Cache the target value

			for (int j=0;j<NNB;j++) {
				if (h_neighbor[i * NNB + j] && (target != j)) {
					if (k<NumNeighborMax){
						targetNeighborList[k++] = j;
						}
					else {
						fprintf(stderr, "Number of neighbors exceeds the maximum number of neighbors %d\n", k);
						exit(1);
						}
				}
			}
			NumNeighbor[i + TargetStart] = k; // h_num_neighbor[i];

		}
		#ifdef NSIGHT
		nvtxRangePop();
		#endif

	}

	// out data
	for (int i=0; i<NumTargetTotal; i++) {
		acc[i][0]  = h_result[_six*i];
		acc[i][1]  = h_result[_six*i+1];
		acc[i][2]  = h_result[_six*i+2];
		adot[i][0] = h_result[_six*i+3];
		adot[i][1] = h_result[_six*i+4];
		adot[i][2] = h_result[_six*i+5];

		fprintf(stderr, "%d (%d) neighbors of %d = ", i, h_target_list[i], NumNeighbor[i]);
		for (int j=0;j<NumNeighbor[i];j++) {
			fprintf(stderr, "%d, ", NeighborList[i][j]);
		}
		fprintf(stderr, "\n");

	}


	cublasDestroy(handle);
	/*
	my_free(h_background , d_background);
	my_free(h_result     , d_result);
	my_free(h_target     , d_target);
	my_free(h_neighbor   , d_neighbor);
	*/
	//cudaStreamDestroy(stream);
	//my_free_d(do_neighbor);
	//printf("CUDA: done?\n");
}







/*************************************************************************
 *	 Communication with HOST
 *************************************************************************/
void _ReceiveFromHost(
		int _NNB,
		CUDA_REAL m[],
		CUDA_REAL x[][3],
		CUDA_REAL v[][3],
		CUDA_REAL r2[],
		CUDA_REAL mdot[]
		){
	//time_send -= get_wtime();
	nbodymax       = 100000000;
	NNB            = _NNB;
	//NumNeighborMax = _NumNeighborMax;
	isend++;
	assert(NNB <= nbodymax);
	cudaError_t cudaStatus;

	//printf("CUDA: receive starts\n");
	//my_allocate(&h_background, &d_background_tmp, new_size(NNB));
	//cudaMemcpyToSymbol(d_background, &d_background_tmp, new_size(NNB)*sizeof(BackgroundParticle));


	if ((first) || (new_size(NNB) > variable_size )) {
		variable_size = new_size(NNB);
		target_size = ((NNB > nbodymax/NNB) ? int(pow(2,ceil(log(nbodymax/NNB)/log(2.0)))) : NNB);
		fprintf(stderr, "variable_size=%d, target_size=%d\n", variable_size, target_size);

		if (!first) {
			my_free(h_ptcl				 , d_ptcl);
			my_free(h_result       , d_result);
			my_free(h_neighbor     , d_neighbor);
			// my_free(h_num_neighbor , d_num_neighbor);
			cudaFreeHost(h_num_neighbor);
			cudaFree(d_target);
			cudaFree(d_r2);
			cudaFree(d_diff);
			cudaFree(d_magnitudes);

		}
		else {
			first = false;
		}
		my_allocate(&h_ptcl         , &d_ptcl        ,         _seven*variable_size); // x,v,m
		my_allocate(&h_result       , &d_result      ,           _six*variable_size);
		// my_allocate(&h_num_neighbor , &d_num_neighbor,                variable_size);
		// my_allocate(&h_neighbor     , &d_neighbor    , NumNeighborMax*variable_size);
		cudaMalloc((void**)&d_r2        ,        variable_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_target    ,        variable_size * sizeof(int));
		cudaMalloc((void**)&d_diff      , _six * variable_size * target_size * sizeof(CUDA_REAL));
		cudaMalloc((void**)&d_magnitudes, _two * variable_size * target_size * sizeof(CUDA_REAL));
		//cudaMallocHost((void**)&h_diff          , _six * variable_size * variable_size * sizeof(CUDA_REAL));
		//cudaMallocHost((void**)&h_magnitudes    , _two * variable_size * variable_size * sizeof(CUDA_REAL));
		#ifdef TEST_CUBLAS
		my_allocate(&h_neighbor     , &d_neighbor    , variable_size * target_size);
		cudaMallocHost((void**)&h_num_neighbor, variable_size * sizeof(int));
		#endif
		
	}


	for (int j=0; j<NNB; j++) {
		for (int dim=0; dim<Dim; dim++) {
			h_ptcl[_seven*j+dim]   = x[j][dim];
			h_ptcl[_seven*j+dim+3] = v[j][dim];
		}
		h_ptcl[_seven*j+6] = m[j];
		//h_particle[j].setParticle(m[j], x[j], v[j], r2[j], mdot[j]);
	}

	//toDevice(h_background,d_background,variable_size);
	toDevice(h_ptcl,d_ptcl, _seven*NNB, stream);
	toDevice(r2    ,d_r2  ,        NNB, stream);
	//fprintf(stdout, "CUDA: receive done\n");
}



void _InitializeDevice(int irank){

	std::cout << "Initializing CUDA ..." << std::endl;
	// Select CUDA device (optional)
	int device = 0; // Choose GPU device 0
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);

	cudaStreamCreate(&stream);

	std::cout << "There are " << deviceCount << " GPUs." << std::endl;
	if (device < 0 || device >= deviceCount) {
		    // Handle invalid device index
	}

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");


	char hostname[150];
	memset(hostname,0,150);
	gethostname(hostname,150);
	

	fprintf(stderr, "# GPU initialization - rank: %d; HOST %s; NGPU %d; device: %d %s\n", irank, hostname,numGPU, devid, prop.name);


	cudaSetDevice(device);

	// Initialize CUDA context
	/*
	cudaError_t cudaStatus = cudaFree(0);
	if (cudaStatus != cudaSuccess) {
		std::cerr << "CUDA initialization failed: " << cudaGetErrorString(cudaStatus) << std::endl;
		return;
	}
	*/

	is_open = true;
	// CUDA is now initialized and ready to be used
	std::cout << "CUDA initialized successfully!" << std::endl;

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
	void SendToDevice(int *_NNB, CUDA_REAL m[], CUDA_REAL x[][3], CUDA_REAL v[][3], CUDA_REAL r2[], CUDA_REAL mdot[]) {
		_ReceiveFromHost(*_NNB, m, x, v, r2, mdot);
	}
	void ProfileDevice(int *irank){
		_ProfileDevice(*irank);
	}
	void CalculateAccelerationOnDevice(int *NumTargetTotal, int *h_target_list, CUDA_REAL acc[][3], CUDA_REAL adot[][3], int NumNeighbor[], int **NeighborList) {
		GetAcceleration(*NumTargetTotal, h_target_list, acc, adot, NumNeighbor, NeighborList);
	}
}

