#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include <cassert>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "cuda_types.h"
#include "../defs.h"
//#include "cuda_global.h"
//#include "cuda_functions.h"
//#include "cuda_routines.h"

#define _PROFILE

//#define THREAD 1024 // 2048 for A100
//#define BLOCK 32    // 32 for A100 

//#define THREAD 128 // 2048 for A100
#define THREAD 1 // 2048 for A100
#define BLOCK 1024    // 32 for A100 

#define _six 6
#define _two 2
#define _seven 7

#define ESP2 1e-4
#define new_size(A) ((A > 1024) ? int(pow(2,ceil(log(A)/log(2.0)))) : 1024)


static int NNB;
//static int NumNeighborMax;
static REAL time_send, time_grav, time_out, time_nb;
static long long numInter;
static int icall,ini,isend;
static int nbodymax;
static int devid, numGPU;
static bool is_open = false;
static bool devinit = false;
static bool first   = true;
static int variable_size;
//const int memory_size = 512;
//BackgroundParticle *h_background, *d_background;

extern REAL *h_ptcl, *d_ptcl; //, *background;
extern REAL *h_result, *d_result;
extern REAL *d_diff, *d_magnitudes, *d_r2;
extern int *h_neighbor, *d_neighbor, *h_num_neighbor, *d_num_neighbor;
extern int *d_target;

REAL *h_ptcl=nullptr, *d_ptcl=nullptr;; //, *background;
REAL *h_result=nullptr, *d_result=nullptr;
REAL *d_diff=nullptr,*d_magnitudes=nullptr, *d_r2=nullptr;
int *h_neighbor=nullptr, *d_neighbor=nullptr, *d_num_neighbor=nullptr, *h_num_neighbor=nullptr;
int *d_target=nullptr;


extern cudaStream_t stream;
cudaStream_t stream;

extern REAL *h_diff, *h_magnitudes;
REAL *h_diff, *h_magnitudes;

/*************************************************************************
 *	 Computing Acceleration
 *************************************************************************/
__global__	void initialize(REAL* result, int* neighbor, int* num_neighbor, REAL* diff, REAL *magnitudes, int n, int m, int* subset);
__global__ void compute_pairwise_diff_subset(const REAL* ptcl, REAL* diff, int n, int m, const int* subset);
__global__ void compute_magnitudes_subset(const REAL *r2, const REAL* diff, REAL* magnitudes, int n, int m, int* subset);
__global__ void compute_forces_subset(const REAL* ptcl, REAL *diff, const REAL* magnitudes, int n, int m, const int* subset);
__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const REAL* r2, const REAL* magnitudes, int n, int m, const int *subset);
__global__ void reduce_forces(const REAL *diff, REAL *result, int n, int m);

__device__ void _addition(Result &result, const Result res);
__device__ void _copy(Result &result, const Result res);


// CUDA kernel to compute the forces for a subset of particles
__global__ void print_forces_subset(REAL* result, int m) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m) {
		printf("acc: (%d) = %e\n", idx, result[_six*idx]);
				/*
				atomicAdd(&result[i+1], scale * diff[six_idx + 1]);
				atomicAdd(&result[i+2], scale * diff[six_idx + 2]);

				atomicAdd(&result[i+3], scale * (diff[six_idx + 3] - magnitudes[idx+1]*diff[six_idx    ]/magnitudes[idx]));
				atomicAdd(&result[i+4], scale * (diff[six_idx + 4] - magnitudes[idx+1]*diff[six_idx + 1]/magnitudes[idx]));
				atomicAdd(&result[i+5], scale * (diff[six_idx + 5] - magnitudes[idx+1]*diff[six_idx + 2]/magnitudes[idx]));
				*/
	}
}

void initializeCudaAndCublas(cublasHandle_t* handle);

void GetAcceleration(
		int NumTarget,
		int h_target_list[],
		REAL acc[][3],
		REAL adot[][3],
		int NumNeighbor[],
		int **NeighborList
		) {

	assert(is_open);
	assert((NumTarget > 0) && (NumTarget <= NNB));

	cudaError_t cudaStatus;
	cudaError_t error;

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


	//toDevice(h_target_list, d_target, NumTarget, stream);
	toDevice(h_target_list, d_target, NumTarget);

	// Kernel launch parameters
	dim3 blockSize(variable_size);
	dim3 gridSize(NumTarget);
	//dim3 gridSize((NumTarget * NNB + blockSize.x - 1) / blockSize.x);

	// Compute pairwise differences for the subset
	initialize<<<gridSize, blockSize, 0, stream>>>\
		(d_result, d_neighbor, d_num_neighbor, d_diff, d_magnitudes, NNB, NumTarget, d_target);
	cudaDeviceSynchronize();

	compute_pairwise_diff_subset<<<gridSize, blockSize, 0, stream>>>\
		(d_ptcl, d_diff, NNB, NumTarget, d_target);
	cudaDeviceSynchronize();

	compute_magnitudes_subset<<<gridSize, blockSize, 0, stream>>>\
		(d_r2, d_diff, d_magnitudes, NNB, NumTarget, d_target);
	cudaDeviceSynchronize();

	// Compute gravitational forces for the subset
	compute_forces_subset<<<gridSize, blockSize, 0, stream>>>\
		(d_ptcl, d_diff, d_magnitudes, NNB, NumTarget, d_target);
	assign_neighbor<<<gridSize, blockSize, 0, stream>>>\
		(d_neighbor, d_num_neighbor, d_r2, d_magnitudes, NNB, NumTarget, d_target);
	cudaDeviceSynchronize();

	reduce_forces<<<gridSize, blockSize, 0, stream>>>\
		(d_diff, d_result, NNB, NumTarget);
	cudaDeviceSynchronize();

	/*
	print_forces_subset<<<gridSize, blockSize>>>\
		(d_result, NumTarget);
		*/




	/*
	toHost(h_diff, d_diff, _six*NumTarget*NNB);
	for (int i = 0; i < NumTarget; ++i) {
		//std::cerr << "PID=" << h_target_list[i] << std::endl;
		for (int j = 0; j < NNB; ++j) {
			std::cerr << h_diff[_six*(i * NNB + j)] << " ";
		}
		std::cerr << std::endl;
	}
	*/

	//toHost(h_result  , d_result  , variable_size, stream);
	//toHost(h_neighbor, d_neighbor, variable_size, stream);


	cudaStreamSynchronize(stream); // Wait for all operations to finish

	toHost(h_result      , d_result      ,           _six*NumTarget);
	toHost(h_neighbor    , d_neighbor    , NumNeighborMax*NumTarget);
	toHost(h_num_neighbor, d_num_neighbor,                NumTarget);
	//printf("CUDA: transfer to host done\n");


	//cudaStreamSynchronize(stream); // Wait for all operations to finish


	for (int i=0;i<NumTarget;i++) {
		for (int j=0;j<h_num_neighbor[i];j++) {
			NeighborList[i][j] = h_neighbor[NumNeighborMax*i+j];
			//fprintf(stderr, "%d, ", NeighborList[i][j]);
		}
		//fprintf(stdout, "%d (%d) neighbors of %d", h_target_list[i], i, h_num_neighbor[i]);
		/*
		for (int j=0;j<h_num_neighbor[i];j++) {
			fprintf(stdout, "%d, ", NeighborList[i][j]);
		}
		*/
		//fprintf(stdout, "\n");
		NumNeighbor[i] = h_num_neighbor[i];

		/*
		fprintf(stderr, "PID=%d: a=(%.4e,%.4e,%.4e), adot=(%.4e,%.4e,%.4e)\n",
				h_target_list[i],
				h_result[_six*i],
				h_result[_six*i+1],
				h_result[_six*i+2],
				h_result[_six*i+3],
				h_result[_six*i+4],
				h_result[_six*i+5]
				);
				*/
	}

	// out data
	for (int i=0; i<NumTarget; i++) {
		acc[i][0]  = h_result[_six*i];
		acc[i][1]  = h_result[_six*i+1];
		acc[i][2]  = h_result[_six*i+2];
		adot[i][0] = h_result[_six*i+3];
		adot[i][1] = h_result[_six*i+4];
		adot[i][2] = h_result[_six*i+5];
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


__global__	void initialize(REAL* result, int* neighbor, int* num_neighbor, REAL* diff, REAL *magnitudes, int n, int m, int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	diff[_six*idx    ] = 0.;
	diff[_six*idx + 1] = 0.;
	diff[_six*idx + 2] = 0.;
	diff[_six*idx + 3] = 0.;
	diff[_six*idx + 4] = 0.;
	diff[_six*idx + 5] = 0.;

	magnitudes[_two*idx    ] = 0.;
	magnitudes[_two*idx + 1] = 0.;

	if (idx < m * n) {
		int i = idx / n;
		int j = idx % n;

		if (j == 0) {
			result[_six*i] = 0.;
			result[_six*i + 1] = 0.;
			result[_six*i + 2] = 0.;
			result[_six*i + 3] = 0.;
			result[_six*i + 4] = 0.;
			result[_six*i + 5] = 0.;
			num_neighbor[i] = 0;
			/*
			for (j=0; j<NumNeighborMax; j++)
				neighbor[NumNeighborMax*i+j] = 0;
				*/
		}
	}
}

// CUDA kernel to compute pairwise differences for a subset of particles
__global__ void compute_pairwise_diff_subset(const REAL* ptcl, REAL* diff, int n, int m, const int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		int i = subset[idx / n];
		int j = idx % n;
		idx *= _six;
		i *= _seven;
		j *= _seven;

		diff[idx]   = ptcl[j]   - ptcl[i];
		diff[idx+1] = ptcl[j+1] - ptcl[i+1];
		diff[idx+2] = ptcl[j+2] - ptcl[i+2];
		diff[idx+3] = ptcl[j+3] - ptcl[i+3];
		diff[idx+4] = ptcl[j+4] - ptcl[i+4];
		diff[idx+5] = ptcl[j+5] - ptcl[i+5];

		//printf("(%d,%d) = %e, %e, %e\n", i/_seven, j/_seven,  ptcl[i], ptcl[j], diff[idx]);
	}
}


__global__ void compute_magnitudes_subset(const REAL *r2, const REAL* diff, REAL* magnitudes, int n, int m, int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < n * m) {
		int i = subset[idx / n];
		int j = idx % n;
		int six_idx = _six*idx;
		int two_idx = _two*idx;


		magnitudes[two_idx]   += diff[(six_idx)]    *diff[(six_idx)];
		magnitudes[two_idx]   += diff[(six_idx) + 1]*diff[(six_idx) + 1];
		magnitudes[two_idx]   += diff[(six_idx) + 2]*diff[(six_idx) + 2];
		magnitudes[two_idx+1] += diff[(six_idx)]    *diff[(six_idx) + 3];
		magnitudes[two_idx+1] += diff[(six_idx) + 1]*diff[(six_idx) + 4];
		magnitudes[two_idx+1] += diff[(six_idx) + 2]*diff[(six_idx) + 5];

		//printf("(%d,%d) = %e, %e\n", i, j,  magnitudes[two_idx], r2[i]);

		if (magnitudes[two_idx] <= r2[i]) {
			//printf("(%d, %d): %e, %e\n",subset[i], j, magnitudes[two_idx], r2[i]);
			magnitudes[two_idx]   = -magnitudes[two_idx];
		}
	}
}


// CUDA kernel to compute the forces for a subset of particles
__global__ void compute_forces_subset(const REAL* ptcl, REAL *diff, const REAL* magnitudes, int n, int m, const int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		//int i = subset[idx / n];
		int i = idx / n;
		int j = idx % n;
		int six_idx = idx*_six;
		idx *= _two;
		REAL acc[Dim], adot[Dim];

		if (magnitudes[idx] <= 0.) {
			acc[0]  = 0.;
			acc[1]  = 0.;
			acc[2]  = 0.;
			adot[0] = 0.;
			adot[1] = 0.;
			adot[2] = 0.;
		}
		else {
			REAL scale = ptcl[_seven*j+6] / (magnitudes[idx] * sqrtf(magnitudes[idx]));
			acc[0]  = scale * diff[six_idx];
			acc[1]  = scale * diff[six_idx + 1];
			acc[2]  = scale * diff[six_idx + 2];

			adot[0] = scale * (diff[six_idx + 3] - 3*magnitudes[idx+1]*diff[six_idx    ]/magnitudes[idx]);
			adot[1] = scale * (diff[six_idx + 4] - 3*magnitudes[idx+1]*diff[six_idx + 1]/magnitudes[idx]);
			adot[2] = scale * (diff[six_idx + 5] - 3*magnitudes[idx+1]*diff[six_idx + 2]/magnitudes[idx]);
		}

		diff[six_idx]   = acc[0];
		diff[six_idx+1] = acc[1];
		diff[six_idx+2] = acc[2];
		diff[six_idx+3] = adot[0];
		diff[six_idx+4] = adot[1];
		diff[six_idx+5] = adot[2];

		//printf("compute_forces: (%d, %d) = %e\n", i, j,  diff[six_idx]);
	}
}


/*
__device__ REAL warpReduce(REAL val) {
	val += __shfl_down_sync(0xffffffff, val, 16);
	val += __shfl_down_sync(0xffffffff, val, 8);
	val += __shfl_down_sync(0xffffffff, val, 4);
	val += __shfl_down_sync(0xffffffff, val, 2);
	val += __shfl_down_sync(0xffffffff, val, 1);
	return val;
}
*/

__inline__ __device__ REAL warpReduce(REAL val)
{
	for (int offset = warpSize/2; offset > 0; offset /= 2) 
		val += __shfl_down_sync(0xffffffff, val, offset);
	return val;
}


__global__ void reduce_forces(const REAL *diff, REAL *result, int n, int m) {
	int idx = blockIdx.x * n + threadIdx.x;
	__shared__ REAL warpSum[64]; // Assumes max 32 warps per block
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;
	REAL sum;
	int i = blockIdx.x;
	int j = threadIdx.x;
	int six_idx = _six*(i*n+j);
	int k;
	

#pragma unroll 
	for (k=0;k<_six;k++) {
		sum = (i < m && j < n) ? diff[six_idx+k] : 0;
		/*
		if (k == 0)
			if (i < m && j < n)
				printf("(%d,%d) = %e\n", blockIdx.x,threadIdx.x, diff[six_idx+k]);
				*/

		// Warp reduce
		sum = warpReduce(sum);

		// Block reduce
		if (lane == 0) warpSum[wid] = sum;
		__syncthreads();

		if (wid == 0)
		{
			sum = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane] : 0;
			sum = warpReduce(sum);
			if (lane == 0) result[_six*i+k] = sum;
		}
	}
}



__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const REAL* r2, const REAL* magnitudes, int n, int m, const int *subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m) {
		int i = subset[idx];
		int k = 0;

		for (int j = 0; j < n; j++) {
			if (i != j) {
				k = _two*(n*idx+j);
				if (magnitudes[k] < 0) {
					//printf("(%d, %d,%d) = %d, %e, %e\n", idx, i, j, num_neighbor[idx], magnitudes[k], r2[i]);
					neighbor[NumNeighborMax*idx+num_neighbor[idx]] = j;
					num_neighbor[idx]++;
					if (num_neighbor[idx] > 100)  {
						//printf("Error: (%d, %d,%d) = %d, %e, %e\n", idx, i, j, num_neighbor[idx], magnitudes[k], r2[i]);
						assert(num_neighbor[idx] < 100);
						return;
					}
				}
			}
		}
	}
}




/*
__global__ void compute_forces_subset(REAL* result, const REAL* ptcl, const REAL *diff, const REAL* magnitudes, int n, int m, const int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//if (idx < m * n) {
	int i = subset[idx / n];
	int j = idx % n;
	int six_idx = idx*_six;
	REAL scale;
	idx *= _two;
	__shared__ REAL res[_six];

	if (threadIdx.x == 0) { 
		res[0]=0;
		res[1]=0;
		res[2]=0;
		res[3]=0;
		res[4]=0;
		res[5]=0;
	}

	if (idx >= m * n || i == j || magnitudes[idx] <= 0.) {
		atomicAdd(&res[0], 0.);
		atomicAdd(&res[1], 0.);
		atomicAdd(&res[2], 0.);

		atomicAdd(&res[3], 0.);
		atomicAdd(&res[4], 0.);
		atomicAdd(&res[5], 0.);
	}
	else  {
		scale = ptcl[_seven*j+6] / (magnitudes[idx] *sqrtf(magnitudes[idx]));
		i *= _six;
		atomicAdd(&res[0], scale * diff[six_idx]);
		atomicAdd(&res[1], scale * diff[six_idx + 1]);
		atomicAdd(&res[2], scale * diff[six_idx + 2]);

		atomicAdd(&res[3], scale * (diff[six_idx + 3] - magnitudes[idx+1]*diff[six_idx    ]/magnitudes[idx]));
		atomicAdd(&res[4], scale * (diff[six_idx + 4] - magnitudes[idx+1]*diff[six_idx + 1]/magnitudes[idx]));
		atomicAdd(&res[5], scale * (diff[six_idx + 5] - magnitudes[idx+1]*diff[six_idx + 2]/magnitudes[idx]));
	}
	__syncthreads();

	if (threadIdx.x == 0) { 
		result[i]   = res[0];
		result[i+1] = res[1];
		result[i+2] = res[2];
		result[i+3] = res[3];
		result[i+4] = res[4];
		result[i+5] = res[5];
	}
}
*/


void initializeCudaAndCublas(cublasHandle_t* handle) {
	cudaError_t cudaStat = cudaSetDevice(0);
	if (cudaStat != cudaSuccess) {
		std::cerr << "cudaSetDevice failed!" << std::endl;
		exit(1);
	}

	cublasStatus_t stat = cublasCreate(handle);
	if (stat != CUBLAS_STATUS_SUCCESS) {
		std::cerr << "CUBLAS initialization failed!" << std::endl;
		exit(1);
	}
}


/*************************************************************************
 *	 Communication with HOST
 *************************************************************************/
void _ReceiveFromHost(
		int _NNB,
		REAL m[],
		REAL x[][3],
		REAL v[][3],
		REAL r2[],
		REAL mdot[]
		){
	//time_send -= get_wtime();
	nbodymax       = 100000000;
	NNB            = _NNB;
	//NumNeighborMax = _NumNeighborMax;
	isend++;
	assert(NNB <= nbodymax);
	cudaError_t cudaStatus;

	printf("CUDA: receive starts\n");
	//my_allocate(&h_background, &d_background_tmp, new_size(NNB));
	//cudaMemcpyToSymbol(d_background, &d_background_tmp, new_size(NNB)*sizeof(BackgroundParticle));
	
	if ((first) || (new_size(NNB) > variable_size )) {
		variable_size = new_size(NNB);
		if (!first) {
			my_free(h_ptcl				 , d_ptcl);
			my_free(h_result       , d_result);
			my_free(h_neighbor     , d_neighbor);
			my_free(h_num_neighbor , d_num_neighbor);
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
		my_allocate(&h_num_neighbor , &d_num_neighbor,                variable_size);
		my_allocate(&h_neighbor     , &d_neighbor    , NumNeighborMax*variable_size);
		cudaMalloc((void**)&d_r2        ,        variable_size * sizeof(REAL));
		cudaMalloc((void**)&d_target    ,        variable_size * sizeof(int));
		cudaMalloc((void**)&d_diff      , _six * variable_size * variable_size * sizeof(REAL));
		cudaMalloc((void**)&d_magnitudes, _two * variable_size * variable_size * sizeof(REAL));
		//cudaMallocHost((void**)&h_diff          , _six * variable_size * variable_size * sizeof(REAL));
		//cudaMallocHost((void**)&h_magnitudes    , _two * variable_size * variable_size * sizeof(REAL));
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
	void SendToDevice(int *_NNB, REAL m[], REAL x[][3], REAL v[][3], REAL r2[], REAL mdot[]) {
		_ReceiveFromHost(*_NNB, m, x, v, r2, mdot);
	}
	void ProfileDevice(int *irank){
		_ProfileDevice(*irank);
	}
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list, REAL acc[][3], REAL adot[][3], int NumNeighbor[], int **NeighborList) {
		GetAcceleration(*NumTarget, h_target_list, acc, adot, NumNeighbor, NeighborList);
	}
}

