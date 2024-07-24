


/*
__global__ void compute_forces_subset(CUDA_REAL* result, const CUDA_REAL* ptcl, const CUDA_REAL *diff, const CUDA_REAL* magnitudes, int n, int m, const int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//if (idx < m * n) {
	int i = subset[idx / n];
	int j = idx % n;
	int six_idx = idx*_six;
	CUDA_REAL scale;
	idx *= _two;
	__shared__ CUDA_REAL res[_six];

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




struct CudaParticle{
  float3 pos;
  float3 vel;
  float  mass;
	float  r2; // for AC neighbor
	float  mdot;

	//BackgroundParticle(int) {}
	CudaParticle(float m, float x[3], float v[3], float _r2, float _mdot){
		mass  = m;
		pos.x = x[0];
    pos.y = x[1];
    pos.z = x[2];
    vel.x = v[0];
    vel.y = v[1];
    vel.z = v[2];
		mdot  = _mdot;
		r2  = _r2;

    NAN_CHECK(x[0]);
    NAN_CHECK(x[1]);
    NAN_CHECK(x[2]);
    NAN_CHECK(m);
    NAN_CHECK(v[0]);
    NAN_CHECK(v[1]);
    NAN_CHECK(v[2]);
    NAN_CHECK(_mdot);
    NAN_CHECK(_r2);
  }

  void setParticle(float m, float x[3], float v[3], float _r2, float _mdot){
		mass  = m;
		pos.x = x[0];
		pos.y = x[1];
		pos.z = x[2];
		vel.x = v[0];
		vel.y = v[1];
		vel.z = v[2];
		mdot  = _mdot;
		r2  = _r2;

		NAN_CHECK(x[0]);
		NAN_CHECK(x[1]);
		NAN_CHECK(x[2]);
		NAN_CHECK(m);
		NAN_CHECK(v[0]);
		NAN_CHECK(v[1]);
		NAN_CHECK(v[2]);
		NAN_CHECK(_mdot);
    NAN_CHECK(_r2);
	}
  //__device__ BackgroundParticle() {}
};



struct Result{
	float3 acc;
	float3 adot;
	//unsigned short num_ac;          //  8 words
	//unsigned short ac_list[MaxNeighbor];

	void clear_h(void) {
		acc.x  = acc.y  = acc.z  = 0.f;
		adot.x = adot.y = adot.z = 0.f;
		//nnb = 0;
	}

	__device__  void clear() {
		acc.x  = acc.y  = acc.z  = 0.f;
		adot.x = adot.y = adot.z = 0.f;
	}

	/*
	__device__ void operator+=(const Result &rhs){
		acc.x  += rhs.acc.x;
		acc.y  += rhs.acc.y;
		acc.z  += rhs.acc.z;
		adot.x += rhs.adot.x;
		adot.y += rhs.adot.y;
		adot.z += rhs.adot.z;
	}
	*/
};

struct  Neighbor{
	int NumNeighbor;
	int NeighborList[NumNeighborMax]; // this needs to be modified.


	__device__ void clear() {
		NumNeighbor = 0;
#pragma unroll
		for (int i=0; i<NumNeighborMax; i++) {
			NeighborList[i] = 0;
		}
	}
	void clear_h() {
		NumNeighbor = 0;
#pragma unroll
		for (int i=0; i<NumNeighborMax; i++) {
			NeighborList[i] = 0;
		}
	}
};


/*
struct  Neighbor_d{
	int NumNeighbor;
	int NeighborList[100]; // this needs to be modified.


	__device__ void clear() {
		NumNeighbor = 0;
#pragma unroll
		for (int i=0; i<100; i++) {
			NeighborList[i] = 0;
		}
	}
	void clear_h() {
		NumNeighbor = 0;
#pragma unroll
		for (int i=0; i<100; i++) {
			NeighborList[i] = 0;
		}
	}
};
*/

/*
struct  Neighbor{
	int width = 2;
	int height = 100;

	int NumNeighbor[2];
	int NeighborList[2][100]; // this needs to be modified.

	int* NeighborList_d;
	int* NumNeighbor_d;
	size_t pitch;

  Neighbor(){
		cudaError_t cudaStatus;
		cudaStatus = cudaMallocPitch(&NeighborList_d, &pitch, width * sizeof(int), height);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMallocPitch failed: %s\n", cudaGetErrorString(cudaStatus));
		}
		cudaStatus = cudaMalloc(&NumNeighbor_d, width * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
		}
  }
	void toHost() {
		cudaError_t cudaStatus;
		cudaStatus = cudaMemcpy(NumNeighbor, NumNeighbor_d, width * sizeof(int), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMemcpy failed: %s\n", cudaGetErrorString(cudaStatus));
		}
		cudaStatus = cudaMemcpy2D(NeighborList, width * sizeof(int), NeighborList_d, pitch, width * sizeof(int), height, cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaMemcpy2d failed: %s\n", cudaGetErrorString(cudaStatus));
		}
	}

	void toDevice() {
		cudaError_t cudaStatus;
		cudaStatus = cudaMemcpy(NumNeighbor_d, NumNeighbor, width * sizeof(int), cudaMemcpyHostToDevice);
		cudaStatus = cudaMemcpy2D(NeighborList_d, pitch, NeighborList, width * sizeof(int), width * sizeof(int), height, cudaMemcpyHostToDevice);
	}

	void free() {
		cudaError_t cudaStatus;
		cudaStatus = cudaFree(NeighborList_d);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaFree list failed: %s\n", cudaGetErrorString(cudaStatus));
		}
		cudaStatus = cudaFree(NumNeighbor_d);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Neighbor cudaFree num failed: %s\n", cudaGetErrorString(cudaStatus));
		}
	}
};
*



