// ============================================================================
// CUDA Kernels for N-body Force Calculation
// ============================================================================

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cassert>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "cuda_defs.h"
#include "../def.h"
#include "cuda_kernels.h"

// ============================================================================
// Debug/utility kernels
// ============================================================================

__global__ void print_forces_subset(cuda_real_t* result, int m, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int m_start = max(0, m - 5);
    if ((idx < m) && (idx >= m_start)) {
        for (int j = 0; j < n; j++) {
            printf("(%d %d) = %e\n", idx, j, result[j * m + idx]);
        }
    }
}

// NTHREAED = 64;
// NJBlock = 28
// NNB_PER_BLOCK = 256;

__global__ void compute_forces(const i_particle_t* __restrict__ d_Ip, const j_particle_t* __restrict__ d_Jp, cuda_real_t* __restrict__ acc, int* __restrict__ neighbor, int* num_neighbor,
	 							int m, int n, int i_start){
	// define i and j. in this code, grid is 2D and block is 1D
    int i = threadIdx.x + blockIdx.x * blockDim.x; // Unique thread index across all blocks
	int tid = threadIdx.x;
	int I;
	// int BATCH_SIZE = blockDim.x;
	int idx_save_size = gridDim.y * m;

	int j_begin = blockIdx.y * n / gridDim.y;
	int j_end = (blockIdx.y + 1) * n / gridDim.y;
	if (blockIdx.y == gridDim.y - 1) j_end = n;  // Ensure the last block covers all remaining elements
	
	while (i < m + BATCH_SIZE){ // even with i > m, the last block needs to assign the shared memory for each tid
		int i_ptcl = (i < m) ? i + i_start : m - 1 + i_start; //assign dummy values for the last block	
		i_particle_t Ip = d_Ip[i_ptcl];
		// cuda_real_t i_r2 = Ip.radius_sq;
		
		int NumNeighbor = 0;
		int idx_save = i * gridDim.y + blockIdx.y;
		// cuda_real_t save_acc[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		cuda_real_t ax=0, ay=0, az=0;
		cuda_real_t jx=0, jy=0, jz=0;
		int* BlockNeighbor = &neighbor[NNB_PER_BLOCK*idx_save]; // Pointer to the neighbor list of the current block

		for (int j=j_begin; j < j_end; j+=BATCH_SIZE){ // total particles
			int current_batch_size = min(BATCH_SIZE, j_end - j);
			// assing shared particles for BATCH_SIZE particles to each block
			__shared__ j_particle_t Jp_sh[BATCH_SIZE];

			__syncthreads();
			if (tid < current_batch_size) {
				Jp_sh[tid] = d_Jp[j + tid];
			}
			__syncthreads();

            #pragma unroll 4
			for (int jj=0; jj<current_batch_size; jj++){
				if (i<m){
					// Calculate forces
					cuda_real_t dx = Jp_sh[jj].pos_x - Ip.pos_x;
					cuda_real_t dy = Jp_sh[jj].pos_y - Ip.pos_y;
					cuda_real_t dz = Jp_sh[jj].pos_z - Ip.pos_z;
					cuda_real_t d2 = dx*dx + dy*dy + dz*dz;

					if (d2 > Ip.radius_sq) {
						// Calculate velocity differences
						cuda_real_t dvx = Jp_sh[jj].vel_x - Ip.vel_x;
						cuda_real_t dvy = Jp_sh[jj].vel_y - Ip.vel_y;
						cuda_real_t dvz = Jp_sh[jj].vel_z - Ip.vel_z;
						cuda_real_t inv_sqrt_d2 = rsqrt(d2);
						cuda_real_t inv_d2 = inv_sqrt_d2 * inv_sqrt_d2; // or 1 / magnitude0
						cuda_real_t scale = Jp_sh[jj].mass * inv_sqrt_d2 * inv_d2;
						// Calculate adot_temp
						cuda_real_t common_factor = 3.0 * (dx*dvx + dy*dvy + dz*dvz) * inv_d2;
						
						ax += scale * dx;
						ay += scale * dy;
						az += scale * dz;
						jx += scale * (dvx - common_factor * dx);
						jy += scale * (dvy - common_factor * dy);
						jz += scale * (dvz - common_factor * dz);
						

					}
					else { //if (i_ptcl != j_start + j + jj) 
						// BlockNeighbor[NumNeighbor++] = j_start + j + jj;
						BlockNeighbor[NumNeighbor++] = Jp_sh[jj].index;
						assert (NumNeighbor < NNB_PER_BLOCK);
					}

				} // end of if (i < m)
			} // end of jj loop
		}// end of j loop
		if (i < m){
			// printf("i, bIdx.y, i_ptcl: (ax, adotx): %d, %d, %d, %.3e, %.3e, %.3e, %d %d %d\n", i, blockIdx.y, i_ptcl, save_acc[2], save_acc[5], pi_x, NumNeighbor, j_begin, j_end);
			acc[idx_save] = ax;
			acc[idx_save + idx_save_size] = ay;
			acc[idx_save + 2 * idx_save_size] = az;
			acc[idx_save + 3 * idx_save_size] = jx;
			acc[idx_save + 4 * idx_save_size] = jy;
			acc[idx_save + 5 * idx_save_size] = jz;
			num_neighbor[idx_save] = NumNeighbor;
		}
		i += gridDim.x * blockDim.x;
	} //end of i loop
}


__global__ void reduce_forces_kernel(const cuda_real_t *diff,  // [6 * m * n] total
                                     cuda_real_t       *result, // [6 * m] output
                                     int n, // "rows" in each component
                                     int m  // "columns"
                                    )
{
    // Each thread handles one (component, column).
    // We'll use a 2D block/grid: x -> column j, y -> component c

    int col = blockIdx.x * blockDim.x + threadIdx.x; // j in [0..m-1]
    int comp = blockIdx.y * blockDim.y + threadIdx.y; // c in [0..5]

    // We only have 6 components total:
    if (comp >= 6 || col >= m) return;

    cuda_real_t sumVal = 0.0;
    // sum over the "row" dimension i in [0..n-1]
    // component c is offset by comp*m*n
    // column j is offset by col*n
    // row i is just + i
    for (int i = 0; i < n; i++){
        sumVal += diff[comp * m * n + (col * n) + i];
    }

    // store interleaved in the output: result[j*6 + comp]
    // matching the cublasDgemv style: incY = 6
    result[col * 6 + comp] = sumVal;
}


// using uint16 type?
__global__ void gather_neighbor(const int* neighbor_block, const int* num_neighbor, int* gathered_neighbor, int m) {
    int i = blockIdx.x;  // index for m (target)
    int b = threadIdx.x; // index for block within m
    int t = threadIdx.y; // index for thread within block
    
    const int num_blocks_per_m = GRID_DIM_Y; // gridDim.y;
    const int num_neighbors_per_block = NNB_PER_BLOCK; // Assuming NNB_PER_BLOCK = blockDim.y
    
    if (i >= m) return;

    // Calculate where in the final gathered array this thread should start writing
    int neighbor_start_index = 0;
    for (int j = 0; j < b; j++) {
        neighbor_start_index += num_neighbor[i * num_blocks_per_m + j];
    }

    int local_neighbor_count = num_neighbor[i * num_blocks_per_m + b];
    assert (neighbor_start_index + local_neighbor_count < MAX_NUM_NEIGHBOR);

    for (int n = 0; n < local_neighbor_count; n++) {
        gathered_neighbor[i * MAX_NUM_NEIGHBOR + neighbor_start_index + n] =
            neighbor_block[(i * num_blocks_per_m + b) * num_neighbors_per_block + n];
    }
}

__global__ void gather_numneighbor(const int* numneighbor_block, int* gathered_numneighbor, int m) {
    int i = threadIdx.x + blockIdx.x * blockDim.x; // Unique thread index across all blocks    
	// NumTarget * GRID_DIM_Y
    if (i >= m) return;

	int temp = 0;
	for (int j = 0; j < GRID_DIM_Y; j++) {
		temp += numneighbor_block[i * GRID_DIM_Y + j];
	}
	gathered_numneighbor[i] = temp;
}

__global__	void initialize(cuda_real_t* result, cuda_real_t* diff, int n, int m, int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < m*n) {
		diff[NUM_FORCE_COMPONENTS*idx    ] = 0.;
		diff[NUM_FORCE_COMPONENTS*idx + 1] = 0.;
		diff[NUM_FORCE_COMPONENTS*idx + 2] = 0.;
		diff[NUM_FORCE_COMPONENTS*idx + 3] = 0.;
		diff[NUM_FORCE_COMPONENTS*idx + 4] = 0.;
		diff[NUM_FORCE_COMPONENTS*idx + 5] = 0.;
	}

#ifdef old
	if (idx < m * n) {
		int i = idx / n;
		int j = idx % n;

		if (j == 0) {
			result[NUM_FORCE_COMPONENTS*i] = 0.;
			result[NUM_FORCE_COMPONENTS*i + 1] = 0.;
			result[NUM_FORCE_COMPONENTS*i + 2] = 0.;
			result[NUM_FORCE_COMPONENTS*i + 3] = 0.;
			result[NUM_FORCE_COMPONENTS*i + 4] = 0.;
			result[NUM_FORCE_COMPONENTS*i + 5] = 0.;
			// num_neighbor[i] = 0;
			/*
			for (j=0; j<MAX_NUM_NEIGHBOR; j++)
				neighbor[MAX_NUM_NEIGHBOR*i+j] = 0;
				*/
		}
	}
#endif
}

// CUDA kernel to compute pairwise differences for a subset of particles
__global__ void compute_pairwise_diff_subset(const cuda_real_t* ptcl, cuda_real_t* diff, int n, int m, const int* subset, int start) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		int i = subset[idx / n + start];
		int j = idx % n;
		idx *= NUM_FORCE_COMPONENTS;
		i *= NUM_RESULT_COMPONENTS;
		j *= NUM_RESULT_COMPONENTS;

		diff[idx]   = ptcl[j]   - ptcl[i];
		diff[idx+1] = ptcl[j+1] - ptcl[i+1];
		diff[idx+2] = ptcl[j+2] - ptcl[i+2];
		diff[idx+3] = ptcl[j+3] - ptcl[i+3];
		diff[idx+4] = ptcl[j+4] - ptcl[i+4];
		diff[idx+5] = ptcl[j+5] - ptcl[i+5];

		//printf("(%d,%d) = %e, %e, %e\n", i/NUM_RESULT_COMPONENTS, j/NUM_RESULT_COMPONENTS,  ptcl[i], ptcl[j], diff[idx]);
	}
}

// n: NNB, m: NumTarget
__global__ void compute_magnitudes_subset(const cuda_real_t *r2, const cuda_real_t* diff, cuda_real_t* magnitudes, int n, int m, int* subset, bool* neighbor2, int start) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < n * m) {
		int i = subset[idx / n + start];
		int j = idx % n;
		int six_idx = NUM_FORCE_COMPONENTS*idx;
		int two_idx = NUM_POS_COMPONENTS*idx;


		magnitudes[two_idx]   += diff[(six_idx)]    *diff[(six_idx)];
		magnitudes[two_idx]   += diff[(six_idx) + 1]*diff[(six_idx) + 1];
		magnitudes[two_idx]   += diff[(six_idx) + 2]*diff[(six_idx) + 2];
		magnitudes[two_idx+1] += diff[(six_idx)]    *diff[(six_idx) + 3];
		magnitudes[two_idx+1] += diff[(six_idx) + 1]*diff[(six_idx) + 4];
		magnitudes[two_idx+1] += diff[(six_idx) + 2]*diff[(six_idx) + 5];

		//printf("(%d,%d) = %e, %e\n", i, j,  magnitudes[two_idx], r2[i]);

		if (magnitudes[two_idx] <= r2[i]) {
			//printf("(%d, %d, %d): %e, %e\n", idx/n, i, j, magnitudes[two_idx], r2[i]);
			magnitudes[two_idx]   = -magnitudes[two_idx];
			neighbor2[idx] = true;
		}
		else {
			neighbor2[idx] = false;
		}
	}
}


// CUDA kernel to compute the forces for a subset of particles
__global__ void compute_forces_subset(const cuda_real_t* ptcl, cuda_real_t *diff, const cuda_real_t* magnitudes, int n, int m, const int* subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m * n) {
		//int i = subset[idx / n];
		int i = idx / n;
		int j = idx % n;
		int six_idx = idx*NUM_FORCE_COMPONENTS;
		idx *= NUM_POS_COMPONENTS;
		cuda_real_t acc[DIM], adot[DIM];

		if (magnitudes[idx] <= 0.) {
			acc[0]  = 0.;
			acc[1]  = 0.;
			acc[2]  = 0.;
			adot[0] = 0.;
			adot[1] = 0.;
			adot[2] = 0.;
		}
		else {
			cuda_real_t scale = ptcl[NUM_RESULT_COMPONENTS*j+6] / (magnitudes[idx] * sqrtf(magnitudes[idx]));
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
__device__ cuda_real_t warpReduce(cuda_real_t val) {
	val += __shfl_down_sync(0xffffffff, val, 16);
	val += __shfl_down_sync(0xffffffff, val, 8);
	val += __shfl_down_sync(0xffffffff, val, 4);
	val += __shfl_down_sync(0xffffffff, val, 2);
	val += __shfl_down_sync(0xffffffff, val, 1);
	return val;
}
*/

__inline__ __device__ cuda_real_t warpReduce(cuda_real_t val)
{
	for (int offset = warpSize/2; offset > 0; offset /= 2) 
		val += __shfl_down_sync(0xffffffff, val, offset);
	return val;
}


#define NEW_FORCE
#ifdef NEW_FORCE
__global__ void reduce_forces(const cuda_real_t *diff, cuda_real_t *result, int n, int m) {

	__shared__ cuda_real_t warpSum[64]; // Assumes max 32 warps per block
	__shared__ cuda_real_t res[NUM_FORCE_COMPONENTS]; //  this is for storing the results
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;
	int bdim = blockDim.x;
	cuda_real_t sum;
	int six_idx;
	int k,l, i = blockIdx.x, j;
	int a = (n+bdim-1)/(bdim);

	if (threadIdx.x < NUM_FORCE_COMPONENTS) 
		res[threadIdx.x] = 0.;
	__syncthreads();


	for (l=0;l<a*bdim;l+=bdim) {
		j = threadIdx.x + l;
		//printf("(%d,%d,%d)\n", i, j, l);
		six_idx = NUM_FORCE_COMPONENTS*(i*n+j);
		warpSum[wid] = 0.;
		__syncthreads();
		#pragma unroll
		for (k=0;k<NUM_FORCE_COMPONENTS;k++) { // ax ay az adotx adoty adotz

			sum = (i < m && j < n) ? diff[six_idx+k] : 0;

			/*
			if (k == 0)
				if (i < m && j < n)
					printf("(%d,%d) = %e, %e\n", i, j, sum, diff[six_idx+k]);
					*/

			// Warp reduce
			sum = warpReduce(sum);
			/*
			if (k == 0)
				if (i < m && j < n)
					printf("first reduction (%d,%d) = %e\n", i, j, sum);
					*/

			// Block reduce
			if (lane == 0) warpSum[wid] = sum;
			__syncthreads();

			if (wid == 0)
			{
				sum = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane] : 0;
				/*
				if (k == 0)
					if (i < m && j < n)
						printf("before second reduction (%d,%d) = %e\n", i, j, sum);
						*/

				sum = warpReduce(sum);

				/*
				if (k == 0)
					if (i < m && j < n)
						printf("second reduction (%d,%d) = %e\n", i, j, sum);
						*/

				if (lane == 0 && i < m) {
					res[k] += sum;
					//printf("%d = (%e,%e)\n", i, sum, res[k]);
				}
			}
			__syncthreads();
		} // reduce across threads
	}
	if (wid == 0 && lane == 0 && i < m) {
		#pragma unroll
		for (k=0; k<NUM_FORCE_COMPONENTS;k++) {
			//printf("%d = (%e)\n", threadIdx.x, res[k]);
			result[NUM_FORCE_COMPONENTS*i+k] = res[k];
		}
	}
	__syncthreads();
}

#else
__global__ void reduce_forces(const cuda_real_t *diff, cuda_real_t *result, int n, int m) {
	int idx = blockIdx.x * n + threadIdx.x;
	__shared__ cuda_real_t warpSum[64]; // Assumes max 32 warps per block
	int lane = threadIdx.x % warpSize;
	int wid = threadIdx.x / warpSize;
	cuda_real_t sum;
	int i = blockIdx.x;
	int j = threadIdx.x;
	int six_idx = NUM_FORCE_COMPONENTS*(i*n+j);
	int k;

	//printf("old version\n");
	#pragma unroll 
	for (k=0;k<NUM_FORCE_COMPONENTS;k++) {
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
			if (lane == 0) result[NUM_FORCE_COMPONENTS*i+k] = sum;
		}
	}
}
#endif




#define MAX_SIZE 9 // maximum size of int array  blockDim.x*MaxSize = total size of int array 
#define NEW_V2 
#ifdef NEW_V2 // this works fine as :)
__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const cuda_real_t* r2, const cuda_real_t* magnitudes, int n, int m, const int *subset) {

	//int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	int bid = blockIdx.x;

		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bid,n,blockDim.x);
	if (tid < n && bid < m) {
		extern __shared__ int sdata[];
		__shared__ int offset;

		sdata[tid+1]=0;
		if (tid==0)  {
			sdata[tid]=0;
			offset = 0;
		}
		__syncthreads();

		int gdim = min(m,gridDim.x);
		int bdim = min(n,blockDim.x);
		int list[MAX_SIZE];
		int n_num;
		int a = (n+bdim*MAX_SIZE-1)/(bdim*MAX_SIZE);
		int start, end;
		int i = subset[bid]; //target particle id
		int j, k;
		int idx = 0;

		//printf("assign_neighbor: %d\n",l);

		// background particles I
		for (k=0; k<a; k++) {

			start = tid+(k*MAX_SIZE*bdim);
			end   = min((k+1)*MAX_SIZE*bdim, n);
			sdata[tid+1] = 0;
			n_num = 0;

			//printf("tid=%d: (start, end) =(%d, %d)\n", tid, start, end);

			// background particles II
			for (j=start; j<end; j+=bdim) {
				if (i != j) {
					//printf("(l,j)=(%d,%d)\n",l,j);
					idx = NUM_POS_COMPONENTS*(n*bid+j);
					//printf("(%d, %d,%d) = %d, %e, %e\n", l, i, j, num_neighbor[l], magnitudes[idx], r2[l]);
					if (magnitudes[idx] < 0) {
						list[n_num] = j;
						n_num++;
						//printf("(%d,%d,%d) = %d, %e, %e\n", l, i, j, n_num, magnitudes[idx], r2[l]);
					}
				}
			} // endfor bk ptcl II
			sdata[tid+1] = n_num;

			__syncthreads();

			if (tid == 0) {

				for (j=2; j<=bdim; j++)
					sdata[j] += sdata[j-1];

				if ((offset+sdata[bdim]) > MAX_NUM_NEIGHBOR) {
					printf("blockid=%d, Too many neighbors (%d, %d)\n", bid, offset, sdata[bdim]);
					assert(offset+sdata[bdim] < MAX_NUM_NEIGHBOR);
				}

				/*
					 if (l == 11) {
					 printf("\n(%d, %d) = sdata[bdim]=%d\n", l, i, sdata[bdim]);
				//for (j=0; j<=bdim; j++) 
				//printf("%d, ",sdata[j]);
				//printf("\n");
				}
				 */
			}
			__syncthreads();

			/*
				 if (l==0)
				 printf("(%d,%d), (num, sdata) =%d, %d\n", l, tid, n_num, sdata[tid]);
			 */

			for (j=0;j<n_num;j++) {
				neighbor[MAX_NUM_NEIGHBOR*bid+offset+sdata[tid]+j] = list[j];
				//printf("(%d,%d), j=%d\n", l, tid, list[j]);
			}
			__syncthreads();

			if (tid == 0) {
				//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
				offset += sdata[bdim];
				//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
			}
			__syncthreads();
		} //endfor bk ptcl I

		if (tid == 0) {
			num_neighbor[bid] = offset; // bid shoud be modified
			offset = 0;
			sdata[0] = 0;
		}
	} // m*n stuff
}

#elif defined(NEW_V1) // works well

__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const cuda_real_t* r2, const cuda_real_t* magnitudes, int n, int m, const int *subset) {

	//int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	int bid = blockIdx.x;

		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bid,n,blockDim.x);
	if (tid < n && bid < m) {
		extern __shared__ int sdata[];
		__shared__ int offset;

		sdata[tid+1]=0;
		if (tid==0)  {
			sdata[tid]=0;
			offset = 0;
		}
		__syncthreads();

		int gdim = min(m,gridDim.x);
		int bdim = min(n,blockDim.x);
		int list[MAX_SIZE];
		int n_num;
		int a = (n+bdim*MAX_SIZE-1)/(bdim*MAX_SIZE);
		int start, end;
		int i; //target particle id
		int j, k, l;
		int idx = 0;


		//printf("%d's bdim=%d, n=%d, blockdim=%d\n", tid,bdim,n,blockDim.x);
		// target particles
		for (l=bid; l<m; l+=gdim) {
			i = subset[l];
			//printf("assign_neighbor: %d\n",l);

			// background particles I
			for (k=0; k<a; k++) {

				start = tid+(k*MAX_SIZE*bdim);
				end   = min((k+1)*MAX_SIZE*bdim, n);
				sdata[tid+1] = 0;
				n_num = 0;

				//printf("tid=%d: (start, end) =(%d, %d)\n", tid, start, end);

				// background particles II
				for (j=start; j<end; j+=bdim) {
					if (i != j) {
						//printf("(l,j)=(%d,%d)\n",l,j);
						idx = NUM_POS_COMPONENTS*(n*l+j);
						//printf("(%d, %d,%d) = %d, %e, %e\n", l, i, j, num_neighbor[l], magnitudes[idx], r2[l]);
						if (magnitudes[idx] < 0) {
							list[n_num] = j;
							n_num++;
							//printf("(%d,%d,%d) = %d, %e, %e\n", l, i, j, n_num, magnitudes[idx], r2[l]);
						}
					}
				} // endfor bk ptcl II
				sdata[tid+1] = n_num;

				__syncthreads();

				if (tid == 0) {

					for (j=2; j<=bdim; j++)
						sdata[j] += sdata[j-1];

					if ((offset+sdata[bdim]) > MAX_NUM_NEIGHBOR) {
						printf("blockid=%d, Too many neighbors (%d, %d)\n", bid, offset, sdata[bdim]);
						assert(offset+sdata[bdim] < MAX_NUM_NEIGHBOR);
					}

					/*
					if (l == 11) {
						printf("\n(%d, %d) = sdata[bdim]=%d\n", l, i, sdata[bdim]);
						//for (j=0; j<=bdim; j++) 
							//printf("%d, ",sdata[j]);
						//printf("\n");
					}
					*/
				}
				__syncthreads();

				/*
					 if (l==0)
					 printf("(%d,%d), (num, sdata) =%d, %d\n", l, tid, n_num, sdata[tid]);
				 */

				for (j=0;j<n_num;j++) {
					neighbor[MAX_NUM_NEIGHBOR*l+offset+sdata[tid]+j] = list[j];
					//printf("(%d,%d), j=%d\n", l, tid, list[j]);
				}
				__syncthreads();

				if (tid == 0) {
					//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
					offset += sdata[bdim];
					//printf("(%d, %d), offset=%d, sdata[bdim]=%d\n", l, i, offset, sdata[bdim]);
				}
				__syncthreads();
			} //endfor bk ptcl I
			if (tid == 0) {
				num_neighbor[l] = offset; // bid shoud be modified
				offset = 0;
				sdata[0] = 0;
			}
			__syncthreads();
		} // endfor target paticles
	} // m*n stuff
}





#else

__global__ void assign_neighbor(int *neighbor, int* num_neighbor, const REAL* r2, const REAL* magnitudes, int n, int m, const int *subset) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < m) {
		int i = subset[idx];
		int k = 0;

		for (int j = 0; j < n; j++) {
			if (i != j) {
				k = NUM_POS_COMPONENTS*(n*idx+j);
				if (magnitudes[k] < 0) {
					//printf("(%d, %d,%d) = %d, %e, %e\n", idx, i, j, num_neighbor[idx], magnitudes[k], r2[i]);
					neighbor[MAX_NUM_NEIGHBOR*idx+num_neighbor[idx]] = j;
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

#endif

