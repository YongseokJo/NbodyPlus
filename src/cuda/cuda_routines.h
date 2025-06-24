
#ifndef ROUTINES_H
#define ROUTINES_H
#include <cublas_v2.h>
#include <iostream>
#include <cuda_runtime.h>


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

void checkCudaError(cudaError_t result) {
	if (result != cudaSuccess) {
		std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << std::endl;
		exit(EXIT_FAILURE);
	}
}



template <typename T>
void my_allocate(T **host, T **device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(device, size*sizeof(T));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	//*host = (T*)calloc(size, sizeof(T));
	//*host = (T*)malloc(size*sizeof(T));
	/*
		 if (host == NULL) {
		 fprintf(stderr, "Memory allocation failed\n");
		 }
	 */
	//host = new T[size]();
	cudaStatus = cudaMallocHost(host, size*sizeof(T));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMallocHost failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void my_free(T &host, T &device) {
	cudaError_t cudaStatus;
	cudaStatus = cudaFree(device);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	cudaStatus = cudaFreeHost(host);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));
	}
	//free(host);
}

//	For multiple devices
template <typename T>
void my_allocate(T **host, T **device, const int size, int deviceCount, int targetsize) {
    cudaError_t cudaStatus;
	int temp;

    // Allocate device memory for each device
    for (int i = 0; i < deviceCount; i++) {
        cudaStatus = cudaSetDevice(i); // Set the current device
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaSetDevice failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
            continue; // Skip to the next device if setting the device fails
        }
		if (targetsize == 0) {
	        cudaStatus = cudaMalloc(&device[i], size * sizeof(T));
		}
		else{
			temp = size / targetsize;
			cudaStatus = cudaMalloc(&device[i], temp * ((targetsize + deviceCount - 1) / deviceCount) * sizeof(T));
		}
		
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
        }
    }

    // Allocate pinned host memory for all devices
    cudaStatus = cudaMallocHost(host, size * sizeof(T));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMallocHost failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}
template <typename T>
void my_free(T *host, T **device, int deviceCount) {
    cudaError_t cudaStatus;

    // Free device memory for each device
    for (int i = 0; i < deviceCount; i++) {
        cudaStatus = cudaSetDevice(i); // Set the current device
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaSetDevice failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
            continue; // Skip freeing this device if setting the device fails
        }
		
        cudaStatus = cudaFree(device[i]);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaFree failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
        }
    }

    // Free pinned host memory
    cudaStatus = cudaFreeHost(host);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaFreeHost failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}
template <typename T>
void my_allocate_d(T **device, const int size, int deviceCount, int targetsize) {
    cudaError_t cudaStatus;
	int temp;
	
    // Allocate device memory for each device
    for (int i = 0; i < deviceCount; i++) {
        cudaStatus = cudaSetDevice(i); // Set the current device
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaSetDevice failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
            continue; // Skip to the next device if setting the device fails
        }
		if (targetsize == 0) {
	        cudaStatus = cudaMalloc(&device[i], size * sizeof(T));
		}
		else{
			temp = size / targetsize;
			cudaStatus = cudaMalloc(&device[i], temp * ((targetsize + deviceCount - 1) / deviceCount) * sizeof(T));
		}
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMalloc failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
        }
    }
}

template <typename T>
void my_free_d(T &device, int deviceCount) {
    cudaError_t cudaStatus;

    // Free device memory for each device
    for (int i = 0; i < deviceCount; i++) {
        cudaStatus = cudaSetDevice(i); // Set the current device
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaSetDevice failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
            continue; // Skip freeing this device if setting the device fails
        }
		
        cudaStatus = cudaFree(device[i]);
        if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaFree failed for device %d: %s\n", i, cudaGetErrorString(cudaStatus));
        }
    }
}
// ===============================


template <typename T>
void my_allocate_d(T **device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMalloc(device, size*sizeof(T));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void my_free_d(T &device) {
	cudaError_t cudaStatus;
	cudaStatus = cudaFree(device);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}




template <typename T>
void toDevice(T *host, T *device, const int size, cudaStream_t &stream) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpyAsync(device, host, size * sizeof(T), cudaMemcpyHostToDevice, stream);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyHostToDevice failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void toHost(T *host, T *device, const int size, cudaStream_t &stream) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpyAsync(host, device, size * sizeof(T), cudaMemcpyDeviceToHost, stream);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyDeviceToHost failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}


template <typename T>
void toDevice(T *host, T *device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(device, host, size * sizeof(T), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyHostToDevice failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}

template <typename T>
void toHost(T *host, T *device, const int size) {
	cudaError_t cudaStatus;
	cudaStatus = cudaMemcpy(host, device, size * sizeof(T), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpyDeviceToHost failed: %s\n", cudaGetErrorString(cudaStatus));
	}
}



#endif
