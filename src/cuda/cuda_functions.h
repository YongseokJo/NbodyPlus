#pragma once
#include "../def.h"
extern "C" {
	void InitializeDevice(int *irank);
	void OpenDevice(const int *irank);
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(int *_NNB, CUDA_REAL m[], CUDA_REAL x[][3], CUDA_REAL v[][3], CUDA_REAL r[], CUDA_REAL mdot[]);
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list, CUDA_REAL acc[][3], CUDA_REAL adot[][3], int NumNeighbor[], int *NeighborList);

	void RegularWorker(int NumTargetTotal, int Jstart, int Jend, int gpu_id);
	void AllocateDeviceMemory(int N_i,int N_j,int gpu_id);
	void SendToDeviceMPI(int N_i, int N_j, int gpu_id);
}
