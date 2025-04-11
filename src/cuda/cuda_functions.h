#pragma once
#include "../def.h"
extern "C" {
	void InitializeDevice(int *irank);
	void OpenDevice(const int *irank);
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(int *_NNB, CUDA_REAL m[], CUDA_REAL x[][3], CUDA_REAL v[][3], CUDA_REAL r[], CUDA_REAL mdot[]);
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list, CUDA_REAL acc[][3], CUDA_REAL adot[][3], int NumNeighbor[], int *NeighborList);
#ifdef GPU_INITIALIZATION
	void InitializationOnDevice(int *NumTargetTotal, int *h_target_list, 
		CUDA_REAL areg[][3], CUDA_REAL areg_dot[][3], CUDA_REAL airr[][3], CUDA_REAL airr_dot[][3], 
		CUDA_REAL areg_dotdot[][3], CUDA_REAL areg_dotdotdot[][3], CUDA_REAL airr_dotdot[][3], CUDA_REAL airr_dotdotdot[][3], 
		int NumNeighbor[], int *NeighborList);
#endif
}
