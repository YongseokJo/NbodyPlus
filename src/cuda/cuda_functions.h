#pragma once
#include "../def.h"
extern "C" {
	void InitializeDevice(int *irank);
	void OpenDevice(const int *irank);
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(int *_NNB, CUDA_REAL h_ptcl_j[], CUDA_REAL r[]);
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list);
}
