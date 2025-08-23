#pragma once
#include "../def.h"
extern "C" {
	void InitializeDevice(int *irank);
	void OpenDevice(const int *irank);
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(std::vector<Jparticle> &Jparticles, std::vector<Iparticle> &Iparticles);
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list, CUDA_REAL acc[][3], CUDA_REAL adot[][3], int NumNeighbor[], int *NeighborList, std::vector<int>& RegularListIndices);
}
