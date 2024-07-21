#pragma once
#include "../defs.h"
extern "C" {
	void InitializeDevice(int *irank);
	void OpenDevice(const int *irank);
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(int *_NNB, REAL m[], REAL x[][3], REAL v[][3], REAL r[], REAL mdot[]);
	void CalculateAccelerationOnDevice(int *NumTarget, int *h_target_list, REAL acc[][3], REAL adot[][3], int NumNeighbor[], int **NeighborList);
}
