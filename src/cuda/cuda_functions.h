#pragma once
#include "../def.h"
extern "C" {
	void InitializeDevice(int *irank);
	void OpenDevice(const int *irank);
	void CloseDevice();
	void ProfileDevice(int *irank);
	void SendToDevice(std::vector<Jparticle> &Jparticles, std::vector<Iparticle> &Iparticles);
	void CalculateAccelerationOnDevice(int *NumTarget, std::vector<int>& RegularList);
}
