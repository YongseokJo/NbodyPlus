#pragma once
#include "../def.h"
extern "C" {
	void InitializeDevice();
	void OpenDevice();
	void CloseDevice();
	void ProfileDevice();
	void SendToDevice(std::vector<Jparticle> &Jparticles, std::vector<Iparticle> &Iparticles);
	void CalculateAccelerationOnDevice(int *NumTarget, std::vector<int>& RegularList);
}
