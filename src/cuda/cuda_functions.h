#pragma once
#include "../def.h"
extern "C" {
	void InitializeDevice();
	void OpenDevice();
	void CloseDevice();
	void ProfileDevice();
	void SendToDevice(std::vector<j_particle_t> &j_particle_ts, std::vector<i_particle_t> &i_particle_ts);
	void CalculateAccelerationOnDevice(int *NumTarget, std::vector<int>& RegularList);
}
