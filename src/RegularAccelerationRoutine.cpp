#include <iostream>
#include "global.h"

void UpdateNextRegTime(std::vector<Particle*> &particle);
void CalculateRegAccelerationOnGPU(std::vector<int> IndexList, std::vector<Particle*> &particle);

bool RegularAccelerationRoutine(std::vector<Particle*> &particle)
{
	std::cout << "Calculating regular force ..." << std::endl;

	// Calulating regular acceleration of the particles
#define GPU
#ifdef GPU
	if (RegIndexList.size() > 0) {
		for (int i=0; i<3; i++) {
			fprintf(stdout, "Time = %.2f Myr, dtIrr = %.4e Myr, dtReg = %.4e Myr, a_tot = (%.2e,%.2e,%.2e)\n", 
					particle[i]->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
					particle[i]->TimeStepIrr*EnzoTimeStep*1e10/1e6,
					particle[i]->TimeStepReg*EnzoTimeStep*1e10/1e6,
					particle[i]->a_tot[0][0],
					particle[i]->a_tot[1][0],
					particle[i]->a_tot[2][0]
					);
		}
		CalculateRegAccelerationOnGPU(RegIndexList, particle);
	}
#else
	for (Particle *ptcl : particle) {
		if (ptcl->isRegular)
		{
			fprintf(stdout, "Particle ID=%d, Time = %.4e, dtIrr = %.4e, dtReg = %.4e\n",
				 	ptcl->getPID(), ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->TimeStepReg);
			// This only computes accelerations without updating particles.
			ptcl->calculateRegAccelerationSecondOrder(particle);
		}
	}

	for (Particle *ptcl : particle) {
		if (ptcl->isRegular)
		{
			ptcl->calculateRegAccelerationFourthOrder(particle);
		}
	}

	// Update particles
	for (Particle *ptcl : particle) {
		// update the regular time step
		if (ptcl->isRegular) {
			if ((ptcl->NumberOfAC == 0) && (NextRegTime == ptcl->CurrentTimeReg + ptcl->TimeStepReg))
			{
				ptcl->updateParticle(NextRegTime, ptcl->a_tot);
				ptcl->CurrentTimeReg += ptcl->TimeStepReg;
				ptcl->CurrentTimeIrr  = ptcl->CurrentTimeReg;
			}
			else
			{ // Not sure about it
				ptcl->CurrentTimeReg = ptcl->CurrentTimeIrr;
			}

			ptcl->calculateTimeStepReg(ptcl->a_reg, ptcl->a_reg);
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			//std::cout << "NBODY+: time step = " <<  ptcl->TimeStepReg*EnzoTimeStep << std::endl;
			//std::cout << "NBODY+: time step = " <<  ptcl->TimeStepReg*EnzoTimeStep*time_unit << "yr" << std::endl;
			ptcl->isRegular = false;
		}
	}
#endif



	global_time = NextRegTime;
	// update the next regular time step
	UpdateNextRegTime(particle);
	std::cout << "Finishing regular force ...\n" << std::flush;
	return true;
}



void UpdateNextRegTime(std::vector<Particle*> &particle) {

	double time_tmp, time = 1e20;

	for (Particle *ptcl : particle)
	{
		// Next regular time step
		time_tmp = ptcl->CurrentTimeReg + ptcl->TimeStepReg;

		// Find the minum regular time step
		if (time > time_tmp)
			time = time_tmp;
	}

	NextRegTime = std::min(time,1.0);

	// Set isRegular of the particles that will be updated next to 1
	//std::cerr << "Regular: ";
	RegIndexList.clear();
	int i = 0;
	for (Particle* ptcl: particle) {
		time_tmp = ptcl->CurrentTimeReg + ptcl->TimeStepReg;
		if (NextRegTime == time_tmp) {
			//std::cerr << ptcl->getPID() << ' ';
			ptcl->isRegular = true;
			RegIndexList.push_back(i);
		}
		i++;
	}
	//std::cerr << '\n' << std::flush;
}

