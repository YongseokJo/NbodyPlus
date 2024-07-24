#include <vector>
#include <iostream>
#include <cmath>
#include "../global.h"



void CalculateIrrgularAcceleration(std::vector<Particle*> ComputationList) {


	REAL dt, mdot, epsilon=1e-6;
	REAL new_time; // 0 for current and 1 for advanced times

	REAL x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	REAL r2, vx; // 0 for current and 1 for predicted values
	REAL a_tmp[Dim], adot_tmp[Dim]; // 0 for current and 1 for predicted accelerations

	REAL m_r3;


	// one more check whether every particle is on the same page
	new_time = ComputationList[0]->CurrentTimeIrr + ComputationList[0]->TimeStepIrr; // the time to be advanced to
	for (Particle* ptcl:ComputationList) {
		if (new_time == ptcl->CurrentTimeIrr + ptcl->TimeStepIrr) { // the time to be advanced to
			std::cout << "Something's wrong with ComputationList!" <<std::endl;
			return;
		}
	}

	// advance particles
	for (Particle* ptcl:ComputationList) {
		dt = ptcl->TimeStepIrr*EnzoTimeStep; // interval of time step
		for (int dim=0; dim<Dim; dim++) {
			ptcl->PredPosition[dim] =\
															 ((ptcl->a_tot[dim][1]*dt/3 + ptcl->a_tot[dim][0])*dt/2 + ptcl->Velocity[dim])*dt + ptcl->Position[dim];
			ptcl->PredVelocity[dim] =\
														 	 (ptcl->a_tot[dim][1]*dt/2 + ptcl->a_tot[dim][0])*dt   + ptcl->Velocity[dim];
		}
	}


	// calculate 0,1th derivative acceleration
	for (Particle* ptcl:ComputationList) {

		r2 = 0.0;
		vx = 0.0;

		x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
		v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

	}


	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
	}

}






