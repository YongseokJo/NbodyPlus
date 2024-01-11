#include "Particle.h"
#include <vector>
#include <iostream>
#include <cmath>





/*
 *  Purporse: Predict particle positions and velocities up to second order 
 *            using a_0 and d(a_0)/dt; refer to Nbody6++ manual Eq. 8.3 and 8.4
 *
 *  Date    : 2024.01.09  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *
 */
void Particle::predictParticleSecondOrder(double next_time) {

	// Doubling check
	// temporary variables for calculation
	//if (PredTime == next_time) {
		//return;
	//}

	double dt;

	dt = next_time - CurrentTimeIrr;

	// only predict the positions if necessary 
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((a_tot[1][dim]*dt/6 + a_tot[0][dim])*dt/2 + Velocity[dim])*dt;
		PredVelocity[dim] = (a_tot[1][dim]*dt/2 + a_tot[0][dim])*dt;
	}

	// updated the predicted time
	PredTime = next_time;
}



/*
 *  Purporse: Correct particle positions and velocities up to fourth order 
 *            using a_p and d(a_p)/dt; refer to Nbody6++ manual Eq. 8.9 and 8.10
 *
 *  Date    : 2024.01.10  by Yongseok Jo
 *
 */
void Particle::correctParticleFourthOrder(double next_time) {

	double dt;

	dt = next_time - CurrentTimeIrr;

	// only predict the positions if necessary 
	//for (int dim=0; dim<Dim; dim++) {
		//PredPosition[dim] = ((ForceDot[dim]*dt/6 + Force[dim])*dt/2 + Velocity[dim])*dt + Position[dim];
		//PredVelocity[dim] = (ForceDot[dim]*dt/2 + Force[dim])*dt + Velocity[dim];
	//}

}



void Particle::updateParticle() {


	double dt, dt2, dt3, dt4;

	dt  = TimeStepIrr;
	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;


	// only predict the positions if necessary 
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] += a_tot[dim][1]*dt3/6  + a_tot[dim][0]*dt2/2 + Velocity[dim]*dt;
		Velocity[dim] += a_tot[dim][1]*dt2/2  + a_tot[dim][0]*dt;
		Position[dim] += a_tot[dim][2]*dt4/24 + a_tot[dim][3]*dt4*dt/120;
		Velocity[dim] += a_tot[dim][2]*dt3/6  + a_tot[dim][3]*dt4/24;
	}
}
