#include "Particle.h"
#include <vector>
#include <iostream>
#include <cmath>


void Particle::predictPosAndVel(double next_time){

	// temporary variables for calculation
	if (PredTime == next_time) {
		return;
	}

	double dt,dt1,dt2;

	dt = next_time - CurrentTimeIrr;
	dt1 = 1.5*dt;
	dt2 = 2.0*dt;

	// only predict the positions if necessary 
	// (if the necessary prediction time differs from value of time already predicted)
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((ForceDot[dim]*dt + Force[dim])*dt + Velocity[dim])*dt + Position[dim];
		PredVelocity[dim] = (ForceDot[dim]*dt1 + Force[dim])*dt2 + Velocity[dim];
	}

	// updated the predicted time
	PredTime = next_time;
}


