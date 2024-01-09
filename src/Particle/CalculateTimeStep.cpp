#include "Particle.h"
#include <cmath>


double getNewTimeStep(double *F, double dF[3][4]);
double getBlockTimeStep(double dt);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(double *f, double df[3][4]) {
	TimeStepIrr = getBlockTimeStep(getNewTimeStep(f, df));
	if (CurrentTimeIrr+TimeStepIrr > CurrentTimeReg+TimeStepReg) {
		TimeStepIrr *= 0.5;
		//fprintf(stdout, "CTS, time irr =%e\n", TimeStepIrr);
		//std::cout << std::flush;
	}
}

// Update TimeStepReg
void Particle::calculateTimeStepReg(double *f, double df[3][4]) {
	TimeStepReg = getBlockTimeStep(getNewTimeStep(f, df));
	if (TimeStepIrr > TimeStepReg)
		;
}
