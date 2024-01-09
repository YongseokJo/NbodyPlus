#include "Particle.h"
#include <cmath>


double getNewTimeStep(double *F, double dF[3][4]);
double getBlockTimeStep(double dt);

// Update TimeStepIrr
void Particle::calculateTimeStepIrr(double *f, double df[3][4]) {
	double TimeStepIrrTmp;
	TimeStepIrrTmp = getBlockTimeStep(getNewTimeStep(f, df));
	if (CurrentTimeIrr+TimeStepIrrTmp > CurrentTimeReg+TimeStepReg) {
		TimeStepIrrTmp *= 0.5;
		//fprintf(stdout, "CTS, time irr =%e\n", TimeStepIrr);
		//std::cout << std::flush;
	}
	if (TimeStepIrrTmp > 2*TimeStepIrr) {
		if (fmod(CurrentTimeIrr, 2*TimeStepIrr)==0)
			TimeStepIrrTmp = 2*TimeStepIrr;
		else
			TimeStepIrrTmp = TimeStepIrr;
	} 
	else if (TimeStepIrrTmp < TimeStepIrr) {
		TimeStepIrrTmp = 0.5*TimeStepIrr;
	}
	TimeStepIrr = TimeStepIrrTmp;
}

// Update TimeStepReg
void Particle::calculateTimeStepReg(double *f, double df[3][4]) {
	TimeStepReg = getBlockTimeStep(getNewTimeStep(f, df));
	if (TimeStepIrr > TimeStepReg)
		;
}
