#include "Particle.h"
#include <cmath>


double getNewTimeStep(double *F, double dF[3][4]);
double getBlockTimeStep(double dt);


void Particle::initializeTimeStep() {

	double *FTmp, dtIrr, dtReg;
	int iter;

	FTmp = new double[3];
	for (int i=0; i<Dim; i++) {
		FTmp[i] = FIrr[i] + FReg[i];
	}


	//this->calculateTimeStepReg(FTmp,dFReg);
	//this->calculateTimeStepIrr(FTmp,dFIrr);
	dtIrr = getNewTimeStep(FTmp, dFIrr);
	dtReg = getNewTimeStep(FTmp, dFReg);
	TimeStepIrr = getBlockTimeStep(dtIrr);
	TimeStepReg = getBlockTimeStep(dtReg);

	// if chain, half the time step but not implemented as of now

	CurrentTimeIrr = global_time;
	CurrentTimeReg = global_time;



	if (global_time <= 0) {

		int i = 0;
		while (fmod(global_time, TimeStepIrr) != 0) {
			TimeStepIrr *= 0.5;
			i++;
			if (i < 16 || TimeStepIrr > pow(2,-40))
				continue;
			TimeStepIrr = pow(2,-40);
			std::cout << "WARNING! TimeStepIrr is too big!\n" << std::endl;
		}

		i=0;
		while (fmod(global_time, TimeStepReg) != 0) {
			TimeStepReg *= 0.5;
			i++;
			if (i < 16 || TimeStepReg > pow(2,-40))
				continue;
			TimeStepReg = pow(2,-40);
			std::cout << "WARNING! TimeStepReg is too big!\n" << std::endl;
		}
	}

	TimeStepReg = std::min(TimeStepReg, 1.);
	while (TimeStepReg < TimeStepIrr) 
		TimeStepIrr = 0.5*TimeStepReg;

	PredTimeIrr = CurrentTimeIrr + TimeStepIrr;

}



