#include "Particle.h"
#include "defs.h"
#include <cmath>
#include <algorithm>


double getNewTimeStep(double *F, double dF[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep;

	F2     = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
	Fdot2  = dF[0][0]*dF[0][0] + dF[1][0]*dF[1][0] + dF[2][0]*dF[2][0];
	F2dot2 = dF[0][1]*dF[0][1] + dF[1][1]*dF[1][1] + dF[2][1]*dF[2][1];
	F3dot2 = dF[0][2]*dF[0][2] + dF[1][2]*dF[1][2] + dF[2][2]*dF[2][2];

	TimeStep = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
	TimeStep = std::sqrt(eta*TimeStep);

	return TimeStep;
}

double getBlockTimeStep(double dt) {
	double block_dt;

	block_dt = pow(2, floor(log(dt) / log(2.0)));

	if (block_dt < pow(2,-64)) {
		std::cerr << "The block time step is smaller than the smallest time step of the system. Program terminated.\n";
		exit(EXIT_FAILURE); 
	}

	return block_dt;
}


void Particle::initializeTimeStep() {

	double *FTmp, dtIrr, dtReg;
	int iter;

	FTmp = new double[3];
	for (int i=0; i<Dim; i++) {
		FTmp[i] = FIrr[i] + FReg[i];
	}

	dtIrr = getNewTimeStep(FTmp, dFIrr);
	dtReg = getNewTimeStep(FTmp, dFReg);

	// if chain, half the time step but not implemented as of now

	CurrentTimeIrr = global_time;
	CurrentTimeReg = global_time;

	TimeStepIrr = getBlockTimeStep(dtIrr);
	TimeStepReg = getBlockTimeStep(dtReg);


	if (global_time <= 0) {

		int i = 0;
		while (fmod(global_time, TimeStepIrr) <= 0) {
			TimeStepIrr *= 0.5;
			i++;
			if (i < 16 || TimeStepIrr > pow(2,-40))
				continue;
			TimeStepIrr = pow(2,-40);
			std::cout << "WARNING! TimeStepIrr is too big!\n" << std::endl;
		}

		i=0;
		while (fmod(global_time, TimeStepReg) <= 0) {
			TimeStepReg *= 0.5;
			i++;
			if (i < 16 || TimeStepReg > pow(2,-40))
				continue;
			TimeStepReg = pow(2,-40);
			std::cout << "WARNING! TimeStepReg is too big!\n" << std::endl;
		}
	}

	TimeStepReg = std::min(TimeStepReg, pow(2,-64));
	while (TimeStepReg < TimeStepIrr) 
		TimeStepIrr = 0.5*TimeStepReg;
	
	PredTimeIrr = CurrentTimeIrr + TimeStepIrr;

	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] =  Position[dim];
		Force[dim] = 0.5*Force[dim];
		ForceDot[dim] = ForceDot[dim]/6;

	}
}



