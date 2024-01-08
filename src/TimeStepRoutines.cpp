#include "defs.h"
#include <cmath>
#include <algorithm>
#include <iostream>


double getNewTimeStep(double *F, double dF[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep;

	F2     = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
	Fdot2  = dF[0][0]*dF[0][0] + dF[1][0]*dF[1][0] + dF[2][0]*dF[2][0];
	F2dot2 = dF[0][1]*dF[0][1] + dF[1][1]*dF[1][1] + dF[2][1]*dF[2][1];
	F3dot2 = dF[0][2]*dF[0][2] + dF[1][2]*dF[1][2] + dF[2][2]*dF[2][2];

	TimeStep = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
	TimeStep = std::sqrt(eta*TimeStep);
	//std::cout<< TimeStep << " ";
	//exit(EXIT_FAILURE); 
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



