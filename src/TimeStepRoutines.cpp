#include "defs.h"
#include <cmath>
#include <algorithm>
#include <iostream>


double getNewTimeStep(double F[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep;

	Fdot2  = F[0][0]*F[0][0] + F[1][0]*F[1][0] + F[2][0]*F[2][0];
	Fdot2  = F[0][1]*F[0][1] + F[1][1]*F[1][1] + F[2][1]*F[2][1];
	F2dot2 = F[0][2]*F[0][2] + F[1][2]*F[1][2] + F[2][2]*F[2][2];
	F3dot2 = F[0][3]*F[0][3] + F[1][3]*F[1][3] + F[2][3]*F[2][3];

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



