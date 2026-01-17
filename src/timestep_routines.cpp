#include "global.h"
#include "def.h"
#include "particle_data.h"
#include <cmath>
#include <algorithm>
#include <iostream>




double getNewTimeStepReg(double v[3], double df[3][4]) {

	double v2, F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     = df[0][0]*df[0][0] + df[1][0]*df[1][0] + df[2][0]*df[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];
	v2     = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];


	double dta;
  dta  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
	dta  = std::sqrt(eta*dta);
	return dta; //dta;
}


double getNewTimeStepIrr(double f[3][4], double df[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     =   f[0][0]*f[0][0] +  f[1][0]*f[1][0]  + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];


	/*
	else if  (F2 == 0 && F2dot2 != 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
	else if (F2 != 0 && F2dot2 == 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
													 */
	//fprintf(stdout, "NBODY+: f dots: %e, %e, %e, %e\n", F2, Fdot2, F2dot2, F3dot2);
	/*
	if (F2 != 0 && F2dot2 != 0 && dt != 0) {
		DivergentPrevent = std::sqrt(F2/Fdot2)/dt;
		fprintf(stdout, "F2 = %e, F2dot2 = %e, dt = %e\n", F2, F2dot2, dt);
		fprintf(stdout, "DivergentPrevent = %e\n", DivergentPrevent);
	}
	else {
		DivergentPrevent = 1;
	}
	*/

  TimeStep  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
  //TimeStep  = F2/Fdot2;
	//TimeStep  = std::sqrt(DivergentPrevent*eta*TimeStep);
	TimeStep  = std::sqrt(eta*TimeStep);
	//std::cout<< TimeStep << " ";
	//exit(EXIT_FAILURE); 
	return TimeStep;
}

double getNewTimeStep(double f[3][4], double df[3][4]) {

	double F2, Fdot2, F2dot2, F3dot2, TimeStep, DivergentPrevent;

	F2     =   f[0][0]*f[0][0] +  f[1][0]*f[1][0]  + f[2][0]*f[2][0];
	Fdot2  = df[0][1]*df[0][1] + df[1][1]*df[1][1] + df[2][1]*df[2][1];
	F2dot2 = df[0][2]*df[0][2] + df[1][2]*df[1][2] + df[2][2]*df[2][2];
	F3dot2 = df[0][3]*df[0][3] + df[1][3]*df[1][3] + df[2][3]*df[2][3];


	/*
	else if  (F2 == 0 && F2dot2 != 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
	else if (F2 != 0 && F2dot2 == 0)
		DivergentPrevent = 1 ; //std::sqrt(F2/Fdot2/dt);
													 */
	//fprintf(stdout, "NBODY+: f dots: %e, %e, %e, %e\n", F2, Fdot2, F2dot2, F3dot2);
	/*
	if (F2 != 0 && F2dot2 != 0 && dt != 0) {
		DivergentPrevent = std::sqrt(F2/Fdot2)/dt;
		fprintf(stdout, "F2 = %e, F2dot2 = %e, dt = %e\n", F2, F2dot2, dt);
		fprintf(stdout, "DivergentPrevent = %e\n", DivergentPrevent);
	}
	else {
		DivergentPrevent = 1;
	}
	*/

  TimeStep  = (std::sqrt(F2*F2dot2)+Fdot2)/(std::sqrt(Fdot2*F3dot2)+F2dot2);
  //TimeStep  = F2/Fdot2;
	//TimeStep  = std::sqrt(DivergentPrevent*eta*TimeStep);
	TimeStep  = std::sqrt(eta*TimeStep);
	//std::cout<< TimeStep << " ";
	//exit(EXIT_FAILURE); 
	return TimeStep;
}

void getBlockTimeStep(double dt, int &TimeLevel, ull_t &TimeBlock, double &TimeStep) {
	TimeLevel = static_cast<int>(floor(log(dt/enzo_time_step)/log(2.0)));
	//TimeLevel = static_cast<int>(ceil(log(dt/enzo_time_step)/log(2.0)));
	//std::cout << "NBODY+: TimeLevel = " << TimeLevel << std::endl;
	//std::cout << "NBODY+: TimeStep = " << TimeStep << std::endl;

	if (TimeLevel < time_block) {
		//std::cerr << "TimeLevel is less than time block!!" << std::endl;
		TimeLevel = time_block;
	}

	TimeStep = static_cast<double>(pow(2, TimeLevel));
	TimeBlock = static_cast<ull_t>(pow(2, TimeLevel-time_block));
}


// ============================================================================
// SoA-compatible overloads: accept ParticleData& and particle index
// ============================================================================

// Helper: extract velocity as contiguous array from ParticleData
static void extract_velocity(const ParticleData& data, size_t i, double v[3]) {
	v[0] = data.get_vel_x(i);
	v[1] = data.get_vel_y(i);
	v[2] = data.get_vel_z(i);
}

// Helper: extract acceleration array [3][4] from ParticleData for specified type
static void extract_acc_total(const ParticleData& data, size_t i, double a[3][4]) {
	for (int d = 0; d < 3; d++) {
		for (int o = 0; o < 4; o++) {
			a[d][o] = data.get_acc_total(i, d, o);
		}
	}
}

static void extract_acc_reg(const ParticleData& data, size_t i, double a[3][4]) {
	for (int d = 0; d < 3; d++) {
		for (int o = 0; o < 4; o++) {
			a[d][o] = data.get_acc_reg(i, d, o);
		}
	}
}

static void extract_acc_irr(const ParticleData& data, size_t i, double a[3][4]) {
	for (int d = 0; d < 3; d++) {
		for (int o = 0; o < 4; o++) {
			a[d][o] = data.get_acc_irr(i, d, o);
		}
	}
}

// SoA-compatible: getNewTimeStepReg using ParticleData
double getNewTimeStepReg(const ParticleData& data, size_t i) {
	double v[3], df[3][4];
	extract_velocity(data, i, v);
	extract_acc_reg(data, i, df);
	return getNewTimeStepReg(v, df);
}

// SoA-compatible: getNewTimeStepIrr using ParticleData
double getNewTimeStepIrr(const ParticleData& data, size_t i) {
	double f[3][4], df[3][4];
	extract_acc_total(data, i, f);
	extract_acc_irr(data, i, df);
	return getNewTimeStepIrr(f, df);
}

// SoA-compatible: getNewTimeStep using ParticleData
double getNewTimeStep(const ParticleData& data, size_t i) {
	double f[3][4], df[3][4];
	extract_acc_total(data, i, f);
	extract_acc_irr(data, i, df);
	return getNewTimeStep(f, df);
}


