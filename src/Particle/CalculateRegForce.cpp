#include "Particle.h"
#include "../defs.h"
#include <vector>
#include <iostream>
#include <cmath>


/*
 *  Purporse: Particle Class
 *
 *  Date    : 2024.01.02  by Yongseok Jo
 *
 */

void direct_sum(double *x, double *v, double r2, double vx,
		double mass, double a[3], double adot[3]) {

	double m_r3;

	r2 += EPS2;  // add softening length
	m_r3 = mass/r2/sqrt(r2);

	for (int dim=0; dim<Dim; dim++){
		a[dim]    += m_r3*x[dim];
		adot[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
	}
}



/*
 *  Purporse: Calculate regular forces
 *
 *  Date    : 2024.01.09  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *
 */
void Particle::calculateRegForce(std::vector<Particle*> &particle, double MinRegTime) {


	int maxNeighborNum = 250;
	double rs2 = 0.1;  // scale radius...? need to check


	double dt, r[2], r2, vx;
	double a0_reg[2][Dim], a0dot_reg[2][Dim], x[Dim], v[Dim];
	double a0_irr[2][Dim], a0dot_irr[2][Dim];
	double treg[2];
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4;

	treg[0] = CurrentTimeReg; // current time of particle
	treg[1] = CurrentTimeReg + TimeStepReg; // the time to be advanced to
	dt = TimeStepReg;
	ACList.clear();
	NumberOfAC = 0;

	// scan over all single particles to find the regular force components
	std::cout <<  "Looping single particles for force...\n" << std::flush;

	for (int i=0; i<2; i++) {
		for (Particle* ptcl: particle) {
			if (ptcl != this)
				continue;

			r2 = 0;
			vx = 0;


			ptcl->predictParticleSecondOrder(treg[i]);
			this->predictParticleSecondOrder(treg[i]);
			for (int dim=0; dim<Dim; dim++) {
				// When particles are not at the current time, extrapolate up to 2nd order
				x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
				v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

				r2 += x[dim]*x[dim];
				vx += v[dim]*x[dim];
			}
			r[i]  = sqrt(r2);

			if (NNB >100 && r[0] < ptcl->RadiusOfAC) {
				NumberOfAC++;
				this->ACList.push_back(ptcl);
				// irregular acceleration
				direct_sum(x ,v, r2, vx, ptcl->Mass, a0_irr[i], a0dot_irr[i]);
			}
			else {
				// regular acceleration
				direct_sum(x ,v, r2, vx, ptcl->Mass, a0_reg[i], a0dot_reg[i]);
			}
			if (i == 0) 
				for (int dim=0; dim<Dim; dim++) {
					a_tot[dim][0] = a0_reg[0][dim] + a0_irr[0][dim];
					a_tot[dim][1] = a0_reg[1][dim] + a0_irr[1][dim];
				}
		} // endfor ptcl
	}// endfor i, 0 for current time, 1 for provisional



	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	// calculated the final corrected forces
	for (int dim=0; dim<Dim; dim++) {
		da_dt2  = (   a0_reg[0][dim] - a0_reg[1][dim]   ) / dt2;
		adot_dt = (a0dot_reg[0][dim] + a0dot_reg[1][dim]) / dt;

		a2 =  12*da_dt2 + 6*adot_dt;
		a3 = (-6*da_dt2 - 2*adot_dt)/dt;
		//Position[dim] += a2*dt4/24 + a3*dt4*dt/120;
		//Velocity[dim] += a2*dt3/6  + a3*dt4/24;

		a_reg[dim][0] = a0_reg[dim][0]; 
		a_reg[dim][1] = a0dot_reg[dim][0]; 
		a_reg[dim][2] = a2; 
		a_reg[dim][3] = a3;
	}

	for (int dim=0; dim<Dim; dim++) {
		da_dt2  = (   a0_irr[0][dim] - a0_irr[1][dim]   ) / dt2;
		adot_dt = (a0dot_irr[0][dim] + a0dot_irr[1][dim]) / dt;

		a2 =  12*da_dt2 + 6*adot_dt;
		a3 = (-6*da_dt2 - 2*adot_dt)/dt;
		//Position[dim] += a2*dt4/24 + a3*dt4*dt/120;
		//Velocity[dim] += a2*dt3/6  + a3*dt4/24;

		a_irr[dim][0] = a0_irr[dim][0]; 
		a_irr[dim][1] = a0dot_irr[dim][0]; 
		a_irr[dim][2] = a2; 
		a_irr[dim][3] = a3;

		a_tot[dim][0] = a_reg[dim][0] + a_irr[dim][0];
		a_tot[dim][1] = a_reg[dim][1] + a_irr[dim][1];
		a_tot[dim][2] = a_reg[dim][2] + a_irr[dim][2];
		a_tot[dim][3] = a_reg[dim][3] + a_irr[dim][3];
	}


	// update the regular time step
	CurrentTimeReg += TimeStepReg;
	if (NumberOfAC == 0) 
		updateParticle(CurrentTimeReg, a_tot);
	this->calculateTimeStepReg(a_reg, a_reg);
	this->isRegular = 0;
}



