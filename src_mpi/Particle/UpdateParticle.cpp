#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include "../global.h"
#include "../def.h"






void Particle::predictParticleSecondOrder(double dt, CUDA_REAL pos[], CUDA_REAL vel[]) {
	// Doubling check
	// temporary variables for calculation

	// only predict the positions if necessary
	// how about using polynomial correction here?
	
	dt = dt*EnzoTimeStep;

	if (dt == 0) {
		for (int dim=0; dim<Dim; dim++) {
			pos[dim] = (CUDA_REAL)Position[dim];
			vel[dim] = (CUDA_REAL)Velocity[dim];
		}
	}
	else {
		for (int dim=0; dim<Dim; dim++) {
			pos[dim] = (CUDA_REAL) ((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim];
			vel[dim] = (CUDA_REAL) (a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim];
		}
	}
	return;
}

void Particle::predictParticleSecondOrder(double dt, double pos[], double vel[]) {
	// Doubling check
	// temporary variables for calculation

	// only predict the positions if necessary
	// how about using polynomial correction here?
	
	dt = dt*EnzoTimeStep;

	if (dt == 0) {
		for (int dim=0; dim<Dim; dim++) {
			pos[dim] = Position[dim];
			vel[dim] = Velocity[dim];
		}
	}
	else {
		for (int dim=0; dim<Dim; dim++) {
			pos[dim] = ((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim];
			vel[dim] =  (a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim];
		}
	}
	return;
}


/*
 *  Purporse: Correct particle positions and velocities up to fourth order
 *            using a_p and d(a_p)/dt; refer to Nbody6++ manual Eq. 8.9 and 8.10
 *
 *  Date    : 2024.01.10  by Yongseok Jo
 *  Modified: 2024.01.11  by Seoyoung Kim
 *
 */
void Particle::correctParticleFourthOrder(double dt, double pos[], double vel[], double a[3][4]) {
	double dt3,dt4,dt5;

	dt = dt*EnzoTimeStep;

	dt3 = dt*dt*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	// correct the predicted values positions and velocities at next_time
	// and save the corrected values to particle positions and velocities
	// the latest values of a2dot and a3dots (obtained from hermite method) are used
	for (int dim=0; dim<Dim; dim++) {
		NewPosition[dim] = pos[dim]+ a[dim][2]*dt4/24 + a[dim][3]*dt5/120;
		NewVelocity[dim] = vel[dim]+ a[dim][2]*dt3/6  + a[dim][3]*dt4/24;
	}
}


/*
void Particle::polynomialPrediction(double current_time) {

}
*/


void Particle::updateParticle() {
	
	for (int dim=0; dim<Dim; dim++) {
		this->Position[dim] = this->NewPosition[dim];
		this->Velocity[dim] = this->NewVelocity[dim];
	}
	
	//updateTimeStep();
}



void Particle::updateRadius() {

	/* exponential (aggressive) */
	/*
		 const double c = 0.5;
		 const double b = std::log(2) / (NumNeighborMax);  // ln(2) / 40
		 double exp = a * (std::exp(b * NumberOfAC) - 1);
		 */

	/* n=2 polynomial (mild) as n increases it grows mild */

	if (this->NumberOfNeighbor > FixNumNeighbor) {
		const int n = 2;
		const double c = (NumNeighborMax-FixNumNeighbor);
		const double b = 0.9 / std::pow(c,n);  // ln(2) / 40
		double x = this->NumberOfNeighbor-FixNumNeighbor;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		this->RadiusOfNeighbor *= (1.-a);
	}
	else if (this->NumberOfNeighbor < FixNumNeighbor) {
		const int n = 3;
		const double c = (NumNeighborMax-FixNumNeighbor);
		const double b = 0.5 / std::pow(c,n);  // ln(2) / 40
		double x = this->NumberOfNeighbor-FixNumNeighbor;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		this->RadiusOfNeighbor *= (1.-a);
	}
}







double getNewTimeStepReg(double v[3], double f[3][4]);
double getNewTimeStepIrr(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, ULL &TimeBlock, double &TimeStep);

void Particle::calculateTimeStepIrr() {
	double TimeStepTmp;
	int TimeLevelTmp, TimeLevelTmp0;
	ULL TimeBlockTmp;

	if (this->NumberOfNeighbor == 0) {
		TimeLevelIrr = TimeLevelReg;
		TimeStepIrr = static_cast<double>(pow(2, TimeLevelReg));
		TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelReg-time_block));
		return;
	}

	getBlockTimeStep(getNewTimeStepIrr(a_tot, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	/*
	if (NextRegTimeBlock < TimeBlockReg && isRegular) {
		fprintf(stderr, "PID=%d, NextRegTimeBlock=%llu, TimeBlockReg=%llu\n", PID, NextRegTimeBlock, TimeBlockReg);
	}
	*/

	if (TimeLevelTmp > TimeLevelIrr) {
		if (fmod(NewCurrentBlockIrr, 2*TimeBlockIrr)==0) {
			TimeLevelTmp = TimeLevelIrr+1;
			TimeBlockTmp = 2*TimeBlockIrr;
		}
		else {
			TimeLevelTmp = TimeLevelIrr;
			TimeBlockTmp = TimeBlockIrr;
		}
	}
	else if (TimeLevelTmp < TimeLevelIrr) {
		if (TimeLevelTmp < TimeLevelIrr-1) {
			TimeLevelTmp = TimeLevelIrr - 2;
			TimeBlockTmp = TimeBlockIrr/4;
		}
		else {
			TimeLevelTmp = TimeLevelIrr - 1;
			TimeBlockTmp = TimeBlockIrr/2;
		}
	} else {
		TimeLevelTmp = TimeLevelIrr;
		TimeBlockTmp = TimeBlockIrr;
	}


	while (((NewCurrentBlockIrr < CurrentBlockReg+TimeBlockReg) && (NewCurrentBlockIrr+TimeBlockTmp > CurrentBlockReg+TimeBlockReg)) || (TimeLevelTmp >= TimeLevelReg)) {
		/*
		fprintf(stderr,"CurrentBlockIrr = %llu\n",
				CurrentBlockIrr);
		fprintf(stderr,"CurrentBlockReg = %llu\n",
				CurrentBlockReg);
		fprintf(stderr,"TimeBlockReg    = %llu\n",
				TimeBlockReg);
		fprintf(stderr,"TimeBlockIrr    = %llu\n",
				TimeBlockIrr);
				*/
		//if (TimeLevelTmp > TimeLevelReg)
			//fprintf(stderr, "PID=%d, Irr=%d, Reg=%d\n", PID, TimeLevelTmp0, TimeLevelReg);

		TimeLevelTmp--;
		TimeBlockTmp *= 0.5;
	}

	TimeLevelIrr = TimeLevelTmp;

	if (TimeLevelIrr < time_block) {
		//std::cerr << "TimeLevelIrr is too small" << std::endl;
		TimeLevelIrr = std::max(time_block, TimeLevelIrr);
	}


	TimeStepIrr = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));

	if (TimeStepIrr*EnzoTimeStep*1e4<1e-11) {
		fprintf(stderr, "Too small TimeStepIrr! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				PID, TimeStepIrr*EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*EnzoTimeStep*1e4);
		while (TimeStepIrr*EnzoTimeStep*1e4<1e-11) {
			TimeLevelIrr++;
			TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
			TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));
		}
	}

	if (TimeStepIrr > 1) {
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}


	/*
	if (TimeStepIrr*EnzoTimeStep*1e4 < KSTime && (this->isCMptcl == false))
		BinaryCandidateList.push_back(this);
	if (PID == 430) {
		std::cerr << "After TimeLevelIrr=" << TimeLevelIrr << std::endl;
		std::cerr << std::endl;
	}
	*/
}

void Particle::calculateTimeStepIrr2() {
	double TimeStepTmp;
	int TimeLevelTmp, TimeLevelTmp0;
	ULL TimeBlockTmp;

	if (this->NumberOfNeighbor == 0) {
		TimeLevelIrr = TimeLevelReg;
		TimeStepIrr = static_cast<double>(pow(2, TimeLevelReg));
		TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelReg-time_block));
		return;
	}

	getBlockTimeStep(getNewTimeStepIrr(a_tot, a_irr), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	/*
	if (NextRegTimeBlock < TimeBlockReg && isRegular) {
		fprintf(stderr, "PID=%d, NextRegTimeBlock=%llu, TimeBlockReg=%llu\n", PID, NextRegTimeBlock, TimeBlockReg);
	}
	*/


	while (((NewCurrentBlockIrr < CurrentBlockReg+TimeBlockReg) && (NewCurrentBlockIrr+TimeBlockTmp > CurrentBlockReg+TimeBlockReg)) || (TimeLevelTmp >= TimeLevelReg)) {
		/*
		fprintf(stderr,"CurrentBlockIrr = %llu\n",
				CurrentBlockIrr);
		fprintf(stderr,"CurrentBlockReg = %llu\n",
				CurrentBlockReg);
		fprintf(stderr,"TimeBlockReg    = %llu\n",
				TimeBlockReg);
		fprintf(stderr,"TimeBlockIrr    = %llu\n",
				TimeBlockIrr);
				*/
		//if (TimeLevelTmp > TimeLevelReg)
			//fprintf(stderr, "PID=%d, Irr=%d, Reg=%d\n", PID, TimeLevelTmp0, TimeLevelReg);

		TimeLevelTmp--;
		TimeBlockTmp *= 0.5;
	}

	TimeLevelIrr = TimeLevelTmp;

	if (TimeLevelIrr < time_block) {
		//std::cerr << "TimeLevelIrr is too small" << std::endl;
		TimeLevelIrr = std::max(time_block, TimeLevelIrr);
	}


	TimeStepIrr = static_cast<double>(pow(2, TimeLevelIrr));
	TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));

	if (TimeStepIrr*EnzoTimeStep*1e4<1e-11) {
		fprintf(stderr, "Too small TimeStepIrr! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				PID, TimeStepIrr*EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*EnzoTimeStep*1e4);
		while (TimeStepIrr*EnzoTimeStep*1e4<1e-11) {
			TimeLevelIrr++;
			TimeStepIrr  = static_cast<double>(pow(2, TimeLevelIrr));
			TimeBlockIrr = static_cast<ULL>(pow(2, TimeLevelIrr-time_block));
		}
	}

	if (TimeStepIrr > 1) {
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}


	/*
	if (TimeStepIrr*EnzoTimeStep*1e4 < KSTime && (this->isCMptcl == false))
		BinaryCandidateList.push_back(this);
	if (PID == 430) {
		std::cerr << "After TimeLevelIrr=" << TimeLevelIrr << std::endl;
		std::cerr << std::endl;
	}
	*/
}



// Update TimeStepReg // need to review
void Particle::calculateTimeStepReg() {
	//fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepTmp;
	ULL TimeBlockTmp;
	int TimeLevelTmp, TimeLevelTmp0;

	getBlockTimeStep(getNewTimeStepReg(Velocity, a_reg), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);

	//fprintf(stderr, "in CalReg, raw time step=%.2eMyr, ", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);


	//std::cout << "NBODY+: TimeStepRegTmp = " << TimeStepTmp << std::endl;

	TimeLevelTmp0 = TimeLevelTmp;

	if (TimeLevelTmp >= TimeLevelReg+1) {
		if (fmod(CurrentBlockReg, 2*TimeBlockReg)==0 \
				&& CurrentTimeReg != 0) {
			TimeLevelTmp   = TimeLevelReg+1;
			TimeBlockTmp   = TimeBlockReg*2;

			while ((TimeLevelTmp0 > TimeLevelTmp) \
					&& (fmod(CurrentBlockReg, 2*TimeBlockTmp) == 0)) {
				TimeBlockTmp = 2*TimeBlockTmp;
				TimeLevelTmp   = TimeLevelTmp + 1;
			}
		}
		else {
			TimeLevelTmp   = TimeLevelReg;
		}
	}
	else if (TimeLevelTmp < TimeLevelReg) {
		TimeLevelTmp = TimeLevelReg - 1;
		if (TimeLevelTmp0 < TimeLevelTmp)
			TimeLevelTmp--;
	}
	else {
		TimeLevelTmp = TimeLevelReg;
	}


	// update needed. regcor_gpu.for:725 (Makino, ApJ, 369)
	/*
	if (TimeStepRegTmp > 0.1 && TimeStepRegTmp > TimeStepReg) {
		double v2 = 0., a2=0., dt;
		for (int dim=0; dim<Dim; dim++) {
			v2 += (PredVelocity[dim]-NewVelocity[dim])*(PredVelocity[dim]-NewVelocity[dim]);
			a2 += a_reg[dim][0]*a_reg[dim][0];
		}
		dt = TimeStepReg*std::pow((1e-4*TimeStepReg*TimeStepReg*a2/v2),0.1);
		if (dt < TimeStepRegTmp) {
			TimeStepRegTmp = TimeStepReg;
		}	
	}
	*/

	//fprintf(stderr, " final time step=%.2eMyr\n", TimeStepRegTmp*EnzoTimeStep*1e10/1e6);

	TimeLevelReg = std::max(time_block,TimeLevelTmp);
	//TimeLevelReg = std::max(time_block, TimeLevelReg);

	if (this->NumberOfNeighbor == 0) {
		TimeLevelIrr = TimeLevelReg;
	}

	TimeStepReg  = static_cast<double>(pow(2, TimeLevelReg));
	TimeBlockReg = static_cast<ULL>(pow(2, TimeLevelReg-time_block));

	if (TimeStepReg*EnzoTimeStep*1e4 < 1e-7) {
		fprintf(stderr, "PID: %d, TimeStep = %.3e, TimeStepTmp0 = %.3e\n",
			 	PID, TimeStepReg*EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*EnzoTimeStep*1e4);
		fflush(stderr);
	}

	if (CurrentTimeReg+TimeStepReg > 1 && CurrentTimeReg != 1.0) {
		TimeStepReg = 1 - CurrentTimeReg;
		TimeBlockReg = block_max-CurrentBlockReg;
	}
	/*
	if (TimeStepReg*EnzoTimeStep*1e4<1e-9) {
		fprintf(stderr, "PID: %d, TimeStep = %.3e, TimeStepTmp0 = %.3e\n",
			 	PID, TimeStepReg*EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*EnzoTimeStep*1e4);
		throw std::runtime_error("TimeStepReg is too small.");
	}
	*/
	if (TimeStepReg*EnzoTimeStep*1e4<1e-9) {
		fprintf(stderr, "Too small TimeStepReg! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				PID, TimeStepReg*EnzoTimeStep*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*EnzoTimeStep*1e4);
		while (TimeStepReg*EnzoTimeStep*1e4<1e-9) {
			TimeLevelReg++;
			TimeStepReg  = static_cast<double>(pow(2, TimeLevelReg));
			TimeBlockReg = static_cast<ULL>(pow(2, TimeLevelReg-time_block));
		}
	}
	if (TimeStepReg > 1) {
		fprintf(stderr, "TimeStepReg=%e, TimeLevelReg=%d, TimeLevelTmp0=%d\n",TimeStepReg, TimeLevelReg, TimeLevelTmp0);
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",TimeStepIrr, TimeLevelIrr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}

	//std::cout << "NBODY+: TimeStepReg = " << TimeStepReg << std::endl;
}










