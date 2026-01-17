#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include "../global.h"
#include "../def.h"

/* // commented out by EW 2025.7.15 to define this function as inline function
template <typename T>
void Particle::predictParticleSecondOrder(double dt, T pos[], T vel[]) {
	// Doubling check
	// temporary variables for calculation

	// only predict the positions if necessary
	// how about using polynomial correction here?
	
	dt = dt*enzo_time_step;

	if (dt == 0) {
		for (int dim=0; dim<DIM; dim++) {
			pos[dim] = static_cast<T>(Position[dim]);
			vel[dim] = static_cast<T>(Velocity[dim]);
		}
	}
	else {
		for (int dim=0; dim<DIM; dim++) {
			pos[dim] = static_cast<T>(((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim]);
			vel[dim] = static_cast<T>((a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim]);
		}
	}
	return;
}
*/


/*
 *  Purporse: Correct particle positions and velocities up to fourth order
 *            using a_p and d(a_p)/dt; refer to Nbody6++ manual Eq. 8.9 and 8.10
 *
 *  Date    : 2024.01.10  by Yongseok Jo
 *  Modified: 2024.01.11  by Seoyoung Kim
 *
 */
void Particle::correct_particle_fourth_order(double dt, double pos[], double vel[], double a[3][4]) {
	double dt3,dt4,dt5;

	dt = dt*enzo_time_step;

	dt3 = dt*dt*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	// correct the predicted values positions and velocities at next_time
	// and save the corrected values to particle positions and velocities
	// the latest values of a2dot and a3dots (obtained from hermite method) are used
	for (int dim=0; dim<DIM; dim++) {
		new_position[dim] = pos[dim]+ a[dim][2]*dt4/24 + a[dim][3]*dt5/120;
		new_velocity[dim] = vel[dim]+ a[dim][2]*dt3/6  + a[dim][3]*dt4/24;
	}
}


/*
void Particle::polynomialPrediction(double current_time) {

}
*/

/* // commented out by EW 2025.3.3 to define this function as inline function
void Particle::updateParticle() {
	
	for (int dim=0; dim<DIM; dim++) {
		this->position[dim] = this->new_position[dim];
		this->velocity[dim] = this->new_velocity[dim];
	}
	
	//updateTimeStep();
}
*/


void Particle::update_radius() {

	/* exponential (aggressive) */
	/*
	const double c = 0.5;
	const double b = std::log(2) / (MAX_NUM_NEIGHBOR);  // ln(2) / 40
	double exp = a * (std::exp(b * NumberOfAC) - 1);
	*/

	/* n=2 polynomial (mild) as n increases it grows mild */

	const double MaxRadius2 = MAX_NEIGHBOR_RADIUS*MAX_NEIGHBOR_RADIUS/(position_unit*position_unit);

	if (this->num_neighbors > fixed_num_neighbors) {
		const int n = 2;
		const double c = (MAX_NUM_NEIGHBOR-fixed_num_neighbors);
		const double b = 0.9 / std::pow(c,n);  // ln(2) / 40
		double x = this->num_neighbors-fixed_num_neighbors;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		this->neighbor_radius_sq *= (1.-a);
		if (this->neighbor_radius_sq > MaxRadius2)
			this->neighbor_radius_sq = MaxRadius2;
	}
	else if (this->num_neighbors < fixed_num_neighbors) {
		const int n = 3;
		const double c = (MAX_NUM_NEIGHBOR-fixed_num_neighbors);
		const double b = 0.5 / std::pow(c,n);  // ln(2) / 40
		double x = this->num_neighbors-fixed_num_neighbors;
		double a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
		//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
		//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
		this->neighbor_radius_sq *= (1.-a);
		if (this->neighbor_radius_sq > MaxRadius2)
			this->neighbor_radius_sq = MaxRadius2;
	}
}







double getNewTimeStepReg(double v[3], double f[3][4]);
double getNewTimeStepIrr(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, ull_t &TimeBlock, double &TimeStep);

void Particle::calculate_time_step_irr() {
	double TimeStepTmp;
	int TimeLevelTmp, TimeLevelTmp0;
	ull_t TimeBlockTmp;

	if (this->num_neighbors == 0) {
		time_level_irr = time_level_reg;
		time_step_irr = static_cast<double>(pow(2, time_level_reg));
		time_block_irr = static_cast<ull_t>(pow(2, time_level_reg-time_block));
		return;
	}

	getBlockTimeStep(getNewTimeStepIrr(acc_total, acc_irregular), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	/*
	if (next_reg_time_block < TimeBlockReg && isRegular) {
		fprintf(stderr, "PID=%d, next_reg_time_block=%llu, TimeBlockReg=%llu\n", PID, next_reg_time_block, TimeBlockReg);
	}
	*/

	if (TimeLevelTmp > time_level_irr) {
		if (fmod(new_current_block_irr, 2*time_block_irr)==0) {
			TimeLevelTmp = time_level_irr+1;
			TimeBlockTmp = 2*time_block_irr;
		}
		else {
			TimeLevelTmp = time_level_irr;
			TimeBlockTmp = time_block_irr;
		}
	}
	else if (TimeLevelTmp < time_level_irr) {
		if (TimeLevelTmp < time_level_irr-1) {
			TimeLevelTmp = time_level_irr - 2;
			TimeBlockTmp = time_block_irr/4;
		}
		else {
			TimeLevelTmp = time_level_irr - 1;
			TimeBlockTmp = time_block_irr/2;
		}
	} else {
		TimeLevelTmp = time_level_irr;
		TimeBlockTmp = time_block_irr;
	}


	while (((new_current_block_irr < current_block_reg+time_block_reg) && (new_current_block_irr+TimeBlockTmp > current_block_reg+time_block_reg)) || (TimeLevelTmp >= time_level_reg)) {
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

	time_level_irr = TimeLevelTmp;

	if (time_level_irr < time_block) {
		//std::cerr << "TimeLevelIrr is too small" << std::endl;
		time_level_irr = std::max(time_block, time_level_irr);
	}


	time_step_irr = static_cast<double>(pow(2, time_level_irr));
	time_block_irr = static_cast<ull_t>(pow(2, time_level_irr-time_block));

	if (time_step_irr*enzo_time_step*1e4<1e-11) {
		fprintf(stdout, "Too small TimeStepIrr! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				pid, time_step_irr*enzo_time_step*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*enzo_time_step*1e4);
		while (time_step_irr*enzo_time_step*1e4<1e-11) {
			time_level_irr++;
			time_step_irr  = static_cast<double>(pow(2, time_level_irr));
			time_block_irr = static_cast<ull_t>(pow(2, time_level_irr-time_block));
		}
	}

	if (time_step_irr > 1) {
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",time_step_irr, time_level_irr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}


	/*
	if (TimeStepIrr*enzo_time_step*1e4 < KSTime && (this->is_cm_particle == false))
		BinaryCandidateList.push_back(this);
	if (PID == 430) {
		std::cerr << "After TimeLevelIrr=" << TimeLevelIrr << std::endl;
		std::cerr << std::endl;
	}
	*/
}

void Particle::calculate_time_step_irr_v2() {
	double TimeStepTmp;
	int TimeLevelTmp, TimeLevelTmp0;
	ull_t TimeBlockTmp;

	if (this->num_neighbors == 0) {
		time_level_irr = time_level_reg;
		time_step_irr = static_cast<double>(pow(2, time_level_reg));
		time_block_irr = static_cast<ull_t>(pow(2, time_level_reg-time_block));
		return;
	}

	getBlockTimeStep(getNewTimeStepIrr(acc_total, acc_irregular), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);
	TimeLevelTmp0 = TimeLevelTmp;

	/*
	if (next_reg_time_block < TimeBlockReg && isRegular) {
		fprintf(stderr, "PID=%d, next_reg_time_block=%llu, TimeBlockReg=%llu\n", PID, next_reg_time_block, TimeBlockReg);
	}
	*/


	while (((new_current_block_irr < current_block_reg+time_block_reg) && (new_current_block_irr+TimeBlockTmp > current_block_reg+time_block_reg)) || (TimeLevelTmp >= time_level_reg)) {
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

	time_level_irr = TimeLevelTmp;

	if (time_level_irr < time_block) {
		//std::cerr << "TimeLevelIrr is too small" << std::endl;
		time_level_irr = std::max(time_block, time_level_irr);
	}


	time_step_irr = static_cast<double>(pow(2, time_level_irr));
	time_block_irr = static_cast<ull_t>(pow(2, time_level_irr-time_block));

	if (time_step_irr*enzo_time_step*1e4<1e-11) {
		fprintf(stdout, "Too small TimeStepIrr! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				pid, time_step_irr*enzo_time_step*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*enzo_time_step*1e4);
		while (time_step_irr*enzo_time_step*1e4<1e-11) {
			time_level_irr++;
			time_step_irr  = static_cast<double>(pow(2, time_level_irr));
			time_block_irr = static_cast<ull_t>(pow(2, time_level_irr-time_block));
		}
	}

	if (time_step_irr > 1) {
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",time_step_irr, time_level_irr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}


	/*
	if (TimeStepIrr*enzo_time_step*1e4 < KSTime && (this->is_cm_particle == false))
		BinaryCandidateList.push_back(this);
	if (PID == 430) {
		std::cerr << "After TimeLevelIrr=" << TimeLevelIrr << std::endl;
		std::cerr << std::endl;
	}
	*/
}



// Update TimeStepReg // need to review
void Particle::calculate_time_step_reg() {
	//fprintf(stdout, "Number of AC=%d\n", NumberOfAC);
	//std::cout << NumberOfAC << std::flush;
	double TimeStepTmp;
	ull_t TimeBlockTmp;
	int TimeLevelTmp, TimeLevelTmp0;

	getBlockTimeStep(getNewTimeStepReg(velocity, acc_regular), TimeLevelTmp, TimeBlockTmp, TimeStepTmp);

	//fprintf(stderr, "in CalReg, raw time step=%.2eMyr, ", TimeStepRegTmp*enzo_time_step*1e10/1e6);


	//std::cout << "NBODY+: TimeStepRegTmp = " << TimeStepTmp << std::endl;

	TimeLevelTmp0 = TimeLevelTmp;

	if (TimeLevelTmp >= time_level_reg+1) {
		if (fmod(current_block_reg, 2*time_block_reg)==0 \
				&& current_time_reg != 0) {
			TimeLevelTmp   = time_level_reg+1;
			TimeBlockTmp   = time_block_reg*2;

			while ((TimeLevelTmp0 > TimeLevelTmp) \
					&& (fmod(current_block_reg, 2*TimeBlockTmp) == 0)) {
				TimeBlockTmp = 2*TimeBlockTmp;
				TimeLevelTmp   = TimeLevelTmp + 1;
			}
		}
		else {
			TimeLevelTmp   = time_level_reg;
		}
	}
	else if (TimeLevelTmp < time_level_reg) {
		TimeLevelTmp = time_level_reg - 1;
		if (TimeLevelTmp0 < TimeLevelTmp)
			TimeLevelTmp--;
	}
	else {
		TimeLevelTmp = time_level_reg;
	}


	// update needed. regcor_gpu.for:725 (Makino, ApJ, 369)
	/*
	if (TimeStepRegTmp > 0.1 && TimeStepRegTmp > TimeStepReg) {
		double v2 = 0., a2=0., dt;
		for (int dim=0; dim<DIM; dim++) {
			v2 += (PredVelocity[dim]-NewVelocity[dim])*(PredVelocity[dim]-NewVelocity[dim]);
			a2 += a_reg[dim][0]*a_reg[dim][0];
		}
		dt = TimeStepReg*std::pow((1e-4*TimeStepReg*TimeStepReg*a2/v2),0.1);
		if (dt < TimeStepRegTmp) {
			TimeStepRegTmp = TimeStepReg;
		}	
	}
	*/

	//fprintf(stderr, " final time step=%.2eMyr\n", TimeStepRegTmp*enzo_time_step*1e10/1e6);

	time_level_reg = std::max(time_block,TimeLevelTmp);
	//TimeLevelReg = std::max(time_block, TimeLevelReg);

	if (this->num_neighbors == 0) {
		time_level_irr = time_level_reg;
	}

	time_step_reg  = static_cast<double>(pow(2, time_level_reg));
	time_block_reg = static_cast<ull_t>(pow(2, time_level_reg-time_block));
// /* // original code by EW 2025.3.17
	if (time_step_reg*enzo_time_step*1e4 < 1e-7) {
		fprintf(stdout, "PID: %d, TimeStep = %.3e, TimeStepTmp0 = %.3e\n",
			 	pid, time_step_reg*enzo_time_step*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*enzo_time_step*1e4);
		fflush(stderr);
	}
// */
/* // modified code by EW 2025.3.17
	if (TimeStepReg*enzo_time_step*1e4 < 1e-5) {
		while (TimeStepReg*enzo_time_step*1e4<1e-5) {
			TimeLevelReg++;
			TimeStepReg  = static_cast<double>(pow(2, TimeLevelReg));
			TimeBlockReg = static_cast<ull_t>(pow(2, TimeLevelReg-time_block));
		}
	}
*/

	if (current_time_reg+time_step_reg > 1 && current_time_reg != 1.0) {
		time_step_reg = 1 - current_time_reg;
		time_block_reg = block_max-current_block_reg;
	}
	/*
	if (TimeStepReg*enzo_time_step*1e4<1e-9) {
		fprintf(stderr, "PID: %d, TimeStep = %.3e, TimeStepTmp0 = %.3e\n",
			 	PID, TimeStepReg*enzo_time_step*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*enzo_time_step*1e4);
		throw std::runtime_error("TimeStepReg is too small.");
	}
	*/
	if (time_step_reg*enzo_time_step*1e4<1e-9) {
		fprintf(stdout, "Too small TimeStepReg! PID: %d, TimeStep = %e, TimeStepTmp0 = %e\n",
				pid, time_step_reg*enzo_time_step*1e4, static_cast<double>(pow(2, TimeLevelTmp0))*enzo_time_step*1e4);
		while (time_step_reg*enzo_time_step*1e4<1e-9) {
			time_level_reg++;
			time_step_reg  = static_cast<double>(pow(2, time_level_reg));
			time_block_reg = static_cast<ull_t>(pow(2, time_level_reg-time_block));
		}
	}
	if (time_step_reg > 1) {
		fprintf(stderr, "TimeStepReg=%e, TimeLevelReg=%d, TimeLevelTmp0=%d\n",time_step_reg, time_level_reg, TimeLevelTmp0);
		fprintf(stderr, "TimeStepIrr=%e, TimeLevelIrr=%d, TimeLevelTmp0=%d\n",time_step_irr, time_level_irr, TimeLevelTmp0);
		fflush(stderr);
		throw std::runtime_error("");
	}

	//std::cout << "NBODY+: TimeStepReg = " << TimeStepReg << std::endl;
}










