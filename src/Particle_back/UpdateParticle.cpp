#include "../global.h"
#include "../defs.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>





/*
 *  Purporse: Predict particle positions and velocities up to second order 
 *            using a_0 and d(a_0)/dt; refer to Nbody6++ manual Eq. 8.3 and 8.4
 *
 *  Date    : 2024.04.29  by Seoyoung Kim
 *  Date    : 2024.01.09  by Seoyoung Kim
 *  Modified: 2024.01.10  by Yongseok Jo
 *
 */
void Particle::predictParticleSecondOrder(REAL time) {

	// Doubling check
	// temporary variables for calculation
	if (this->CurrentTimeReg == time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}

	REAL dt;
	//dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;
	dt = (time - this->CurrentTimeReg)*EnzoTimeStep;
	if (dt < 0) {
		fprintf(stderr, "PID=%d, dt=%e",PID, dt);
		throw std::runtime_error("SecondOrder");
	}

	// only predict the positions if necessary
	// how about using polynomial correction here?
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim];
		PredVelocity[dim] =  (a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim];
	}

	return;
}



void Particle::predictParticleSecondOrderIrr(REAL time) {

	// Doubling check
	// temporary variables for calculation
	if (this->CurrentTimeIrr == time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}

	REAL dt;
	//dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;
	dt = (time - this->CurrentTimeIrr)*EnzoTimeStep;
	/*
	if (dt < 0 || this->CurrentTimeIrr>1) {
		fprintf(stderr, "PID=%d, time=%e, timeirr=%lf, timereg=%lf, dt=%lf, dtirr=%lf, dtreg=%lf\n",
			 	this->PID, time, this->CurrentTimeIrr, this->CurrentTimeReg, dt, this->TimeStepIrr, this->TimeStepReg);
		fprintf(stderr, "PID=%d, time=%e, timeirr=%llu, timereg=%llu, dt=%lf, dtirr=%llu, dtreg=%llu\n",
			 	PID, time, CurrentBlockIrr, CurrentBlockReg, dt, TimeBlockIrr, TimeBlockReg);
		fflush(stderr);
		fflush(stdout);
		//throw std::runtime_error("SecondOrderIrr");
	}
	*/
	// only predict the positions if necessary
	// how about using polynomial correction here?
	for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = ((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim];
		PredVelocity[dim] =  (a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim];
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
void Particle::correctParticleFourthOrder(REAL current_time, REAL next_time, REAL a[3][4]) {

	if (current_time == next_time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}

	REAL dt;
	REAL dt3,dt4,dt5;

	dt = (next_time - current_time)*EnzoTimeStep;

	dt3 = dt*dt*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;

	// correct the predicted values positions and velocities at next_time
	// and save the corrected values to particle positions and velocities
	// the latest values of a2dot and a3dots (obtained from hermite method) are used
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] = PredPosition[dim]+ a[dim][2]*dt4/24 + a[dim][3]*dt5/120;
		Velocity[dim] = PredVelocity[dim]+ a[dim][2]*dt3/6  + a[dim][3]*dt4/24;
	}
}



void Particle::polynomialPrediction(REAL current_time) {

}


void Particle::updateParticle() {
	for (int dim=0; dim<Dim; dim++) {
		Position[dim] = NewPosition[dim];
		Velocity[dim] = NewVelocity[dim];
	}
}

/*
void Particle::updateParticle(REAL current_time, REAL next_time, REAL a[3][4]) {

	if (current_time == next_time) {
		for (int dim=0; dim<Dim; dim++) {
			PredPosition[dim] = Position[dim];
			PredVelocity[dim] = Velocity[dim];
		}
		return;
	}
	predictParticleSecondOrder(current_time, next_time, a);
	correctParticleFourthOrder(current_time, next_time, a);
	Mass += evolveStarMass(current_time, next_time); // CurrentTimeIrr
}*/



void Particle::UpdateRadius() {

	if (LocalDensity == 0) {
		RadiusOfAC = 0.11;
		//RadiusOfAC = InitialRadiusOfAC;
		//RadiusOfAC = 1.00;
		LocalDensity = 10;
	}
	else {
		/* exponential (aggressive) */
		/*
			 const REAL c = 0.5;
			 const REAL b = std::log(2) / (NumNeighborMax);  // ln(2) / 40
			 REAL exp = a * (std::exp(b * NumberOfAC) - 1);
			 */

		/* n=2 polynomial (mild) as n increases it grows mild */
		
		if (NumberOfAC > FixNumNeighbor) {
			const int n = 2;
			const REAL c = (NumNeighborMax-FixNumNeighbor);
			const REAL b = 0.9 / std::pow(c,n);  // ln(2) / 40
			REAL x = NumberOfAC-FixNumNeighbor;
			REAL a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
			//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
				 	//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
			RadiusOfAC *= (1.-a);
		}
		else if (NumberOfAC < FixNumNeighbor) {
			const int n = 3;
			const REAL c = (NumNeighborMax-FixNumNeighbor);
			const REAL b = 0.5 / std::pow(c,n);  // ln(2) / 40
			REAL x = NumberOfAC-FixNumNeighbor;
			REAL a = n%2==0 ? b*std::abs(x)*std::pow(x,n-1) : b*std::pow(x,n);
			//fprintf(stdout, "PID=%d, NumberOfAC=%d, 1-a=%e, R0=%e(%e), R=%e(%e)\n",
				 	//PID,NumberOfAC, 1.-a, RadiusOfAC, RadiusOfAC*RadiusOfAC, RadiusOfAC*(1-a),RadiusOfAC*(1-a)*RadiusOfAC*(1-a));
			RadiusOfAC *= (1.-a);
		}
	}

	/*
	REAL MeanRadius=0, TotalMass=0, LocalDensity0=0;
	for (Particle *neighbor:ACList) {
		MeanRadius += neighbor->Mass*dist(Position,neighbor->Position);
		TotalMass  += neighbor->Mass;
	}
	MeanRadius        /= TotalMass/1.2;  // this is the mean radius with some factor
	LocalDensity0      = TotalMass/std::pow(MeanRadius,3.0);

	if (LocalDensity == 0) {
		this->RadiusOfAC   = std::pow(LocalDensity0/(this->Mass*FixNumNeighbor), -1.0/3.0); 

		if (LocalDensity0 > 1.3*LocalDensity) 
			this->RadiusOfAC *= 0.9;
		else if (LocalDensity0 < 0.7*LocalDensity) 
			this->RadiusOfAC *= 1.1;

		if (NumberOfAC < FixNumNeighbor/2)
			this->RadiusOfAC *= 1.2;

		if (NumberOfAC > FixNumNeighbor*2)
			this->RadiusOfAC *= 0.8;
	}
	*/


	//this->LocalDensity = LocalDensity0;

	//this->RadiusOfAC   = std::pow(LocalDensity0/(this->Mass*FixNumNeighbor), -1.0/3.0); 
	//fprintf(stdout, "PID %d : TWR = %.3e, TM=%.3e, LD0=%.3e, RAC=%.3e, mass=%.3e\n", 
			//PID, MeanRadius, TotalMass, LocalDensity0, RadiusOfAC, Mass);
	//fflush(stdout); 
}


void Particle::UpdateNeighbor(std::vector<Particle*> &particle) {
	bool isExcess=false;
	ACList.clear();
	NumberOfAC = 0;
	for (Particle* ptcl:particle) {

		if (ptcl->PID==PID)
			continue;

		if (dist(Position, ptcl->Position)<this->RadiusOfAC) {
			ACList.push_back(ptcl);	
			NumberOfAC++;
			if (NumberOfAC >= NumNeighborMax) {
				isExcess = true; 
				break;
			}
		}
	}

	if (isExcess) {
		std::cerr << "Number of neighbors exceeded." << std::endl;
		this->RadiusOfAC *= 0.8;
		UpdateNeighbor(particle);	
	}
	else {
		std::sort(ACList.begin(),ACList.end(),
				[](Particle* p1, Particle* p2) { return p1->ParticleOrder < p2->ParticleOrder;});
	}
}






