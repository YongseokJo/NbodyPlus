#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include "../global.h"
#include "../def.h"
#include <unordered_set>


void Particle::computeAccelerationIrr() {

	this->NewNumberOfNeighbor = 0; // for Few-body Search by EW 2025.3.1

	if (this->NumberOfNeighbor == 0) {
		for (int dim=0; dim<Dim; dim++){
			this->NewPosition[dim] = this->Position[dim];
			this->NewVelocity[dim] = this->Velocity[dim];
		}
		return;
	}

	double dt, mdot, epsilon=1e-6;
	double new_time; // 0 for current and 1 for advanced times

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[Dim], adot_tmp[Dim]; // 0 for current and 1 for predicted accelerations
	double pos[Dim], vel[Dim];
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double m_r3;
	Particle* ptcl;
	new_time = this->CurrentTimeIrr + this->TimeStepIrr; // the time to be advanced to
	dt       = this->TimeStepIrr*EnzoTimeStep; // interval of time step


	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion Calculation
	 ********************************************************/
	this->predictParticleSecondOrder(this->TimeStepIrr, pos, vel);

	for (int i=0; i<this->NumberOfNeighbor; i++) {

		ptcl = &particles[Neighbors[this->ParticleIndex * MaxNumNeighbor + i]];

		if (!ptcl->isActive) {
			if (ptcl->CMPtclIndex != -1) {
				CMPtclsSet.insert(ptcl->CMPtclIndex);
			}
			continue;
		}


		/*
		if (ptcl->isCMptcl) {
			fprintf(stderr, "my = %d , pid of cm = %d\n", this->PID, ptcl->PID);
			fflush(stderr);
		}

	 if (ptcl->Position[0]!=ptcl->Position[0]) {
			fprintf(stderr, "Nan occurs, %lf", ptcl->Position[0]);
			fflush(stderr);
			assert(this->Position[0] ==  this->Position[0]);
			exit(EXIT_FAILURE);
	 }
		if (ptcl->PID == this->PID)  {
			fprintf(stderr, "Myself in neighbor (%d)", PID);
			fflush(stderr);
			exit(EXIT_FAILURE);
			continue;
		}
		*/

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		if (sqrt(r2) < RSEARCH/position_unit && vx < 0) {
			NewNeighbors[this->ParticleIndex * MaxNumNeighbor + this->NewNumberOfNeighbor] = Neighbors[this->ParticleIndex * MaxNumNeighbor + i];
			this->NewNumberOfNeighbor++;
		}

		//mdot = ptcl->evolveStarMass(CurrentTimeIrr,
				//CurrentTimeIrr+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
																													 //
																													 // add the contribution of jth particle to acceleration of current and predicted times

		m_r3 = ptcl->Mass/(r2*sqrt(r2));

		for (int dim=0; dim<Dim; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	} // endfor ptcl

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (this->PID == ptcl->PID) {
			continue;
		}

		if (!ptcl->isActive) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", this->PID, ptcl->PID);
			assert(ptcl->isActive);
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		if (sqrt(r2) < RSEARCH/position_unit && vx < 0) {
			NewNeighbors[this->ParticleIndex * MaxNumNeighbor + this->NewNumberOfNeighbor] = i;
			this->NewNumberOfNeighbor++;
		}

		//mdot = ptcl->evolveStarMass(CurrentTimeIrr,
				//CurrentTimeIrr+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
																													//
																													// add the contribution of jth particle to acceleration of current and predicted times

		m_r3 = ptcl->Mass/(r2*sqrt(r2));

		for (int dim=0; dim<Dim; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	}


	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;
	double dt_ex = (new_time - this->CurrentTimeReg)*EnzoTimeStep;

	double A, B, C;

	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt; // VERY IMPORTANT BUG FIXED by EW 2025.7.18

	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	for (int dim=0; dim<Dim; dim++) {


#define noOptimization // I don't think it works 
#ifdef Optimization
		// do the higher order correcteion


		A = -12   *( this->a_irr[dim][0]   - a_tmp[dim] );
		B = -4*dt *( 2*this->a_irr[dim][1] + adot_tmp[dim] );
		C =  6*dt *( this->a_irr[dim][1]   + adot_tmp[dim] );

		a2 = dt *(A+B)/48;
		a3 = dt *(C-A)/120;

		//fprintf(stderr, "da_dt2=%.2e, adot_dt=%.2e, a2=%.2e, a3=%.2e\n", da_dt2, adot_dt, a2, a3);

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt + a3*dt;
		this->NewVelocity[dim] = vel[dim] + 4*a2  + 5*a3;


		// note that these higher order terms and lowers have different neighbors
		this->a_irr[dim][0] = a_tmp[dim];
		this->a_irr[dim][1] = adot_tmp[dim];
		this->a_irr[dim][2] = a2*24/dt3;
		this->a_irr[dim][3] = a3*120/dt4;
#else
		// do the higher order correcteion
		da_dt2  = (this->a_irr[dim][0] - a_tmp[dim]) / dt2; 
		adot_dt = (this->a_irr[dim][1] + adot_tmp[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*this->a_irr[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->NewVelocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		// note that these higher order terms and lowers have different neighbors
		this->a_irr[dim][0] = a_tmp[dim];
		this->a_irr[dim][1] = adot_tmp[dim];
		this->a_irr[dim][2] = a2;
		this->a_irr[dim][3] = a3;
#endif
	}



	for (int dim=0; dim<Dim; dim++) {
		this->a_tot[dim][0] = this->a_reg[dim][0] + this->a_irr[dim][0] + this->a_reg[dim][1]*dt_ex; // affect the next
		this->a_tot[dim][1] = this->a_reg[dim][1] + this->a_irr[dim][1];
		this->a_tot[dim][2] = this->a_reg[dim][2] + this->a_irr[dim][2];
		this->a_tot[dim][3] = this->a_reg[dim][3] + this->a_irr[dim][3];
	}
}




// Modified by EW 2025.3.3 for no CUDA version

void Particle::computeAccelerationReg() {

	double dt, mdot, epsilon=1e-6;
	double new_time; // 0 for current and 1 for advanced times

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a[Dim], adot[Dim], a_new[Dim], adot_new[Dim]; // 0 for current and 1 for predicted accelerations
	double pos[Dim], vel[Dim];
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double m_r3;
	int    j=0;
	Particle* ptcl;
	new_time = this->CurrentTimeReg+this->TimeStepReg; // the time to be advanced to
	dt       = this->TimeStepReg*EnzoTimeStep; // interval of time step
	if (dt == 0.0) {
		dt  = this->TimeBlockReg*time_step*EnzoTimeStep; // interval of time step
		assert(dt != 0.0);
	}
	this->NewNumberOfNeighbor = 0;

	std::unordered_set<int> RealNeighbors; // Neighbors containing CM ptcls, not members
	RealNeighbors.reserve(this->NumberOfNeighbor);

	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a[dim]        = 0.0;
		adot[dim]     = 0.0;
		a_new[dim]    = 0.0;
		adot_new[dim] = 0.0;
		this->a_irr[dim][0] = 0.;
		this->a_irr[dim][1] = 0.;
	}

	for (int i=0; i<this->NumberOfNeighbor; i++) {
		ptcl = &particles[Neighbors[this->ParticleIndex * MaxNumNeighbor + i]];
		if (ptcl->isActive)
			RealNeighbors.insert(Neighbors[this->ParticleIndex * MaxNumNeighbor + i]);
		else if (ptcl->CMPtclIndex != -1)
			RealNeighbors.insert(ptcl->CMPtclIndex);
	}
	if (this->isCMptcl)
		RealNeighbors.erase(this->ParticleIndex);


	/*******************************************************
	 * Regular Acceleartion Calculation
	 ********************************************************/
	if (this->NumberOfNeighbor == 0)
		this->predictParticleSecondOrder(this->TimeStepReg, pos, vel);
	else
		this->predictParticleSecondOrder(0, pos, vel);


	for (int i=0; i<=global_variable->LastParticleIndex; i++) {
		ptcl = &particles[i];

		if (!ptcl->isActive || this->PID == ptcl->PID)
			continue;

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, pos_neighbor, vel_neighbor);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		//mdot = ptcl->evolveStarMass(CurrentTimeIrr,
		//CurrentTimeIrr+TimeStepIrr*1.01)/TimeStepIrr*1e-2; // derivative can be improved
		//
		// add the contribution of jth particle to acceleration of current and predicted times

		m_r3 = ptcl->Mass/r2/sqrt(r2);


		//std::cout << "PIDs=" <<  this->Neighbors[j] << ', ' << ptcl->PID << NumberOfNeighbor<< std::endl;
		//if (this->Neighbors[j] == ptcl->PID) {
		if (RealNeighbors.find(ptcl->ParticleIndex) != RealNeighbors.end()) {
			//std::cout << this->PID << ", PIDs=" <<  this->Neighbors[j] << ", " << ptcl->PID << std::endl;
			j++;
		} 
		else {
			for (int dim=0; dim<Dim; dim++){
				a[dim]    += m_r3*x[dim];
				adot[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}


		if (r2 < this->RadiusOfNeighbor) {
			if (!ptcl->isCMptcl) {
				NewNeighbors[this->ParticleIndex * MaxNumNeighbor + this->NewNumberOfNeighbor] = ptcl->ParticleIndex;
				this->NewNumberOfNeighbor++;
			}
			else {
				for (int k=0; k<ptcl->NumberOfMember; k++) {
					NewNeighbors[this->ParticleIndex * MaxNumNeighbor + this->NewNumberOfNeighbor] = ptcl->Members[k];
					this->NewNumberOfNeighbor++;
				}
			}
			for (int dim=0; dim<Dim; dim++){
				this->a_irr[dim][0] += m_r3*x[dim];
				this->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++){
				a_new[dim]    += m_r3*x[dim];
				adot_new[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
	} // endfor ptcl
	assert(j == RealNeighbors.size());

	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;

	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt;
	for (int dim=0; dim<Dim; dim++) {

		// do the higher order correcteion
		da_dt2  = (this->a_reg[dim][0] - a[dim]) / dt2;
		adot_dt = (this->a_reg[dim][1] + adot[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*this->a_reg[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		//fprintf(stderr, "DIM=%d, pid=%d, a=%.2e, da_dt2=%.2e, adot_dt=%.2e, a2=%.2e, a3=%.2e\n", dim,this->PID,a[dim], da_dt2, adot_dt, a2, a3);

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->NewVelocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;


		this->a_reg[dim][2] = a2;
		this->a_reg[dim][3] = a3;
	}


	for (int dim=0; dim<Dim; dim++) {
		this->a_reg[dim][0] = a_new[dim];
		this->a_reg[dim][1] = adot_new[dim];
		this->a_tot[dim][0] = this->a_reg[dim][0] + this->a_irr[dim][0]; 
		this->a_tot[dim][1] = this->a_reg[dim][1] + this->a_irr[dim][1];
		if (this->NewNumberOfNeighbor == 0) {
			this->a_tot[dim][2] = this->a_reg[dim][2];
			this->a_tot[dim][3] = this->a_reg[dim][3];
		}
	}
}




// Modified by EW 2025.1.30

void Particle::updateRegularParticleCuda() {

	double new_a[Dim], new_adot[Dim];
	double new_time = this->CurrentTimeReg+this->TimeStepReg;
	double pos[Dim], vel[Dim];
	if (this->NumberOfNeighbor == 0)
		this->predictParticleSecondOrder(this->TimeStepReg, pos, vel);
	else
		this->predictParticleSecondOrder(0, pos, vel);

	double a_tmp[Dim], adot_tmp[Dim];

	for (int dim=0; dim<Dim; dim++) {
		new_a[dim]			= this->a_irr[dim][0];
		new_adot[dim]		= this->a_irr[dim][1];
		this->a_irr[dim][0] = 0.;
		this->a_irr[dim][1] = 0.;

		a_tmp[dim]          = 0.;
		adot_tmp[dim]       = 0.;
	}

	std::unordered_set<int> hashTableOld;
	hashTableOld.reserve(this->NumberOfNeighbor);
	std::unordered_set<int> hashTableNew;
	hashTableNew.reserve(this->NewNumberOfNeighbor); // We are including myself in neighbor from GPU kernel by EW 2025.8.26


	hashTableNew.insert(NewNeighbors + this->ParticleIndex * MaxNumNeighbor, NewNeighbors + this->ParticleIndex * MaxNumNeighbor + this->NewNumberOfNeighbor);
	hashTableNew.erase(this->ParticleIndex);
	this->NewNumberOfNeighbor--;

	for (int i = 0; i < this->NumberOfNeighbor; i++) {
		if (particles[Neighbors[this->ParticleIndex * MaxNumNeighbor + i]].isActive)
			hashTableOld.insert(Neighbors[this->ParticleIndex * MaxNumNeighbor + i]);
		else if (particles[Neighbors[this->ParticleIndex * MaxNumNeighbor + i]].isCMptcl)
			hashTableOld.insert(particles[Neighbors[this->ParticleIndex * MaxNumNeighbor + i]].CMPtclIndex);
	}

	Particle* ptcl;
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double dx[Dim], dv[Dim];
	double dr2;
	double dxdv;
	double m_r3;

	for (int _OldNeighborIndex: hashTableOld) {
		ptcl = &particles[_OldNeighborIndex];

		// neighbor in old but not in new
		if ( hashTableNew.find(ptcl->ParticleIndex) == hashTableNew.end() ) {

			if (ptcl->NumberOfNeighbor == 0)
				ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, pos_neighbor, vel_neighbor);
			else
				ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

			dr2  = 0.0;
			dxdv = 0.0;
			for (int dim=0; dim<Dim; dim++) {
				dx[dim] = pos_neighbor[dim] - pos[dim];
				dv[dim] = vel_neighbor[dim] - vel[dim];
				dr2    += dx[dim]*dx[dim];
				dxdv   += dx[dim]*dv[dim];
			}

			m_r3 = ptcl->Mass/dr2/sqrt(dr2);

			for (int dim=0; dim<Dim; dim++){
				a_tmp[dim]    -= m_r3*dx[dim];
				adot_tmp[dim] -= m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
			}
		}
	}

	for (int _NewNeighborIndex: hashTableNew) {
		ptcl = &particles[_NewNeighborIndex];

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, pos_neighbor, vel_neighbor);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		dr2  = 0.0;
		dxdv = 0.0;
		for (int dim=0; dim<Dim; dim++) {
			dx[dim] = pos_neighbor[dim] - pos[dim];
			dv[dim] = vel_neighbor[dim] - vel[dim];
			dr2    += dx[dim]*dx[dim];
			dxdv   += dx[dim]*dv[dim];
		}

		m_r3 = ptcl->Mass/dr2/sqrt(dr2);

		for (int dim=0; dim<Dim; dim++){
			a_irr[dim][0] += m_r3*dx[dim];
			a_irr[dim][1] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
		}

		// neighbor in new but not in old
		if ( hashTableOld.find(ptcl->ParticleIndex) == hashTableOld.end() ) {
			for (int dim=0; dim<Dim; dim++){
				a_tmp[dim]    += m_r3*dx[dim];
				adot_tmp[dim] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
			}
		}
	}


	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	double dt  = this->TimeStepReg*EnzoTimeStep;  // unit conversion
	if (dt == 0.0) {
		dt  = this->TimeBlockReg*time_step*EnzoTimeStep; // unit conversion
		assert(dt != 0.0);
	}
	double dt2 = dt*dt;
	double dt3 = dt2*dt;
	double dt4 = dt3*dt;
	double dt5 = dt4*dt;
	double da_dt2, adot_dt, a2, a3;

	//fprintf(stdout, "PID=%d\n", ptcl->PID);
	for (int dim=0; dim<Dim; dim++) {
		da_dt2  = (this->a_reg[dim][0] - new_a[dim]   - a_tmp[dim]    ) / dt2;
		adot_dt = (this->a_reg[dim][1] + new_adot[dim] + adot_tmp[dim]) / dt;


		a2 =  -6*da_dt2 - 2*adot_dt - 2*this->a_reg[dim][1]/dt;
		a3 = (12*da_dt2 + 6*adot_dt)/dt;

		// note that these higher order terms and lowers have different neighbors

		// 4th order correction
		// save the values in the temporary variables
		this->NewPosition[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->NewVelocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		this->a_reg[dim][2] = a2;
		this->a_reg[dim][3] = a3;

		// reset for future use
		a_tmp[dim]    = 0.;
		adot_tmp[dim] = 0.;
	}

	int _NewNumberOfNeighbor = 0;
	for (int _NewNeighborIndex: hashTableNew) {
		ptcl = &particles[_NewNeighborIndex];
		if (ptcl->isCMptcl) {
			for (int j=0; j<ptcl->NumberOfMember; j++) {
				NewNeighbors[this->ParticleIndex * MaxNumNeighbor + _NewNumberOfNeighbor] = ptcl->Members[j];
				_NewNumberOfNeighbor++;
			}
		}
		else {
			NewNeighbors[this->ParticleIndex * MaxNumNeighbor + _NewNumberOfNeighbor] = ptcl->ParticleIndex;
			_NewNumberOfNeighbor++;
		}
	}
	this->NewNumberOfNeighbor = _NewNumberOfNeighbor;


	for (int dim=0; dim<Dim; dim++) {
		this->a_reg[dim][0] = new_a[dim];
		this->a_reg[dim][1] = new_adot[dim];
		this->a_tot[dim][0] = this->a_reg[dim][0] + this->a_irr[dim][0];
		this->a_tot[dim][1] = this->a_reg[dim][1] + this->a_irr[dim][1];
		if (this->NewNumberOfNeighbor == 0) {
			this->a_tot[dim][2] = this->a_reg[dim][2];
			this->a_tot[dim][3] = this->a_reg[dim][3];
		}
	}
}

