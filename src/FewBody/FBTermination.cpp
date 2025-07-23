#ifdef FEWBODY
#include "../global.h"
#include <unordered_set>

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void computeMemberAcceleration01(Particle* members);
void computeMemberAcceleration23(Particle* members);
void computeMemberAccelerationIrr(Particle* members, double new_time);

void FBTermination(Particle* ptclCM) {

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination.cpp... (CM PID: %d)\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);
	fprintf(binout, "CurrentTimeIrr of the first member (Myr): %e\n", particles[ptclCM->Members[0]].CurrentTimeIrr*EnzoTimeStep*1e4);
	fprintf(binout, "N_member: %d\n", ptclCM->NumberOfMember);

	Particle* members;

	ptclCM->isActive = false;
	for (int i = 0; i < ptclCM->NumberOfMember; i++) {
		members = &particles[ptclCM->Members[i]];
		members->CMPtclIndex = -1;
		if (members->Mass > 0.0)
			members->isActive = true;
	}

	NumberOfParticle--; // CM particle should be inactive by EW 2025.1.20

	// Set member information from CM ptcl; isActive = true & CMPtclIndex = -1 in RootRoutines.cpp
	for (int i = 0; i < ptclCM->NumberOfMember; i++) {
		members = &particles[ptclCM->Members[i]];
		if (members->Mass < 0.0)
			continue;
		
		NumberOfParticle++;

		if (ptclCM->NumberOfMember == 2)
			members->setBinaryInterruptState(BinaryInterruptState::none);
		else if (ptclCM->NumberOfMember > 3)
			members->setBinaryInterruptState(BinaryInterruptState::manybody);

		// members->CurrentTimeIrr		= ptclCM->CurrentTimeIrr;
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;
		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr;
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		members->NewCurrentBlockIrr	= ptclCM->NewCurrentBlockIrr;

		/* newly added by EW 2025.7.18 */
		members->TimeStepIrr     = ptclCM->TimeStepIrr;
		members->TimeBlockIrr    = ptclCM->TimeBlockIrr;
		members->TimeLevelIrr    = ptclCM->TimeLevelIrr;

		members->TimeStepReg     = ptclCM->TimeStepReg;
		members->TimeBlockReg    = ptclCM->TimeBlockReg;
		members->TimeLevelReg    = ptclCM->TimeLevelReg;

		/* original */
		members->TimeLevelIrr		= ptclCM->TimeLevelIrr; // test by EW 2025.1.29
		members->TimeLevelReg		= ptclCM->TimeLevelReg;

		// members->RadiusOfNeighbor = ACRadius*ACRadius; // added by EW 2025.1.16
		members->RadiusOfNeighbor = ptclCM->RadiusOfNeighbor; // modified by EW 2025.1.30

		for (int dim = 0; dim < Dim; dim++) {
			for (int order = 0; order < HERMITE_ORDER; order++)
				members->a_reg[dim][order] = ptclCM->a_reg[dim][order];
		}
		members->NumberOfNeighbor = ptclCM->NumberOfNeighbor;
		std::memcpy(members->Neighbors, ptclCM->Neighbors, sizeof(int)*ptclCM->NumberOfNeighbor);
		for (int j = 0; j < ptclCM->NumberOfMember; j++) {
			Particle* members_members = &particles[ptclCM->Members[j]];
			if (members_members->Mass > 0.0 && members_members->PID != members->PID) 
				members->Neighbors[members->NumberOfNeighbor++] = members_members->ParticleIndex;
		}
		if (members->CurrentTimeIrr == ptclCM->CurrentTimeIrr)
			computeMemberAcceleration01(members);
		else {
			for (int dim = 0; dim < Dim; dim++) {
				for (int order = 0; order < 2; order++) {
					members->a_irr[dim][order] = ptclCM->a_irr[dim][order];
					members->a_tot[dim][order] = ptclCM->a_tot[dim][order];
				}
			}
		}
	}
	for (int i = 0; i < ptclCM->NumberOfMember; i++) {
		members = &particles[ptclCM->Members[i]];
		if (members->Mass < 0.0)
			continue;

		if (members->CurrentTimeIrr == ptclCM->CurrentTimeIrr)
			computeMemberAcceleration23(members);
		else {
			for (int dim = 0; dim < Dim; dim++) {
				for (int order = 2; order < 4; order++) {
					members->a_irr[dim][order] = ptclCM->a_irr[dim][order];
					members->a_tot[dim][order] = ptclCM->a_tot[dim][order];
				}
			}
		}
	}
	
	for (int i = 0; i < ptclCM->NumberOfMember; i++) {
		members = &particles[ptclCM->Members[i]];
		if (members->Mass < 0.0)
			continue;

		members->CurrentTimeIrr = ptclCM->CurrentTimeIrr;

		// Update members' remaining properties...
		/* // Currently, we use a_reg of ptcl, so we don't need to newly set TimeLevelReg, TimeStepReg, TimeBlockReg by EW 2025.7.18
		members->calculateTimeStepReg();
		if (members->TimeLevelReg <= ptclCM->TimeLevelReg-1 
				&& members->TimeBlockReg/2+members->CurrentBlockReg >= NextRegTimeBlock)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
			members->TimeLevelReg = ptclCM->TimeLevelReg-1;
		}
		else if  (members->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
			members->TimeLevelReg = ptclCM->TimeLevelReg+1;
		}
		else 
			members->TimeLevelReg = ptclCM->TimeLevelReg;
		
		members->TimeStepReg  = static_cast<double>(pow(2, members->TimeLevelReg));
		members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block));
		*/

		if (members->NumberOfNeighbor != 0) {
			members->calculateTimeStepIrr2();
			// members->calculateTimeStepIrr();

			if (ptclCM->NumberOfMember > 2) {
				members->TimeLevelIrr--;
				members->TimeStepIrr = static_cast<double>(pow(2, members->TimeLevelIrr));
				members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
			}
			members->NewCurrentBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
			members->NextBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
		} 
		else {
			members->TimeStepReg -= members->CurrentTimeIrr - members->CurrentTimeReg;

			members->TimeStepIrr = members->TimeStepReg;
			members->NewCurrentBlockIrr = members->CurrentBlockReg + members->TimeBlockReg;
			members->NextBlockIrr = members->CurrentBlockReg + members->TimeBlockReg;
			members->TimeBlockIrr = members->NextBlockIrr - members->CurrentBlockIrr;

			members->CurrentBlockReg = members->CurrentBlockIrr;
			// I think fourth order correction is already done in computeAccelerationIrr function! by EW 2025.7.18
			// members->correctParticleFourthOrder(members->CurrentTimeIrr - members->CurrentTimeReg, members->Position, members->Velocity, members->a_tot);
			// members->updateParticle();
			members->CurrentTimeReg = members->CurrentTimeIrr;
		}

		members->printParticleInfo(binout);
	}
	ptclCM->clear();

	fflush(binout);
}

void computeMemberAcceleration01(Particle* members) {

	double new_time = members->CurrentTimeIrr; // the time to be advanced to

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[Dim], adot_tmp[Dim]; // 0 for current and 1 for predicted accelerations
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double m_r3;
	Particle* ptcl;

	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
		members->a_irr[dim][0] = 0.0;
		members->a_irr[dim][1] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion 01 Calculation
	 ********************************************************/

	for (int i = 0; i < members->NumberOfNeighbor; i++) {

		ptcl = &particles[members->Neighbors[i]];

		if (!ptcl->isActive) {
			if (ptcl->CMPtclIndex != -1)
				CMPtclsSet.insert(ptcl->CMPtclIndex);
			continue;
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - members->Position[dim];
			v[dim] = vel_neighbor[dim] - members->Velocity[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->Mass/(r2*sqrt(r2));

		for (int dim=0; dim<Dim; dim++) {
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	}

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->isActive) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", members->PID, ptcl->PID);
			assert(ptcl->isActive);
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - members->Position[dim];
			v[dim] = vel_neighbor[dim] - members->Velocity[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->Mass/(r2*sqrt(r2));

		for (int dim=0; dim<Dim; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	}

	double dt_ex = (members->NumberOfNeighbor == 0) ? 0.0 : (new_time - members->CurrentTimeReg)*EnzoTimeStep;
	for (int dim=0; dim<Dim; dim++) {
		members->a_irr[dim][0] = a_tmp[dim];
		members->a_irr[dim][1] = adot_tmp[dim];
		members->a_tot[dim][0] = members->a_reg[dim][0] + members->a_irr[dim][0] + members->a_reg[dim][1]*dt_ex; // affect the next
		members->a_tot[dim][1] = members->a_reg[dim][1] + members->a_irr[dim][1];
	}
}

void computeMemberAcceleration23(Particle* members) {

	double new_time = members->CurrentTimeIrr; // the time to be advanced to

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double m_r3;
	double a21[Dim], a21dot[Dim], a1[Dim], a2[Dim], a1dot[Dim], a2dot[Dim];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r3;
	double adot2, adot3;
	Particle* ptcl;

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion 23 Calculation
	 ********************************************************/

	 for (int dim=0; dim<Dim; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = members->a_tot[dim][0];
		a1dot[dim]  = members->a_tot[dim][1];
		members->a_irr[dim][2] = 0.0;
		members->a_irr[dim][3] = 0.0;
	}

	for (int i = 0; i < members->NumberOfNeighbor; i++) {

		ptcl = &particles[members->Neighbors[i]];

		if (!ptcl->isActive) {
			if (ptcl->CMPtclIndex != -1)
				CMPtclsSet.insert(ptcl->CMPtclIndex);
			continue;
		}

		r2 = 0;
		v2 = 0;
		vx = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		for (int dim=0; dim<Dim; dim++) {

			double dt = (new_time - ptcl->CurrentTimeIrr)*EnzoTimeStep;
			a2[dim] = (dt == 0) ? ptcl->a_tot[dim][0] : ptcl->a_tot[dim][1]*dt + ptcl->a_tot[dim][0];
			a2dot[dim] = ptcl->a_tot[dim][1];

			x[dim]     = pos_neighbor[dim] - members->Position[dim];
			v[dim]     = vel_neighbor[dim] - members->Velocity[dim];
			r2        += x[dim]*x[dim];
			vx        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl->Mass/r3;

		for (int dim=0; dim<Dim; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vx/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vx/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);

		for (int dim=0; dim<Dim; dim++) {
			adot2 = -ptcl->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
			adot3 = -ptcl->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
			members->a_irr[dim][2] += adot2;
			members->a_irr[dim][3] += adot3;
		}
	}

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->isActive) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", members->PID, ptcl->PID);
			assert(ptcl->isActive);
		}

		r2 = 0;
		v2 = 0;
		vx = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {

			double dt = (new_time - ptcl->CurrentTimeIrr)*EnzoTimeStep;
			a2[dim] = (dt == 0) ? ptcl->a_tot[dim][0] : ptcl->a_tot[dim][1]*dt + ptcl->a_tot[dim][0];
			a2dot[dim] = ptcl->a_tot[dim][1];

			x[dim]     = pos_neighbor[dim] - members->Position[dim];
			v[dim]     = vel_neighbor[dim] - members->Velocity[dim];
			r2        += x[dim]*x[dim];
			vx        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl->Mass/r3; 

		for (int dim=0; dim<Dim; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vx/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vx/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);

		for (int dim=0; dim<Dim; dim++) {
			adot2 = -ptcl->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
			adot3 = -ptcl->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
			members->a_irr[dim][2] += adot2;
			members->a_irr[dim][3] += adot3;
		}
	}

	for (int dim=0; dim<Dim; dim++)	 {
		members->a_tot[dim][2] = members->a_reg[dim][2] + members->a_irr[dim][2];
		members->a_tot[dim][3] = members->a_reg[dim][3] + members->a_irr[dim][3];
	}
}

// Calculate Position, Velocity, a_irr, a_tot at new_time
void computeMemberAccelerationIrr(Particle* members, double new_time) {

	double dt = (new_time - members->CurrentTimeIrr) * EnzoTimeStep;

	double x[Dim], v[Dim]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[Dim], adot_tmp[Dim]; // 0 for current and 1 for predicted accelerations
	double pos[Dim], vel[Dim];
	double pos_neighbor[Dim], vel_neighbor[Dim];
	double m_r3;
	Particle* ptcl;


	// This case should be treated very carefully...
	if (members->NumberOfNeighbor == 0) {
		members->predictParticleSecondOrder(new_time - members->CurrentTimeIrr, pos, vel);
		for (int dim=0; dim<Dim; dim++) {
			members->Position[dim] = pos[dim];
			members->Velocity[dim] = vel[dim];
		}
		return;
	}


	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<Dim; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion Calculation
	 ********************************************************/
	members->predictParticleSecondOrder(new_time - members->CurrentTimeIrr, pos, vel);

	for (int i=0; i<members->NumberOfNeighbor; i++) {

		ptcl = &particles[members->Neighbors[i]];

		if (!ptcl->isActive) {
			if (ptcl->CMPtclIndex != -1) {
				CMPtclsSet.insert(ptcl->CMPtclIndex);
			}
			continue;
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->Mass/(r2*sqrt(r2));

		for (int dim=0; dim<Dim; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	} // endfor ptcl

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->isActive) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", members->PID, ptcl->PID);
			assert(ptcl->isActive);
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->Mass/(r2*sqrt(r2));

		for (int dim=0; dim<Dim; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	}


	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;
	double dt_ex = (new_time - members->CurrentTimeReg)*EnzoTimeStep;

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
		this->Position[dim] = pos[dim] + a2*dt + a3*dt;
		this->Velocity[dim] = vel[dim] + 4*a2  + 5*a3;


		// note that these higher order terms and lowers have different neighbors
		this->a_irr[dim][0] = a_tmp[dim];
		this->a_irr[dim][1] = adot_tmp[dim];
		this->a_irr[dim][2] = a2*24/dt3;
		this->a_irr[dim][3] = a3*120/dt4;
#else
		// do the higher order correcteion
		da_dt2  = (members->a_irr[dim][0] - a_tmp[dim]) / dt2; 
		adot_dt = (members->a_irr[dim][1] + adot_tmp[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*members->a_irr[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		// 4th order correction
		// save the values in the temporary variables
		members->Position[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		members->Velocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		// note that these higher order terms and lowers have different neighbors
		members->a_irr[dim][0] = a_tmp[dim];
		members->a_irr[dim][1] = adot_tmp[dim];
		members->a_irr[dim][2] = a2;
		members->a_irr[dim][3] = a3;
#endif
	}

	for (int dim=0; dim<Dim; dim++) {
		members->a_tot[dim][0] = members->a_reg[dim][0] + members->a_irr[dim][0] + members->a_reg[dim][1]*dt_ex; // affect the next
		members->a_tot[dim][1] = members->a_reg[dim][1] + members->a_irr[dim][1];
		members->a_tot[dim][2] = members->a_reg[dim][2] + members->a_irr[dim][2];
		members->a_tot[dim][3] = members->a_reg[dim][3] + members->a_irr[dim][3];
	}
	members->CurrentTimeIrr = new_time;
}
#endif