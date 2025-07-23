#ifdef FEWBODY
#include "../global.h"
#include "unordered_set"

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void computeCMAcceleration(Particle* ptclCM);

void deleteGroup(Particle* ptclCM) {

	Group* ptclGroup = ptclCM->GroupInfo;

	ptclCM->NewNumberOfNeighbor = ptclGroup->sym_int.particles.getSize();

	assert(!ptclGroup->sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.4

	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		ptclCM->NewNeighbors[i] = members->ParticleIndex;

		for (int dim=0; dim<Dim; dim++) {
			particles[members->ParticleIndex].Position[dim] = ptclCM->Position[dim] + members->Position[dim];
			particles[members->ParticleIndex].Velocity[dim] = ptclCM->Velocity[dim] + members->Velocity[dim];
		}
		particles[members->ParticleIndex].Mass = members->Mass;
	}

/* // original version; this might take so much time by EW 2025.1.4
	for (int dim=0; dim<Dim; dim++) {
		ptclGroup->sym_int.particles.cm.Position[dim] = ptclCM->Position[dim];
		ptclGroup->sym_int.particles.cm.Velocity[dim] = ptclCM->Velocity[dim];
	}
	ptclGroup->sym_int.particles.shiftToOriginFrame();
	group2->sym_int.particles.template writeBackMemberAll<Particle>();
*/

	delete ptclGroup;
}

void Group::initialManager() {

	manager.interaction.gravitational_constant = 1.0;
	manager.time_step_min = 1e-13; // minimum physical time step // 1e-13 in ar.cxx
	manager.ds_scale = 1.0; // step size scaling factor // reference: ar.cxx
	manager.time_error_max = 0.25*manager.time_step_min; // time synchronization absolute error limit for AR, default is 0.25*dt-min
	// reference: ar.cxx
	manager.energy_error_relative_max = 1e-10; // relative energy error limit for AR, phase error requirement
	// 1e-10 in ar.cxx
	// 1e-8 in PeTar
	manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX; // maximum timescale for maximum slowdown factor, time-end
	// if (slowdown_timescale_max.value>0.0) manager.slowdown_timescale_max = slowdown_timescale_max.value;
	// else if (time_end.value>0.0) manager.slowdown_timescale_max = time_end.value;
	// else manager.slowdown_timescale_max = NUMERIC_FLOAT_MAX;
	// should be positive
	manager.slowdown_pert_ratio_ref = 1e-6; // slowdown perturbation ratio reference
	// 1e-6 in ar.cxx
	// 1e-4 in PeTar
	manager.step_count_max = 1000000; // number of maximum (integrate/output) step for AR integration // set symplectic order
	// 1000000 in PeTar & ar.cxx
	manager.step.initialSymplecticCofficients(-6); // Symplectic integrator order, should be even number
	// -6 in PeTar & ar.cxx
	manager.interrupt_detection_option = 2; // modify orbit or check interruption using modifyAndInterruptIter function
											// 0: turn off
											// 1: modify the binary orbits based on detetion criterion
											// 2. modify and also interrupt integrations

	// Eunwoo: it is turned off now but I will turn it on later.
	// Eunwoo: It can be used for merging star (dr < sum of radius) or destroy.
}

void Group::initialIntegrator(int NumMembers) {

	groupCM->NumberOfMember = 0;

	sym_int.manager = &manager;

	sym_int.particles.setMode(COMM::ListMode::copy);
    sym_int.particles.reserveMem(NumMembers);
	sym_int.info.reserveMem(NumMembers);

	fprintf(workerout, "Mem PID:");
    for (int i = 0; i < groupCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[groupCM->NewNeighbors[i]];
		if (!members->isCMptcl) {
			members->CMPtclIndex = groupCM->ParticleIndex; // added for write_out_group function by EW 2025.1.6
			sym_int.particles.addMemberAndAddress(*members);
			fprintf(workerout, " %d", sym_int.particles[groupCM->NumberOfMember].PID);
			groupCM->Members[groupCM->NumberOfMember++] = members->ParticleIndex;
		}
		else {
			for (int j=0; j < members->NewNumberOfNeighbor; j++) {
				Particle* members_members = &particles[members->NewNeighbors[j]];
				members_members->CMPtclIndex = groupCM->ParticleIndex; // added for write_out_group function by EW 2025.1.6
				sym_int.particles.addMemberAndAddress(*members_members);
				fprintf(workerout, " %d", sym_int.particles[groupCM->NumberOfMember].PID);
				groupCM->Members[groupCM->NumberOfMember++] = members_members->ParticleIndex;
			}
		}
    }
	fprintf(workerout, "\n");
	fflush(workerout);

	sym_int.info.r_break_crit = RSEARCH/position_unit; // distance criterion for checking stability
	// more information in symplectic_integrator.h
	// ar.cxx: 1e-3 pc
	// check whether the system is stable for 10000 out period and the apo-center is below break criterion
	// PeTar (hard.hpp): sym_int.info.r_break_crit = std::max(sym_int.info.r_break_crit,ptcl_origin[i].getRGroup());

	// // manager.print(std::cerr); // Eunwoo deleted

    sym_int.reserveIntegratorMem();
	sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);

	sym_int.particles.calcCenterOfMass();

	//! Fix step options for integration with adjusted step (not for time sychronizatio phase)
	// PeTar doesn't set this value explicitly!
	// sym_int.info.fix_step_option = AR::FixStepOption::none; // none: don't fix step
	// sym_int.info.fix_step_option = AR::FixStepOption::always; // always: use the given step without change
	// sym_int.info.fix_step_option = AR::FixStepOption::later; // later: fix step after a few adjustment of initial steps due to energy error

}


// Initialize new Few body group
void NewFBInitialization(Particle* ptclCM) {

	ptclCM->isActive = true;
	ptclCM->isCMptcl = true;

	Group* ptclGroup = new Group();

	ptclCM->GroupInfo = ptclGroup;
	ptclGroup->groupCM = ptclCM;

	// Find member particle with the biggest CurrentTimeIrr
	Particle* ptcl = &particles[ptclCM->NewNeighbors[0]];
	int NumberOfMembers=0;

	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];
		members->isActive = false;
		if (members->CurrentTimeIrr > ptcl->CurrentTimeIrr) {
        	ptcl = members;
    	}
		if (!members->isCMptcl)
			NumberOfMembers++;
		else
			NumberOfMembers += members->NewNumberOfNeighbor;
    }

	fprintf(workerout, "NewFBInitialization. CurrentTimeIrr (Myr): %e\n", ptcl->CurrentTimeIrr*EnzoTimeStep*1e4);

	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];

		double dt = ptcl->CurrentTimeIrr - members->CurrentTimeIrr;
		double pos[Dim], vel[Dim];
		members->predictParticleSecondOrder(dt, pos, vel);
		members->CurrentTimeIrr = ptcl->CurrentTimeIrr;
		
		if (members->isCMptcl) {
			Particle* members_members;

			for (int j=0; j<members->NewNumberOfNeighbor; j++) {
				members_members = &particles[members->NewNeighbors[j]];
				members_members->CurrentTimeIrr = ptcl->CurrentTimeIrr;

				for (int dim=0; dim<Dim; dim++) {
					members_members->Position[dim] += pos[dim] - members->Position[dim];
					members_members->Velocity[dim] += vel[dim] - members->Velocity[dim];
				}
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				members->Position[dim] = pos[dim];
				members->Velocity[dim] = vel[dim];
			}
		}
	}

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(NumberOfMembers); // Binary tree is made and CM particle is made automatically.

	for (int dim=0; dim<Dim; dim++) {
		ptclCM->Position[dim] = ptclGroup->sym_int.particles.cm.Position[dim];
		ptclCM->Velocity[dim] = ptclGroup->sym_int.particles.cm.Velocity[dim];
		ptclCM->Mass = ptclGroup->sym_int.particles.cm.Mass;
	}

	// Set ptcl information like time, PID, etc.
	ptclCM->RadiusOfNeighbor = ptcl->RadiusOfNeighbor; // original by EW 2025.2.4
	// ptclCM->RadiusOfNeighbor = ACRadius*ACRadius;

	ptclCM->CurrentTimeIrr  = ptcl->CurrentTimeIrr;
	ptclCM->CurrentTimeReg  = ptcl->CurrentTimeReg;
	ptclCM->CurrentBlockIrr = ptcl->CurrentBlockIrr; 
	ptclCM->CurrentBlockReg = ptcl->CurrentBlockReg;
	ptclCM->NewCurrentBlockIrr = ptcl->NewCurrentBlockIrr;

	ptclCM->TimeStepIrr     = ptcl->TimeStepIrr;
	ptclCM->TimeBlockIrr    = ptcl->TimeBlockIrr;
	ptclCM->TimeLevelIrr    = ptcl->TimeLevelIrr;

	ptclCM->TimeStepReg     = ptcl->TimeStepReg;
	ptclCM->TimeBlockReg    = ptcl->TimeBlockReg;
	ptclCM->TimeLevelReg    = ptcl->TimeLevelReg;

	for (int dim = 0; dim < Dim; dim++) {
		for (int order = 0; order < HERMITE_ORDER; order++) 
			ptclCM->a_reg[dim][order] = ptcl->a_reg[dim][order];
	}
	ptclCM->NumberOfNeighbor = ptcl->NumberOfNeighbor;
	std::memcpy(ptclCM->Neighbors, ptcl->Neighbors, sizeof(int)*ptcl->NumberOfNeighbor); // this will be adjusted soon! we should delete members...
	computeCMAcceleration(ptclCM); // neighbors are adjusted here!

	fprintf(workerout, "The ID of CM is %d.\n",ptclCM->PID);

	fprintf(workerout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		members->printParticleInfo(workerout);
    }

	ptclGroup->CurrentTime	= ptclCM->CurrentTimeIrr;

	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];

		if (members->isCMptcl)
			members->clear();
	}

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	// CalculateAcceleration01(ptclCM);
	// CalculateAcceleration23(ptclCM);

	for (int dim=0; dim<Dim; dim++) {
        for (int j=0; j<HERMITE_ORDER; j++)
            ptclGroup->sym_int.particles.cm.a_tot[dim][j] = ptclCM->a_tot[dim][j];
    }
	ptclGroup->sym_int.particles.cm.PID = ptclCM->PID; // added for ar_interaction.hpp by EW 2025.7.19
    
    ptclGroup->sym_int.particles.cm.NumberOfNeighbor = ptclCM->NumberOfNeighbor;
	std::memcpy(ptclGroup->sym_int.particles.cm.Neighbors, ptclCM->Neighbors, sizeof(int)*ptclCM->NumberOfNeighbor);

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	/* // Currently, we use a_reg of ptcl, so we don't need to newly set TimeLevelReg, TimeStepReg, TimeBlockReg by EW 2025.7.18
	ptclCM->calculateTimeStepReg();
	if (ptclCM->TimeLevelReg <= ptcl->TimeLevelReg-1 
			&& ptcl->TimeBlockReg/2+ptcl->CurrentBlockReg >= global_variable->NextRegTimeBlock)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg-1;
	}
	else if  (ptclCM->TimeLevelReg >= ptcl->TimeLevelReg+1) {
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg+1;
	}
	else 
		ptclCM->TimeLevelReg = ptcl->TimeLevelReg;

	ptclCM->TimeStepReg  = static_cast<double>(pow(2, ptclCM->TimeLevelReg));
	ptclCM->TimeBlockReg = static_cast<ULL>(pow(2, ptclCM->TimeLevelReg-time_block));
	*/

	if (ptclCM->NumberOfNeighbor != 0) {	

		// ptclCM->calculateTimeStepIrr();
		ptclCM->calculateTimeStepIrr2(); // by EW 2025.1.4
		ptclCM->NewCurrentBlockIrr = ptclCM->CurrentBlockIrr + ptclCM->TimeBlockIrr;
		ptclCM->NextBlockIrr = ptclCM->CurrentBlockIrr + ptclCM->TimeBlockIrr;
	}
	else {
		ptclCM->TimeStepReg -= ptclCM->CurrentTimeIrr - ptclCM->CurrentTimeReg;

		ptclCM->TimeStepIrr = ptclCM->TimeStepReg;
		ptclCM->NewCurrentBlockIrr = ptclCM->CurrentBlockReg + ptclCM->TimeBlockReg;
		ptclCM->NextBlockIrr = ptclCM->CurrentBlockReg + ptclCM->TimeBlockReg;
		ptclCM->TimeBlockIrr = ptclCM->NextBlockIrr - ptclCM->CurrentBlockIrr;

		ptclCM->CurrentBlockReg = ptclCM->CurrentBlockIrr;
		// I think fourth order correction is already done in computeAccelerationIrr function! by EW 2025.7.18
		// ptclCM->correctParticleFourthOrder(ptclCM->CurrentTimeIrr - ptclCM->CurrentTimeReg, ptclCM->Position, ptclCM->Velocity, ptclCM->a_tot);
		// ptclCM->updateParticle();
		ptclCM->CurrentTimeReg = ptclCM->CurrentTimeIrr;
	}

	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		// ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->RadiusOfNeighbor));
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, 1e-3/position_unit); // test12
		fprintf(workerout, "Bound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(workerout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
		fprintf(workerout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(workerout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(workerout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
	else {
		ptclGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
		fprintf(workerout, "Unbound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(workerout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(workerout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(workerout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}

	fprintf(workerout, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization\n");
	ptclCM->printParticleInfo(workerout);

	fprintf(workerout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(workerout);
}

// Use this function when many-body (>3) group breaks during SDAR integration.
void NewFBInitialization3(Group* group) {

	Group* ptclGroup = new Group();

	fprintf(workerout, "NewFBInitialization3. CurrentTimeIrr (Myr): %e\n", group->CurrentTime*EnzoTimeStep*1e4);

	Particle* ptclCM = group->groupCM;

	ptclGroup->groupCM = ptclCM;
	ptclGroup->isTerminate = group->isTerminate;
	ptclGroup->isMerger = group->isMerger;
	ptclGroup->CurrentTime = group->CurrentTime;

	ptclCM->NewNumberOfNeighbor = 0;
	for (int i = 0; i < ptclCM->NumberOfMember; i++) {
		Particle* members = &particles[ptclCM->Members[i]];
		if (members->Mass < 0)
			members->CMPtclIndex = -1;
		else {
			ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor] = ptclCM->Members[i];
			ptclCM->NewNumberOfNeighbor++;
		}
	}

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->NewNumberOfNeighbor); // Binary tree is made and CM particle is made automatically.
	ptclCM->NewNumberOfNeighbor = 0;
	/*
	After NewFBInitialization3, new binary forms and the same particles are detected as new binary members...
	I suspect this error happens because NewNumberOfNeighbor was not set to 0.
	Let's see what happens... by EW 2025.6.25
	It seems that this is right solution! by EW 2025.7.6
	*/

	delete group;
	ptclCM->GroupInfo = ptclGroup;

	fprintf(workerout, "The ID of CM is %d.\n", ptclCM->PID);

	fprintf(workerout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		members->printParticleInfo(workerout);
    }

	for (int dim=0; dim<Dim; dim++) {
        for (int j=0; j<HERMITE_ORDER; j++)
            ptclGroup->sym_int.particles.cm.a_tot[dim][j] = ptclCM->a_tot[dim][j];
    }
	ptclGroup->sym_int.particles.cm.PID = ptclCM->PID; // added for ar_interaction.hpp by EW 2025.7.19
    
    ptclGroup->sym_int.particles.cm.NumberOfNeighbor = ptclCM->NumberOfNeighbor;
	std::memcpy(ptclGroup->sym_int.particles.cm.Neighbors, ptclCM->Neighbors, sizeof(int)*ptclCM->NumberOfNeighbor);

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		// ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->RadiusOfNeighbor));
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, 1e-3/position_unit); // test12
		fprintf(workerout, "Bound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(workerout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
		fprintf(workerout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(workerout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(workerout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
	else {
		ptclGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
		fprintf(workerout, "Unbound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
		fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(workerout, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(workerout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(workerout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}

	fprintf(workerout, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization3\n");
	ptclCM->printParticleInfo(workerout);

	fprintf(workerout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(workerout);
}

void computeCMAcceleration(Particle* ptclCM) {

	double new_time = ptclCM->CurrentTimeIrr; // the time to be advanced to

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
		ptclCM->a_irr[dim][0] = 0.0;
		ptclCM->a_irr[dim][1] = 0.0;
		ptclCM->a_irr[dim][2] = 0.0;
		ptclCM->a_irr[dim][3] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion 01 Calculation
	 ********************************************************/

	for (int i = 0; i < ptclCM->NumberOfNeighbor; i++) {

		ptcl = &particles[ptclCM->Neighbors[i]];

		if (!ptcl->isActive) {
			if (ptcl->CMPtclIndex != -1) {
				if (ptcl->CMPtclIndex == ptclCM->ParticleIndex) {
					if (i != ptclCM->NumberOfNeighbor - 1)
						ptclCM->Neighbors[i] = ptclCM->Neighbors[ptclCM->NumberOfNeighbor - 1];
					ptclCM->NumberOfNeighbor--;
					i--;
				}
				else 
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
			x[dim] = pos_neighbor[dim] - ptclCM->Position[dim];
			v[dim] = vel_neighbor[dim] - ptclCM->Velocity[dim];

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
	assert(ptclCM->NumberOfNeighbor >= 0);

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->isActive) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", ptclCM->PID, ptcl->PID);
			assert(ptcl->isActive);
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predictParticleSecondOrder(new_time - ptcl->CurrentTimeIrr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<Dim; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - ptclCM->Position[dim];
			v[dim] = vel_neighbor[dim] - ptclCM->Velocity[dim];

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

	double dt_ex = (ptclCM->NumberOfNeighbor == 0) ? 0.0 : (new_time - ptclCM->CurrentTimeReg)*EnzoTimeStep;
	for (int dim=0; dim<Dim; dim++) {
		ptclCM->a_irr[dim][0] = a_tmp[dim];
		ptclCM->a_irr[dim][1] = adot_tmp[dim];
		ptclCM->a_tot[dim][0] = ptclCM->a_reg[dim][0] + ptclCM->a_irr[dim][0] + ptclCM->a_reg[dim][1]*dt_ex; // affect the next
		ptclCM->a_tot[dim][1] = ptclCM->a_reg[dim][1] + ptclCM->a_irr[dim][1];
	}

	/*******************************************************
	 * Irregular Acceleartion 23 Calculation
	 ********************************************************/

	double a21[Dim], a21dot[Dim], a1[Dim], a2[Dim], a1dot[Dim], a2dot[Dim];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r3;
	double adot2, adot3;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = ptclCM->a_tot[dim][0];
		a1dot[dim]  = ptclCM->a_tot[dim][1];
	}
 
	for (int i = 0; i < ptclCM->NumberOfNeighbor; i++) {

		ptcl = &particles[ptclCM->Neighbors[i]];

		if (!ptcl->isActive) // CMPtclsSet already contains all active CM particles
			continue;

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

			x[dim]     = pos_neighbor[dim] - ptclCM->Position[dim];
			v[dim]     = vel_neighbor[dim] - ptclCM->Velocity[dim];
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
			ptclCM->a_irr[dim][2] += adot2;
			ptclCM->a_irr[dim][3] += adot3;
		}
	}

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->isActive) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", ptclCM->PID, ptcl->PID);
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

			x[dim]     = pos_neighbor[dim] - ptclCM->Position[dim];
			v[dim]     = vel_neighbor[dim] - ptclCM->Velocity[dim];
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
			ptclCM->a_irr[dim][2] += adot2;
			ptclCM->a_irr[dim][3] += adot3;
		}
	}

	for (int dim=0; dim<Dim; dim++)	 {
		ptclCM->a_tot[dim][2] = ptclCM->a_reg[dim][2] + ptclCM->a_irr[dim][2];
		ptclCM->a_tot[dim][3] = ptclCM->a_reg[dim][3] + ptclCM->a_irr[dim][3];
	}
}
#endif
