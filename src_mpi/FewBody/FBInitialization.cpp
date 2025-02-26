#ifdef FEWBODY
#include "../global.h"

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);

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
	// ptclCM->RadiusOfNeighbor = ptcl->RadiusOfNeighbor; // original by EW 2025.2.4
	ptclCM->RadiusOfNeighbor = ACRadius*ACRadius;

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


	fprintf(workerout, "The ID of CM is %d.\n",ptclCM->PID);

	fprintf(workerout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		fprintf(workerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(workerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(workerout, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
		fprintf(workerout, "PID: %d. Particle type: %d\n", members->PID, members->ParticleType);
        fprintf(workerout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(workerout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(workerout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(workerout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(workerout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(workerout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(workerout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(workerout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(workerout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);

		// fprintf(workerout, "PID: %d. Position - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0], members->Position[1], members->Position[2]);
		// fprintf(workerout, "PID: %d. Velocity - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0], members->Velocity[1], members->Velocity[2]);
		// fprintf(workerout, "PID: %d. Mass - %e, \n", members->PID, members->Mass);
		// fprintf(workerout, "PID: %d. Particle type: %d\n", members->PID, members->ParticleType);
        // fprintf(workerout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(workerout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(workerout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(workerout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		// fprintf(workerout, "PID: %d. Time Steps - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep, members->TimeStepReg*EnzoTimeStep);
		// fprintf(workerout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(workerout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);
    }

	ptclGroup->CurrentTime	= ptclCM->CurrentTimeIrr;

	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];

		if (members->isCMptcl)
			members->clear();
	}

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	CalculateAcceleration01(ptclCM);
	CalculateAcceleration23(ptclCM);

	for (int dim=0; dim<Dim; dim++) {
        for (int j=0; j<HERMITE_ORDER; j++)
            ptclGroup->sym_int.particles.cm.a_irr[dim][j] = ptclCM->a_irr[dim][j];
    }
    
    ptclGroup->sym_int.particles.cm.NumberOfNeighbor = ptclCM->NumberOfNeighbor;
    for (int i=0; i<ptclCM->NumberOfNeighbor; i++)
    	ptclGroup->sym_int.particles.cm.Neighbors[i] = ptclCM->Neighbors[i];

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*EnzoTimeStep);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

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

		ptclCM->correctParticleFourthOrder(ptclCM->CurrentTimeIrr - ptclCM->CurrentTimeReg, ptclCM->Position, ptclCM->Velocity, ptclCM->a_tot);
		ptclCM->updateParticle();
		ptclCM->CurrentTimeReg = ptclCM->CurrentTimeIrr;
	}

/*
	while (ptclCM->CurrentBlockIrr+ptclCM->TimeBlockIrr <= global_time_irr 
			&& ptclCM->TimeLevelIrr <= ptcl->TimeLevelIrr) { //first condition guarantees that ptclcm is small than ptcl
		ptclCM->TimeLevelIrr++;
		ptclCM->TimeStepIrr  = static_cast<double>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}
*/

/* // Eunwoo added
	while (ptclCM->TimeStepIrr*EnzoTimeStep*1e4 < 1e-7) {
		ptclCM->TimeLevelIrr += 1;
		ptclCM->TimeStepIrr = static_cast<REAL>(pow(2, ptclCM->TimeLevelIrr));
		ptclCM->TimeBlockIrr = static_cast<ULL>(pow(2, ptclCM->TimeLevelIrr-time_block));
	}
*/

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

	fprintf(workerout, "Position (pc) - x:%e, y:%e, z:%e, \n", ptclCM->Position[0]*position_unit, ptclCM->Position[1]*position_unit, ptclCM->Position[2]*position_unit);
	fprintf(workerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[1]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[2]*velocity_unit/yr*pc/1e5);
	fprintf(workerout, "Mass (Msol) - %e, \n", ptclCM->Mass*mass_unit);
	fprintf(workerout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(workerout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(workerout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(workerout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	fprintf(workerout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_reg[0][0], ptclCM->a_reg[1][0], ptclCM->a_reg[2][0]);
	// fprintf(workerout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_reg[0][1], ptclCM->a_reg[1][1], ptclCM->a_reg[2][1]);
	// fprintf(workerout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_reg[0][2], ptclCM->a_reg[1][2], ptclCM->a_reg[2][2]);
	// fprintf(workerout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_reg[0][3], ptclCM->a_reg[1][3], ptclCM->a_reg[2][3]);
	fprintf(workerout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_irr[0][0], ptclCM->a_irr[1][0], ptclCM->a_irr[2][0]);
	// fprintf(workerout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_irr[0][1], ptclCM->a_irr[1][1], ptclCM->a_irr[2][1]);
	// fprintf(workerout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_irr[0][2], ptclCM->a_irr[1][2], ptclCM->a_irr[2][2]);
	// fprintf(workerout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_irr[0][3], ptclCM->a_irr[1][3], ptclCM->a_irr[2][3]);
	fprintf(workerout, "Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);
	// fprintf(workerout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclCM->TimeBlockIrr, ptclCM->TimeBlockReg);
	// fprintf(workerout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclCM->CurrentBlockIrr, ptclCM->CurrentBlockReg);

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
		if (members->Mass == 0)
			members->CMPtclIndex = -1;
		else {
			ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor] = ptclCM->Members[i];
			ptclCM->NewNumberOfNeighbor++;
		}
	}

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->NewNumberOfNeighbor); // Binary tree is made and CM particle is made automatically.

	delete group;
	ptclCM->GroupInfo = ptclGroup;

	fprintf(workerout, "The ID of CM is %d.\n", ptclCM->PID);

	fprintf(workerout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];

		fprintf(workerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(workerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(workerout, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
        fprintf(workerout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
        // fprintf(workerout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(workerout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(workerout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(workerout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(workerout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(workerout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(workerout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(workerout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(workerout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(workerout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);
    }

	for (int dim=0; dim<Dim; dim++) {
        for (int j=0; j<HERMITE_ORDER; j++)
            ptclGroup->sym_int.particles.cm.a_irr[dim][j] = ptclCM->a_irr[dim][j];
    }
    
    ptclGroup->sym_int.particles.cm.NumberOfNeighbor = ptclCM->NumberOfNeighbor;
    for (int i=0; i<ptclCM->NumberOfNeighbor; i++)
    	ptclGroup->sym_int.particles.cm.Neighbors[i] = ptclCM->Neighbors[i];

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

	fprintf(workerout, "Position (pc) - x:%e, y:%e, z:%e, \n", ptclCM->Position[0]*position_unit, ptclCM->Position[1]*position_unit, ptclCM->Position[2]*position_unit);
	fprintf(workerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[1]*velocity_unit/yr*pc/1e5, ptclCM->Velocity[2]*velocity_unit/yr*pc/1e5);
	fprintf(workerout, "Mass (Msol) - %e, \n", ptclCM->Mass*mass_unit);
	fprintf(workerout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
	// fprintf(workerout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
	// fprintf(workerout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
	// fprintf(workerout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
	fprintf(workerout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_reg[0][0], ptclCM->a_reg[1][0], ptclCM->a_reg[2][0]);
	// fprintf(workerout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_reg[0][1], ptclCM->a_reg[1][1], ptclCM->a_reg[2][1]);
	// fprintf(workerout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_reg[0][2], ptclCM->a_reg[1][2], ptclCM->a_reg[2][2]);
	// fprintf(workerout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_reg[0][3], ptclCM->a_reg[1][3], ptclCM->a_reg[2][3]);
	fprintf(workerout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_irr[0][0], ptclCM->a_irr[1][0], ptclCM->a_irr[2][0]);
	// fprintf(workerout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_irr[0][1], ptclCM->a_irr[1][1], ptclCM->a_irr[2][1]);
	// fprintf(workerout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_irr[0][2], ptclCM->a_irr[1][2], ptclCM->a_irr[2][2]);
	// fprintf(workerout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_irr[0][3], ptclCM->a_irr[1][3], ptclCM->a_irr[2][3]);
	fprintf(workerout, "Time Steps (Myr) - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr*EnzoTimeStep*1e4, ptclCM->TimeStepReg*EnzoTimeStep*1e4);
	// fprintf(workerout, "Time Blocks - irregular:%llu, regular:%llu \n", ptclCM->TimeBlockIrr, ptclCM->TimeBlockReg);
	// fprintf(workerout, "Current Blocks - irregular: %llu, regular:%llu \n", ptclCM->CurrentBlockIrr, ptclCM->CurrentBlockReg);

	fprintf(workerout, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(workerout);
}
#endif
