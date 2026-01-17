#ifdef FEWBODY
#include "../global.h"
#include "unordered_set"

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void computeCMAcceleration(Particle* ptclCM);

void deleteGroup(Particle* ptclCM) {

	energy_binary -= ptclCM->group_info->sym_int.getEtot();
	energy_binary_sd -= ptclCM->group_info->sym_int.getEtotSlowDown();

	delete ptclCM->group_info;
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

	groupCM->num_members = 0;

	sym_int.manager = &manager;

	sym_int.particles.setMode(COMM::ListMode::copy);
    sym_int.particles.reserveMem(NumMembers);
	sym_int.info.reserveMem(NumMembers);

	fprintf(worker_output_file, "Mem PID:");
    for (int i = 0; i < groupCM->new_num_members; ++i) {
		Particle* members = &particles[groupCM->new_members[i]];
		if (!members->is_cm_particle) {
			members->cm_particle_index = groupCM->particle_index; // added for write_out_group function by EW 2025.1.6
			sym_int.particles.addMemberAndAddress(*members);
			fprintf(worker_output_file, " %d", sym_int.particles[groupCM->num_members].pid);
			groupCM->members[groupCM->num_members++] = members->particle_index;
		}
		else {
			for (int j=0; j < members->num_members; j++) {
				Particle* members_members = &particles[members->members[j]];
				members_members->cm_particle_index = groupCM->particle_index; // added for write_out_group function by EW 2025.1.6
				sym_int.particles.addMemberAndAddress(*members_members);
				fprintf(worker_output_file, " %d", sym_int.particles[groupCM->num_members].pid);
				groupCM->members[groupCM->num_members++] = members_members->particle_index;
			}
		}
    }
	fprintf(worker_output_file, "\n");
	fflush(worker_output_file);

	sym_int.info.r_break_crit = r_search; // distance criterion for checking stability (already in code units)
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

	ptclCM->is_active = true;
	ptclCM->is_cm_particle = true;

	Group* ptclGroup = new Group();

	ptclCM->group_info = ptclGroup;
	ptclGroup->groupCM = ptclCM;

	// Find member particle with the biggest CurrentTimeIrr
	Particle* ptcl = &particles[ptclCM->new_members[0]];
	int NumberOfMembers=0;

	for (int i = 0; i < ptclCM->new_num_members; ++i) {
		Particle* members = &particles[ptclCM->new_members[i]];
		members->is_active = false;
		if (members->current_time_irr > ptcl->current_time_irr) {
        	ptcl = members;
    	}
		if (!members->is_cm_particle)
			NumberOfMembers++;
		else
			NumberOfMembers += members->num_members;
    }

	fprintf(worker_output_file, "NewFBInitialization. CurrentTimeIrr (Myr): %e\n", ptcl->current_time_irr*enzo_time_step*1e4);

	for (int i = 0; i < ptclCM->new_num_members; ++i) {
		Particle* members = &particles[ptclCM->new_members[i]];

		double dt = ptcl->current_time_irr - members->current_time_irr;
		double pos[DIM], vel[DIM];
		members->predict_particle_second_order(dt, pos, vel);
		members->current_time_irr = ptcl->current_time_irr;
		
		if (members->is_cm_particle) {
			Particle* members_members;

			for (int j=0; j<members->num_members; j++) {
				members_members = &particles[members->members[j]];
				members_members->current_time_irr = ptcl->current_time_irr;

				for (int dim=0; dim<DIM; dim++) {
					members_members->position[dim] += pos[dim] - members->position[dim];
					members_members->velocity[dim] += vel[dim] - members->velocity[dim];
				}
			}
		}
		else {
			for (int dim=0; dim<DIM; dim++) {
				members->position[dim] = pos[dim];
				members->velocity[dim] = vel[dim];
			}
		}
	}

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(NumberOfMembers); // Binary tree is made and CM particle is made automatically.

	for (int dim=0; dim<DIM; dim++) {
		ptclCM->position[dim] = ptclGroup->sym_int.particles.cm.position[dim];
		ptclCM->velocity[dim] = ptclGroup->sym_int.particles.cm.velocity[dim];
		ptclCM->mass = ptclGroup->sym_int.particles.cm.mass;
	}

	// Set ptcl information like time, PID, etc.
	ptclCM->neighbor_radius_sq = ptcl->neighbor_radius_sq; // original by EW 2025.2.4
	// ptclCM->neighbor_radius_sq = ACRadius*ACRadius;

	ptclCM->current_time_irr  = ptcl->current_time_irr;
	ptclCM->current_time_reg  = ptcl->current_time_reg;
	ptclCM->current_block_irr = ptcl->current_block_irr; 
	ptclCM->current_block_reg = ptcl->current_block_reg;
	ptclCM->new_current_block_irr = ptcl->current_block_irr; // modified by EW 2025.9.12

	ptclCM->time_step_irr     = ptcl->time_step_irr;
	ptclCM->time_block_irr    = ptcl->time_block_irr;
	ptclCM->time_level_irr    = ptcl->time_level_irr;

	ptclCM->time_step_reg     = ptcl->time_step_reg;
	ptclCM->time_block_reg    = ptcl->time_block_reg;
	ptclCM->time_level_reg    = ptcl->time_level_reg;

	for (int dim = 0; dim < DIM; dim++) {
		for (int order = 0; order < HERMITE_ORDER; order++) 
			ptclCM->acc_regular[dim][order] = ptcl->acc_regular[dim][order];
	}
	ptclCM->num_neighbors = ptcl->num_neighbors;
	std::memcpy(neighbors + ptclCM->neighbors_offset, 
				neighbors + ptcl->neighbors_offset, 
				sizeof(int) * ptcl->num_neighbors); // this will be adjusted soon! we should delete members...
	computeCMAcceleration(ptclCM); // neighbors are adjusted here!

	fprintf(worker_output_file, "The ID of CM is %d.\n",ptclCM->pid);

	fprintf(worker_output_file, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		members->print_particle_info(worker_output_file);
    }

	ptclGroup->CurrentTime	= ptclCM->current_time_irr;

	for (int i = 0; i < ptclCM->new_num_members; ++i) {
		Particle* members = &particles[ptclCM->new_members[i]];

		if (members->is_cm_particle)
			members->clear();
	}

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	// CalculateAcceleration01(ptclCM);
	// CalculateAcceleration23(ptclCM);

	for (int dim=0; dim<DIM; dim++) {
        for (int j=0; j<HERMITE_ORDER; j++)
            ptclGroup->sym_int.particles.cm.acc_total[dim][j] = ptclCM->acc_total[dim][j];
    }
	ptclGroup->sym_int.particles.cm.pid = ptclCM->pid; // added for ar_interaction.hpp by EW 2025.7.19
	ptclGroup->sym_int.particles.cm.particle_index = ptclCM->particle_index; // added for separate shared neighbor array by EW 2025.9.1
	ptclGroup->sym_int.particles.cm.neighbors_offset = ptclCM->neighbors_offset; // added for separate shared neighbor array by EW 2025.9.1
    
    ptclGroup->sym_int.particles.cm.num_neighbors = ptclCM->num_neighbors;

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*enzo_time_step);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	energy_binary += ptclGroup->sym_int.getEtot();
	energy_binary_sd += ptclGroup->sym_int.getEtotSlowDown();

	/* // Currently, we use a_reg of ptcl, so we don't need to newly set TimeLevelReg, TimeStepReg, TimeBlockReg by EW 2025.7.18
	ptclCM->calculate_time_step_reg();
	if (ptclCM->time_level_reg <= ptcl->time_level_reg-1 
			&& ptcl->time_block_reg/2+ptcl->current_block_reg >= g_state->next_reg_time_block)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
		ptclCM->time_level_reg = ptcl->time_level_reg-1;
	}
	else if  (ptclCM->time_level_reg >= ptcl->time_level_reg+1) {
		ptclCM->time_level_reg = ptcl->time_level_reg+1;
	}
	else 
		ptclCM->time_level_reg = ptcl->time_level_reg;

	ptclCM->time_step_reg  = static_cast<double>(pow(2, ptclCM->time_level_reg));
	ptclCM->time_block_reg = static_cast<ull_t>(pow(2, ptclCM->time_level_reg-time_block));
	*/

	if (ptclCM->num_neighbors != 0) {	

		// ptclCM->calculate_time_step_irr();
		ptclCM->calculate_time_step_irr_v2(); // by EW 2025.1.4
		// ptclCM->new_current_block_irr = ptclCM->current_block_irr + ptclCM->time_block_irr; // commented out by EW 2025.9.12
		ptclCM->next_block_irr = ptclCM->current_block_irr + ptclCM->time_block_irr;
	}
	else {
		ptclCM->time_step_reg = (ptclCM->current_block_reg + ptclCM->time_block_reg - ptclCM->current_block_irr) * time_step;

		ptclCM->time_step_irr = ptclCM->time_step_reg;
		// ptclCM->new_current_block_irr = ptclCM->current_block_reg + ptclCM->time_block_reg; // commented out by EW 2025.9.12
		ptclCM->next_block_irr = ptclCM->current_block_reg + ptclCM->time_block_reg;
		ptclCM->time_block_irr = ptclCM->next_block_irr - ptclCM->current_block_irr;

		// ptclCM->current_block_reg = ptclCM->current_block_irr;
		// I think fourth order correction is already done in computeAccelerationIrr function! by EW 2025.7.18
		// ptclCM->correct_particle_fourth_order(ptclCM->current_time_irr - ptclCM->current_time_reg, ptclCM->position, ptclCM->velocity, ptclCM->acc_total);
		// ptclCM->update_particle();
		ptclCM->current_time_reg = ptclCM->current_time_irr;
	}

	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		// ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->neighbor_radius_sq));
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, 1e-3/position_unit); // test12
		fprintf(worker_output_file, "Bound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
		fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(worker_output_file, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
		fprintf(worker_output_file, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(worker_output_file, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(worker_output_file, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
	else {
		ptclGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
		fprintf(worker_output_file, "Unbound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
		fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(worker_output_file, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(worker_output_file, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(worker_output_file, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}

	fprintf(worker_output_file, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization\n");
	ptclCM->print_particle_info(worker_output_file);

	fprintf(worker_output_file, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(worker_output_file);
}

// Use this function when many-body (>3) group breaks during SDAR integration.
void NewFBInitialization3(Group* group) {

	Group* ptclGroup = new Group();

	fprintf(worker_output_file, "NewFBInitialization3. CurrentTimeIrr (Myr): %e\n", group->CurrentTime*enzo_time_step*1e4);

	Particle* ptclCM = group->groupCM;

	ptclGroup->groupCM = ptclCM;
	ptclGroup->isTerminate = group->isTerminate;
	ptclGroup->isMerger = group->isMerger;
	ptclGroup->CurrentTime = group->CurrentTime;

	ptclCM->new_num_members = 0;
	for (int i = 0; i < ptclCM->num_members; i++) {
		Particle* members = &particles[ptclCM->members[i]];
		if (members->mass < 0)
			members->cm_particle_index = -1;
		else {
			ptclCM->new_members[ptclCM->new_num_members++] = ptclCM->members[i];
		}
	}

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->new_num_members); // Binary tree is made and CM particle is made automatically.
	ptclCM->new_num_members = 0;
	/*
	After NewFBInitialization3, new binary forms and the same particles are detected as new binary members...
	I suspect this error happens because NewNumberOfMember was not set to 0.
	Let's see what happens... by EW 2025.6.25
	It seems that this is right solution! by EW 2025.7.6
	*/

	energy_binary -= group->sym_int.getEtot();
	energy_binary_sd -= group->sym_int.getEtotSlowDown();

	delete group;
	ptclCM->group_info = ptclGroup;

	energy_binary += ptclGroup->sym_int.getEtot();
	energy_binary_sd += ptclGroup->sym_int.getEtotSlowDown();

	fprintf(worker_output_file, "The ID of CM is %d.\n", ptclCM->pid);

	fprintf(worker_output_file, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		members->print_particle_info(worker_output_file);
    }

	for (int dim=0; dim<DIM; dim++) {
        for (int j=0; j<HERMITE_ORDER; j++)
            ptclGroup->sym_int.particles.cm.acc_total[dim][j] = ptclCM->acc_total[dim][j];
    }
	ptclGroup->sym_int.particles.cm.pid = ptclCM->pid; // added for ar_interaction.hpp by EW 2025.7.19
	ptclGroup->sym_int.particles.cm.particle_index = ptclCM->particle_index; // added for separate shared neighbor array by EW 2025.9.1
	ptclGroup->sym_int.particles.cm.neighbors_offset = ptclCM->neighbors_offset; // added for separate shared neighbor array by EW 2025.9.1
    
    ptclGroup->sym_int.particles.cm.num_neighbors = ptclCM->num_neighbors;

	ptclGroup->sym_int.initialIntegration(ptclGroup->CurrentTime*enzo_time_step);
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);

	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		// ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->neighbor_radius_sq));
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, 1e-3/position_unit); // test12
		fprintf(worker_output_file, "Bound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
		fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(worker_output_file, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
		fprintf(worker_output_file, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(worker_output_file, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(worker_output_file, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
	else {
		ptclGroup->sym_int.info.r_break_crit = 2*bin_root.semi*(1-bin_root.ecc); // r_break_crit = 2*peri
		fprintf(worker_output_file, "Unbound. separation: %e pc\n\t", bin_root.r*position_unit);
		fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
		fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
		fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
		fprintf(worker_output_file, "period: %e Myr\n\t", bin_root.period*1e4);
		fprintf(worker_output_file, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(worker_output_file, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}

	fprintf(worker_output_file, "\nFBInitialization.cpp: result of CM particle value calculation from function NewFBInitialization3\n");
	ptclCM->print_particle_info(worker_output_file);

	fprintf(worker_output_file, "---------------------END-OF-NEW-GROUP---------------------\n\n");
	fflush(worker_output_file);
}

void computeCMAcceleration(Particle* ptclCM) {

	double new_time = ptclCM->current_time_irr; // the time to be advanced to

	double x[DIM], v[DIM]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[DIM], adot_tmp[DIM]; // 0 for current and 1 for predicted accelerations
	double pos_neighbor[DIM], vel_neighbor[DIM];
	double m_r3;
	Particle* ptcl;

	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<DIM; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
		ptclCM->acc_irregular[dim][0] = 0.0;
		ptclCM->acc_irregular[dim][1] = 0.0;
		ptclCM->acc_irregular[dim][2] = 0.0;
		ptclCM->acc_irregular[dim][3] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion 01 Calculation
	 ********************************************************/

	for (int i = 0; i < ptclCM->num_neighbors; i++) {

		ptcl = &particles[neighbors[ptclCM->neighbors_offset + i]];

		if (!ptcl->is_active) {
			if (ptcl->cm_particle_index != -1) {
				if (ptcl->cm_particle_index == ptclCM->particle_index) {
					if (i != ptclCM->num_neighbors - 1)
						neighbors[ptclCM->neighbors_offset + i] = neighbors[ptclCM->neighbors_offset + ptclCM->num_neighbors - 1];
					ptclCM->num_neighbors--;
					i--;
				}
				else 
					CMPtclsSet.insert(ptcl->cm_particle_index);
			}
			continue;
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		if (ptcl->num_neighbors == 0)
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_reg, pos_neighbor, vel_neighbor);
		else
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<DIM; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - ptclCM->position[dim];
			v[dim] = vel_neighbor[dim] - ptclCM->velocity[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->mass/(r2*sqrt(r2));

		for (int dim=0; dim<DIM; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	}
	assert(ptclCM->num_neighbors >= 0);

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->is_active) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", ptclCM->pid, ptcl->pid);
			assert(ptcl->is_active);
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		if (ptcl->num_neighbors == 0)
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_reg, pos_neighbor, vel_neighbor);
		else
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<DIM; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - ptclCM->position[dim];
			v[dim] = vel_neighbor[dim] - ptclCM->velocity[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->mass/(r2*sqrt(r2));

		for (int dim=0; dim<DIM; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	}

	double dt_ex = (ptclCM->num_neighbors == 0) ? 0.0 : (new_time - ptclCM->current_time_reg)*enzo_time_step;
	for (int dim=0; dim<DIM; dim++) {
		ptclCM->acc_irregular[dim][0] = a_tmp[dim];
		ptclCM->acc_irregular[dim][1] = adot_tmp[dim];
		ptclCM->acc_total[dim][0] = ptclCM->acc_regular[dim][0] + ptclCM->acc_irregular[dim][0] + ptclCM->acc_regular[dim][1]*dt_ex; // affect the next
		ptclCM->acc_total[dim][1] = ptclCM->acc_regular[dim][1] + ptclCM->acc_irregular[dim][1];
	}

	/*******************************************************
	 * Irregular Acceleartion 23 Calculation
	 ********************************************************/

	double a21[DIM], a21dot[DIM], a1[DIM], a2[DIM], a1dot[DIM], a2dot[DIM];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r3;
	double adot2, adot3;

	for (int dim=0; dim<DIM; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = ptclCM->acc_total[dim][0];
		a1dot[dim]  = ptclCM->acc_total[dim][1];
	}
 
	for (int i = 0; i < ptclCM->num_neighbors; i++) {

		ptcl = &particles[neighbors[ptclCM->neighbors_offset + i]];

		if (!ptcl->is_active) // CMPtclsSet already contains all active CM particles
			continue;

		r2 = 0;
		v2 = 0;
		vx = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		if (ptcl->num_neighbors == 0)
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_reg, pos_neighbor, vel_neighbor);
		else
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		for (int dim=0; dim<DIM; dim++) {

			double dt = (new_time - ptcl->current_time_irr)*enzo_time_step;
			a2[dim] = (dt == 0) ? ptcl->acc_total[dim][0] : ptcl->acc_total[dim][1]*dt + ptcl->acc_total[dim][0];
			a2dot[dim] = ptcl->acc_total[dim][1];

			x[dim]     = pos_neighbor[dim] - ptclCM->position[dim];
			v[dim]     = vel_neighbor[dim] - ptclCM->velocity[dim];
			r2        += x[dim]*x[dim];
			vx        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl->mass/r3; 

		for (int dim=0; dim<DIM; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vx/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vx/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);

		for (int dim=0; dim<DIM; dim++) {
			adot2 = -ptcl->mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
			adot3 = -ptcl->mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
			ptclCM->acc_irregular[dim][2] += adot2;
			ptclCM->acc_irregular[dim][3] += adot3;
		}
	}

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->is_active) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", ptclCM->pid, ptcl->pid);
			assert(ptcl->is_active);
		}

		r2 = 0;
		v2 = 0;
		vx = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		if (ptcl->num_neighbors == 0)
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_reg, pos_neighbor, vel_neighbor);
		else
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<DIM; dim++) {

			double dt = (new_time - ptcl->current_time_irr)*enzo_time_step;
			a2[dim] = (dt == 0) ? ptcl->acc_total[dim][0] : ptcl->acc_total[dim][1]*dt + ptcl->acc_total[dim][0];
			a2dot[dim] = ptcl->acc_total[dim][1];

			x[dim]     = pos_neighbor[dim] - ptclCM->position[dim];
			v[dim]     = vel_neighbor[dim] - ptclCM->velocity[dim];
			r2        += x[dim]*x[dim];
			vx        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl->mass/r3; 

		for (int dim=0; dim<DIM; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vx/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vx/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);

		for (int dim=0; dim<DIM; dim++) {
			adot2 = -ptcl->mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
			adot3 = -ptcl->mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
			ptclCM->acc_irregular[dim][2] += adot2;
			ptclCM->acc_irregular[dim][3] += adot3;
		}
	}

	for (int dim=0; dim<DIM; dim++)	 {
		ptclCM->acc_total[dim][2] = ptclCM->acc_regular[dim][2] + ptclCM->acc_irregular[dim][2];
		ptclCM->acc_total[dim][3] = ptclCM->acc_regular[dim][3] + ptclCM->acc_irregular[dim][3];
	}
}
#endif
