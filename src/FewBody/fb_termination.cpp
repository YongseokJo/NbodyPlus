#ifdef FEWBODY
#include "../global.h"
#include <unordered_set>

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void computeMemberAcceleration01(Particle* members);
void computeMemberAcceleration23(Particle* members);
void computeMemberAccelerationIrr(Particle* members, double new_time);

void FBTermination(Particle* ptclCM) {

	fprintf(bin_output_file,"--------------------------------------\n");
	fprintf(bin_output_file,"In FBTermination.cpp... (CM PID: %d)\n", ptclCM->pid);
	fprintf(bin_output_file, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->current_time_irr*enzo_time_step*1e4);
	fprintf(bin_output_file, "CurrentTimeIrr of the first member (Myr): %e\n", particles[ptclCM->members[0]].current_time_irr*enzo_time_step*1e4);
	fprintf(bin_output_file, "N_member: %d\n", ptclCM->num_members);

	Particle* members;

	ptclCM->is_active = false;
	for (int i = 0; i < ptclCM->num_members; i++) {
		members = &particles[ptclCM->members[i]];
		members->cm_particle_index = -1;
		if (members->mass > 0.0)
			members->is_active = true;
	}

	num_particles--; // CM particle should be inactive by EW 2025.1.20

	// Set member information from CM ptcl; isActive = true & CMPtclIndex = -1 in RootRoutines.cpp
	for (int i = 0; i < ptclCM->num_members; i++) {
		members = &particles[ptclCM->members[i]];
		if (members->mass < 0.0)
			continue;
		
		num_particles++;
		/*
		if (ptclCM->num_members == 2)
			members->set_binary_interrupt_state(BinaryInterruptState::none);
		else if (ptclCM->num_members > 3)
			members->set_binary_interrupt_state(BinaryInterruptState::manybody);
		*/
		members->set_binary_interrupt_state(BinaryInterruptState::none);

		// members->current_time_irr		= ptclCM->current_time_irr;
		members->current_time_reg		= ptclCM->current_time_reg;
		members->current_block_irr	= ptclCM->current_block_irr;
		members->current_block_reg	= ptclCM->current_block_reg;
		members->new_current_block_irr	= ptclCM->current_block_irr; // modified by EW 2025.9.12

		/* newly added by EW 2025.7.18 */
		members->time_step_irr     = ptclCM->time_step_irr;
		members->time_block_irr    = ptclCM->time_block_irr;
		members->time_level_irr    = ptclCM->time_level_irr;

		members->time_step_reg     = ptclCM->time_step_reg;
		members->time_block_reg    = ptclCM->time_block_reg;
		members->time_level_reg    = ptclCM->time_level_reg;

		/* original */
		members->time_level_irr		= ptclCM->time_level_irr; // test by EW 2025.1.29
		members->time_level_reg		= ptclCM->time_level_reg;

		// members->neighbor_radius_sq = ACRadius*ACRadius; // added by EW 2025.1.16
		members->neighbor_radius_sq = ptclCM->neighbor_radius_sq; // modified by EW 2025.1.30

		for (int dim = 0; dim < DIM; dim++) {
			for (int order = 0; order < HERMITE_ORDER; order++)
				members->acc_regular[dim][order] = ptclCM->acc_regular[dim][order];
		}
		members->num_neighbors = ptclCM->num_neighbors;
		std::memcpy(neighbors + members->neighbors_offset, neighbors + ptclCM->neighbors_offset, sizeof(int)*ptclCM->num_neighbors);
		for (int j = 0; j < ptclCM->num_members; j++) {
			Particle* members_members = &particles[ptclCM->members[j]];
			if (members_members->mass > 0.0 && members_members->pid != members->pid) {
				neighbors[members->neighbors_offset + members->num_neighbors] = members_members->particle_index;
				members->num_neighbors++;
			}
		}
		if (members->current_time_irr == ptclCM->current_time_irr)
			computeMemberAcceleration01(members);
		else {
			for (int dim = 0; dim < DIM; dim++) {
				for (int order = 0; order < 2; order++) {
					members->acc_irregular[dim][order] = ptclCM->acc_irregular[dim][order];
					members->acc_total[dim][order] = ptclCM->acc_total[dim][order];
				}
			}
		}
	}
	for (int i = 0; i < ptclCM->num_members; i++) {
		members = &particles[ptclCM->members[i]];
		if (members->mass < 0.0)
			continue;

		if (members->current_time_irr == ptclCM->current_time_irr)
			computeMemberAcceleration23(members);
		else {
			for (int dim = 0; dim < DIM; dim++) {
				for (int order = 2; order < 4; order++) {
					members->acc_irregular[dim][order] = ptclCM->acc_irregular[dim][order];
					members->acc_total[dim][order] = ptclCM->acc_total[dim][order];
				}
			}
		}
	}
	
	for (int i = 0; i < ptclCM->num_members; i++) {
		members = &particles[ptclCM->members[i]];
		if (members->mass < 0.0)
			continue;

		members->current_time_irr = ptclCM->current_time_irr;

		// Update members' remaining properties...
		/* // Currently, we use a_reg of ptcl, so we don't need to newly set TimeLevelReg, TimeStepReg, TimeBlockReg by EW 2025.7.18
		members->calculate_time_step_reg();
		if (members->time_level_reg <= ptclCM->time_level_reg-1 
				&& members->time_block_reg/2+members->current_block_reg >= next_reg_time_block)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
			members->time_level_reg = ptclCM->time_level_reg-1;
		}
		else if  (members->time_level_reg >= ptclCM->time_level_reg+1) {
			members->time_level_reg = ptclCM->time_level_reg+1;
		}
		else 
			members->time_level_reg = ptclCM->time_level_reg;
		
		members->time_step_reg  = static_cast<double>(pow(2, members->time_level_reg));
		members->time_block_reg = static_cast<ull_t>(pow(2, members->time_level_reg-time_block));
		*/

		if (members->num_neighbors != 0) {
			members->calculate_time_step_irr_v2();
			// members->calculate_time_step_irr();
			/*
			if (ptclCM->num_members > 2) {
				members->time_level_irr--;
				members->time_step_irr = static_cast<double>(pow(2, members->time_level_irr));
				members->time_block_irr = static_cast<ull_t>(pow(2, members->time_level_irr-time_block));
			}
			*/
			// members->new_current_block_irr = members->current_block_irr + members->time_block_irr; // commented out by EW 2025.9.12
			members->next_block_irr = members->current_block_irr + members->time_block_irr;
		} 
		else {
			members->time_step_reg = (members->current_block_reg + members->time_block_reg - members->current_block_irr) * time_step;

			members->time_step_irr = members->time_step_reg;
			// members->new_current_block_irr = members->current_block_reg + members->time_block_reg; // commented out by EW 2025.9.12
			members->next_block_irr = members->current_block_reg + members->time_block_reg;
			members->time_block_irr = members->next_block_irr - members->current_block_irr;

			// members->current_block_reg = members->current_block_irr;
			// I think fourth order correction is already done in computeAccelerationIrr function! by EW 2025.7.18
			// members->correct_particle_fourth_order(members->current_time_irr - members->current_time_reg, members->position, members->velocity, members->acc_total);
			// members->update_particle();
			members->current_time_reg = members->current_time_irr;
		}

		members->print_particle_info(bin_output_file);
	}
	ptclCM->clear();

	fflush(bin_output_file);
}

void computeMemberAcceleration01(Particle* members) {

	double new_time = members->current_time_irr; // the time to be advanced to

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
		members->acc_irregular[dim][0] = 0.0;
		members->acc_irregular[dim][1] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion 01 Calculation
	 ********************************************************/

	for (int i = 0; i < members->num_neighbors; i++) {

		ptcl = &particles[neighbors[members->neighbors_offset + i]];

		if (!ptcl->is_active) {
			if (ptcl->cm_particle_index != -1)
				CMPtclsSet.insert(ptcl->cm_particle_index);
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
			x[dim] = pos_neighbor[dim] - members->position[dim];
			v[dim] = vel_neighbor[dim] - members->velocity[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->mass/(r2*sqrt(r2));

		for (int dim=0; dim<DIM; dim++) {
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	}

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->is_active) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", members->pid, ptcl->pid);
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
			x[dim] = pos_neighbor[dim] - members->position[dim];
			v[dim] = vel_neighbor[dim] - members->velocity[dim];

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

	double dt_ex = (members->num_neighbors == 0) ? 0.0 : (new_time - members->current_time_reg)*enzo_time_step;
	for (int dim=0; dim<DIM; dim++) {
		members->acc_irregular[dim][0] = a_tmp[dim];
		members->acc_irregular[dim][1] = adot_tmp[dim];
		members->acc_total[dim][0] = members->acc_regular[dim][0] + members->acc_irregular[dim][0] + members->acc_regular[dim][1]*dt_ex; // affect the next
		members->acc_total[dim][1] = members->acc_regular[dim][1] + members->acc_irregular[dim][1];
	}
}

void computeMemberAcceleration23(Particle* members) {

	double new_time = members->current_time_irr; // the time to be advanced to

	double x[DIM], v[DIM]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double pos_neighbor[DIM], vel_neighbor[DIM];
	double m_r3;
	double a21[DIM], a21dot[DIM], a1[DIM], a2[DIM], a1dot[DIM], a2dot[DIM];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r3;
	double adot2, adot3;
	Particle* ptcl;

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion 23 Calculation
	 ********************************************************/

	 for (int dim=0; dim<DIM; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = members->acc_total[dim][0];
		a1dot[dim]  = members->acc_total[dim][1];
		members->acc_irregular[dim][2] = 0.0;
		members->acc_irregular[dim][3] = 0.0;
	}

	for (int i = 0; i < members->num_neighbors; i++) {

		ptcl = &particles[neighbors[members->neighbors_offset + i]];

		if (!ptcl->is_active) {
			if (ptcl->cm_particle_index != -1)
				CMPtclsSet.insert(ptcl->cm_particle_index);
			continue;
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

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		for (int dim=0; dim<DIM; dim++) {

			double dt = (new_time - ptcl->current_time_irr)*enzo_time_step;
			a2[dim] = (dt == 0) ? ptcl->acc_total[dim][0] : ptcl->acc_total[dim][1]*dt + ptcl->acc_total[dim][0];
			a2dot[dim] = ptcl->acc_total[dim][1];

			x[dim]     = pos_neighbor[dim] - members->position[dim];
			v[dim]     = vel_neighbor[dim] - members->velocity[dim];
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
			members->acc_irregular[dim][2] += adot2;
			members->acc_irregular[dim][3] += adot3;
		}
	}

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->is_active) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", members->pid, ptcl->pid);
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

			x[dim]     = pos_neighbor[dim] - members->position[dim];
			v[dim]     = vel_neighbor[dim] - members->velocity[dim];
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
			members->acc_irregular[dim][2] += adot2;
			members->acc_irregular[dim][3] += adot3;
		}
	}

	for (int dim=0; dim<DIM; dim++)	 {
		members->acc_total[dim][2] = members->acc_regular[dim][2] + members->acc_irregular[dim][2];
		members->acc_total[dim][3] = members->acc_regular[dim][3] + members->acc_irregular[dim][3];
	}
}

// Calculate Position, Velocity, a_irr, a_tot at new_time
void computeMemberAccelerationIrr(Particle* members, double new_time) {

	double dt = (new_time - members->current_time_irr) * enzo_time_step;

	double x[DIM], v[DIM]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[DIM], adot_tmp[DIM]; // 0 for current and 1 for predicted accelerations
	double pos[DIM], vel[DIM];
	double pos_neighbor[DIM], vel_neighbor[DIM];
	double m_r3;
	Particle* ptcl;


	// This case should be treated very carefully...
	if (members->num_neighbors == 0) {
		members->predict_particle_second_order(new_time - members->current_time_irr, pos, vel);
		for (int dim=0; dim<DIM; dim++) {
			members->position[dim] = pos[dim];
			members->velocity[dim] = vel[dim];
		}
		return;
	}


	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<DIM; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion Calculation
	 ********************************************************/
	members->predict_particle_second_order(new_time - members->current_time_irr, pos, vel);

	for (int i=0; i<members->num_neighbors; i++) {

		ptcl = &particles[neighbors[members->neighbors_offset + i]];

		if (!ptcl->is_active) {
			if (ptcl->cm_particle_index != -1) {
				CMPtclsSet.insert(ptcl->cm_particle_index);
			}
			continue;
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predict_particle_second_order(new_time - ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<DIM; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
			vx += v[dim]*x[dim];
		}

		m_r3 = ptcl->mass/(r2*sqrt(r2));

		for (int dim=0; dim<DIM; dim++){
			a_tmp[dim]    += m_r3*x[dim];
			adot_tmp[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
		}
	} // endfor ptcl

	for (int i: CMPtclsSet) {
		ptcl = &particles[i];

		if (!ptcl->is_active) {
			fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", members->pid, ptcl->pid);
			assert(ptcl->is_active);
		}

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		ptcl->predict_particle_second_order(new_time - ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<DIM; dim++) {
			// calculate position and velocity differences for current time
			x[dim] = pos_neighbor[dim] - pos[dim];
			v[dim] = vel_neighbor[dim] - vel[dim];

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


	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;
	double dt_ex = (new_time - members->current_time_reg)*enzo_time_step;

	double A, B, C;

	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt; // VERY IMPORTANT BUG FIXED by EW 2025.7.18

	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	for (int dim=0; dim<DIM; dim++) {


#define noOptimization // I don't think it works 
#ifdef Optimization
		// do the higher order correcteion


		A = -12   *( this->acc_irregular[dim][0]   - a_tmp[dim] );
		B = -4*dt *( 2*this->acc_irregular[dim][1] + adot_tmp[dim] );
		C =  6*dt *( this->acc_irregular[dim][1]   + adot_tmp[dim] );

		a2 = dt *(A+B)/48;
		a3 = dt *(C-A)/120;

		//fprintf(stderr, "da_dt2=%.2e, adot_dt=%.2e, a2=%.2e, a3=%.2e\n", da_dt2, adot_dt, a2, a3);

		// 4th order correction
		// save the values in the temporary variables
		this->position[dim] = pos[dim] + a2*dt + a3*dt;
		this->velocity[dim] = vel[dim] + 4*a2  + 5*a3;


		// note that these higher order terms and lowers have different neighbors
		this->acc_irregular[dim][0] = a_tmp[dim];
		this->acc_irregular[dim][1] = adot_tmp[dim];
		this->acc_irregular[dim][2] = a2*24/dt3;
		this->acc_irregular[dim][3] = a3*120/dt4;
#else
		// do the higher order correcteion
		da_dt2  = (members->acc_irregular[dim][0] - a_tmp[dim]) / dt2; 
		adot_dt = (members->acc_irregular[dim][1] + adot_tmp[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*members->acc_irregular[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		// 4th order correction
		// save the values in the temporary variables
		members->position[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		members->velocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		// note that these higher order terms and lowers have different neighbors
		members->acc_irregular[dim][0] = a_tmp[dim];
		members->acc_irregular[dim][1] = adot_tmp[dim];
		members->acc_irregular[dim][2] = a2;
		members->acc_irregular[dim][3] = a3;
#endif
	}

	for (int dim=0; dim<DIM; dim++) {
		members->acc_total[dim][0] = members->acc_regular[dim][0] + members->acc_irregular[dim][0] + members->acc_regular[dim][1]*dt_ex; // affect the next
		members->acc_total[dim][1] = members->acc_regular[dim][1] + members->acc_irregular[dim][1];
		members->acc_total[dim][2] = members->acc_regular[dim][2] + members->acc_irregular[dim][2];
		members->acc_total[dim][3] = members->acc_regular[dim][3] + members->acc_irregular[dim][3];
	}
	members->current_time_irr = new_time;
}
#endif