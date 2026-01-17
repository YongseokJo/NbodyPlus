#include "../global.h"

void CalculateAcceleration01(Particle* ptcl);
void CalculateAcceleration23(Particle* ptcl);



/*******************************************************
 *  Initialize Accelerations
 *******************************************************/



void CalculateAcceleration01(Particle* ptcl1) {

	double x[DIM], v[DIM];
	double m_r3;
	double v2;
	double r2 = 0;
	double vx = 0;

	for (int dim=0; dim<DIM; dim++) {
		x[dim]    = 0.;
		v[dim]    = 0.;
	}

	ptcl1->num_neighbors = 0;
	for(int dim=0; dim<DIM; dim++) {
		for (int order=0; order<HERMITE_ORDER; order++) {
			ptcl1->acc_regular[dim][order] = 0.0;
			ptcl1->acc_irregular[dim][order] = 0.0;
			ptcl1->acc_total[dim][order] = 0.0;
		}
	}

	//std::cout << "nbody+: Entering CalculateInitialAcceleration  ..." << std::endl;
	//ptcl1->predict_particle_second_order(newTime);
	
	//fprintf(stdout, "pid=%d, nn=%d, numpart=%d\n", ptcl1->pid, ptcl1->num_neighbors, last_particle_index);
	Particle *ptcl2;
	for (int i=0; i<=g_state->last_particle_index; i++) {
		ptcl2 = &particles[i];

		if (!ptcl2->is_active || ptcl1->pid == ptcl2->pid) {
			continue;
		}

		r2 = 0;
		vx = 0;
		v2 = 0;

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		//ptcl2->predict_particle_second_order(newTime);
		for (int dim=0; dim<DIM; dim++) {
			x[dim] = ptcl2->position[dim] - ptcl1->position[dim];
			v[dim] = ptcl2->velocity[dim] - ptcl1->velocity[dim];
			r2    += x[dim]*x[dim];
			vx    += v[dim]*x[dim];
			v2    += v[dim]*v[dim];
		}

		m_r3 = ptcl2->mass/r2/sqrt(r2); 

		if (r2 > ptcl1->neighbor_radius_sq) {
			for (int dim=0; dim<DIM; dim++) {
				// Calculate 0th and 1st derivatives of acceleration
				ptcl1->acc_regular[dim][0] += m_r3*x[dim];
				ptcl1->acc_regular[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
		else {
			for (int dim=0; dim<DIM; dim++) {
				ptcl1->acc_irregular[dim][0] += m_r3*x[dim];
				ptcl1->acc_irregular[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
			if (!ptcl2->is_cm_particle) {
				neighbors[ptcl1->neighbors_offset + ptcl1->num_neighbors] = ptcl2->particle_index;
				ptcl1->num_neighbors++;
				assert(ptcl1->num_neighbors < MAX_NUM_NEIGHBOR);
			}
			else {
				for (int j=0; j<ptcl2->num_members; j++) {
					neighbors[ptcl1->neighbors_offset + ptcl1->num_neighbors] = ptcl2->members[j];
					ptcl1->num_neighbors++;
					assert(ptcl1->num_neighbors < MAX_NUM_NEIGHBOR);
				}
			}
			//fprintf(stdout, "pid=%d, nn=%d\n", ptcl1->pid, ptcl1->num_neighbors);
		} // endfor dim
	} // endfor ptcl2

	for (int dim=0; dim<DIM; dim++)	 {
		for (int order=0; order<2; order++) {
			ptcl1->acc_total[dim][order] = ptcl1->acc_regular[dim][order] + ptcl1->acc_irregular[dim][order]; 
		}
	}
	return;
}

void CalculateAcceleration23(Particle* ptcl1) {

	double x[DIM], v[DIM], a21[DIM], a21dot[DIM], a1[DIM], a2[DIM], a1dot[DIM], a2dot[DIM];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r2, r3, vr, m_r3;
	double adot2, adot3;

	for (int dim=0; dim<DIM; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = ptcl1->acc_total[dim][0];
		a1dot[dim]  = ptcl1->acc_total[dim][1];
	}

	Particle *ptcl2;
	for (int i=0; i<=g_state->last_particle_index; i++) {
		ptcl2 = &particles[i];

		if (!ptcl2->is_active || ptcl1->pid == ptcl2->pid) {
			continue;
		}

		r2 = 0;
		r3 = 0;
		v2 = 0;
		vr = 0;
		rdf_r2 = 0;
		vdf_r2 = 0;
		rdfdot_r2 = 0;

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		for (int dim=0; dim<DIM; dim++) {
			a2[dim]    = ptcl2->acc_total[dim][0];
			a2dot[dim] = ptcl2->acc_total[dim][1];
			x[dim]     = ptcl2->position[dim] - ptcl1->position[dim];
			v[dim]     = ptcl2->velocity[dim] - ptcl1->velocity[dim];
			r2        += x[dim]*x[dim];
			vr        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl2->mass/r3; 

		for (int dim=0; dim<DIM; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vr/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vr/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);


		if (r2 > ptcl1->neighbor_radius_sq) {
			for (int dim=0; dim<DIM; dim++) {
				adot2 = -ptcl2->mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				ptcl1->acc_regular[dim][2] += adot2;
				ptcl1->acc_regular[dim][3] += adot3;
			}
		}
		else {
			for (int dim=0; dim<DIM; dim++) {
				adot2 = -ptcl2->mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				ptcl1->acc_irregular[dim][2] += adot2;
				ptcl1->acc_irregular[dim][3] += adot3;
			}
		} // endfor if
	} //endfor ptcl2

	for (int dim=0; dim<DIM; dim++)	 {
		for (int order=2; order<4; order++) {
			ptcl1->acc_total[dim][order] = ptcl1->acc_regular[dim][order] + ptcl1->acc_irregular[dim][order];
		}
	}
}




/*******************************************************
 *  Initialize Time Steps 
 *******************************************************/


double getNewTimeStep(double f[3][4], double df[3][4]);
void getBlockTimeStep(double dt, int& TimeLevel, ull_t &TimeBlock, double &TimeStep);


void Particle::initialize_time_step() {

	//std::cout << "Initializing timesteps ..." << std::endl;

	double dtIrr, dtReg;

	dtReg = getNewTimeStep(this->acc_regular, this->acc_regular);
	//std::cout << "dtReg=" << dtReg << std::endl;
	getBlockTimeStep(dtReg, this->time_level_reg, this->time_block_reg, this->time_step_reg);
	//std::cout << "TimeStepReg=" << this->time_step_reg*enzo_time_step*1e10/1e6 << std::endl;


	this->time_step_reg  = std::min(1.0, this->time_step_reg);
	this->time_block_reg = std::min(block_max, this->time_block_reg);
	this->time_level_reg = std::min(0, this->time_level_reg);

	if (this->num_neighbors != 0) {
		dtIrr = getNewTimeStep(this->acc_total, this->acc_irregular);
		getBlockTimeStep(dtIrr, this->time_level_irr, this->time_block_irr, this->time_step_irr);
	}
	else {
		this->time_block_irr = this->time_block_reg;
		this->time_level_irr = this->time_level_reg;
		this->time_step_irr  = this->time_step_reg;
	}

	this->current_time_irr  = 0;
	this->current_time_reg  = 0;
	this->current_block_irr = 0;
	this->current_block_reg = 0;

}
