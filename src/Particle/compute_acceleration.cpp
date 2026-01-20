#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include "../global.h"
#include "../def.h"
#include "../particle_data.h"
#include "../profiler.h"
#include "../simd_force.h"
#include <unordered_set>

// ============================================================================
// Acceleration accumulator for SoA-compatible force calculations
// Stores temporary acceleration values in contiguous arrays for vectorization
// ============================================================================
struct AccumulatorSoA {
    double a[3];      // acceleration components (x, y, z)
    double adot[3];   // jerk components (x, y, z)

    AccumulatorSoA() {
        a[0] = a[1] = a[2] = 0.0;
        adot[0] = adot[1] = adot[2] = 0.0;
    }

    void reset() {
        a[0] = a[1] = a[2] = 0.0;
        adot[0] = adot[1] = adot[2] = 0.0;
    }

    // Add contribution from a particle at distance dx with velocity difference dv
    void add_contribution(double mass, double r2, const double dx[3], const double dv[3], double dxdv) {
        double m_r3 = mass / (r2 * std::sqrt(r2));
        double coeff_adot = -3.0 * dxdv / r2;
        for (int d = 0; d < 3; d++) {
            a[d] += m_r3 * dx[d];
            adot[d] += m_r3 * (dv[d] + coeff_adot * dx[d]);
        }
    }

    // Subtract contribution (for neighbor transitions)
    void sub_contribution(double mass, double r2, const double dx[3], const double dv[3], double dxdv) {
        double m_r3 = mass / (r2 * std::sqrt(r2));
        double coeff_adot = -3.0 * dxdv / r2;
        for (int d = 0; d < 3; d++) {
            a[d] -= m_r3 * dx[d];
            adot[d] -= m_r3 * (dv[d] + coeff_adot * dx[d]);
        }
    }
};

// ============================================================================
// Helper: Predict neighbor position/velocity using SoA data
// ============================================================================
static inline void predict_neighbor_soa(const ParticleData& data, size_t j, double dt,
                                         double pos_out[3], double vel_out[3]) {
    // Use the free function version
    predict_second_order(data, j, dt, pos_out, vel_out);
}


void Particle::compute_acceleration_irr() {

	this->new_num_members = 0; // for Few-body Search by EW 2025.3.1

	if (this->num_neighbors == 0) {
		for (int dim=0; dim<DIM; dim++){
			this->new_position[dim] = this->position[dim];
			this->new_velocity[dim] = this->velocity[dim];
		}
		// Phase 15: Record zero neighbors with zero time
		PROFILE_NEIGHBOR_TIME(0, 0);
		return;
	}

	// Phase 15: Record initial neighbor count and start timing
	int initial_neighbor_count = this->num_neighbors;
	auto particle_start_time = std::chrono::high_resolution_clock::now();

	// Phase 19: Initialize cache counters on first call
	static bool cache_init_done = false;
	if (!cache_init_done) {
		PROFILE_CACHE_INIT();
		cache_init_done = true;
	}

	double dt, mdot, epsilon=1e-6;
	double new_time; // 0 for current and 1 for advanced times

	double x[DIM], v[DIM]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a_tmp[DIM], adot_tmp[DIM]; // 0 for current and 1 for predicted accelerations
	double pos[DIM], vel[DIM];
	double pos_neighbor[DIM], vel_neighbor[DIM];
	double m_r3;
	Particle* ptcl;
	new_time = this->current_time_irr + this->time_step_irr; // the time to be advanced to
	dt       = this->time_step_irr*enzo_time_step; // interval of time step


	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<DIM; dim++){
		a_tmp[dim]    = 0.0;
		adot_tmp[dim] = 0.0;
	}

	std::unordered_set<int> CMPtclsSet;

	/*******************************************************
	 * Irregular Acceleartion Calculation
	 ********************************************************/
	PROFILE_START(TimerID::IrregularPredict);
	this->predict_particle_second_order(this->time_step_irr, pos, vel);
	PROFILE_STOP(TimerID::IrregularPredict);

	int neighbor_pairs = 0;  // Track neighbor pairs for throughput

	// Phase 19: Start cache counters for neighbor loop
	auto cache_start = std::chrono::high_resolution_clock::now();
	PROFILE_CACHE_START();

	PROFILE_START(TimerID::IrregularNeighborLoop);

	// Pre-gather neighbor data into aligned buffers for SIMD processing
	NeighborBatch batch;
	int cm_indices[NEIGHBOR_BATCH_SIZE];
	int cm_indices_count = 0;

	gather_neighbor_data(particles, neighbors, this->neighbors_offset,
	                     this->num_neighbors, new_time, batch,
	                     cm_indices, cm_indices_count);

	// Compute forces using vectorized kernel
	double a_nb[3] = {0.0, 0.0, 0.0};
	double adot_nb[3] = {0.0, 0.0, 0.0};

	compute_force_vectorized(pos, vel, batch, a_nb, adot_nb);

	// Accumulate into temporary acceleration
	for (int dim = 0; dim < DIM; dim++) {
		a_tmp[dim] += a_nb[dim];
		adot_tmp[dim] += adot_nb[dim];
	}

	neighbor_pairs = batch.count;

	// Check for new members (used for few-body search)
	// This uses the pre-gathered data to avoid redundant prediction
	for (int i = 0; i < batch.count; i++) {
		double dx = batch.pos_x[i] - pos[0];
		double dy = batch.pos_y[i] - pos[1];
		double dz = batch.pos_z[i] - pos[2];
		double dvx = batch.vel_x[i] - vel[0];
		double dvy = batch.vel_y[i] - vel[1];
		double dvz = batch.vel_z[i] - vel[2];
		double r2_check = dx*dx + dy*dy + dz*dz;
		double vx_check = dx*dvx + dy*dvy + dz*dvz;
		if (sqrt(r2_check) < r_search && vx_check < 0) {
			this->new_members[this->new_num_members++] = batch.indices[i];
		}
	}

	// Build CM particle set from gathered indices
	for (int i = 0; i < cm_indices_count; i++) {
		CMPtclsSet.insert(cm_indices[i]);
	}

	PROFILE_STOP(TimerID::IrregularNeighborLoop);
	PROFILE_WORK(TimerID::IrregularPairsEvaluated, neighbor_pairs);

	PROFILE_START(TimerID::IrregularCMLoop);
	int cm_pairs = 0;

	if (!CMPtclsSet.empty()) {
		// Convert set to array for vectorized processing
		int cm_array[NEIGHBOR_BATCH_SIZE];
		int cm_array_count = 0;
		for (int idx : CMPtclsSet) {
			if (cm_array_count < NEIGHBOR_BATCH_SIZE) {
				cm_array[cm_array_count++] = idx;
			}
		}

		// Gather CM particle data
		NeighborBatch cm_batch;
		gather_cm_particle_data(particles, cm_array, cm_array_count,
		                        this->pid, new_time, cm_batch);

		// Compute forces using vectorized kernel
		double a_cm[3] = {0.0, 0.0, 0.0};
		double adot_cm[3] = {0.0, 0.0, 0.0};

		compute_force_vectorized(pos, vel, cm_batch, a_cm, adot_cm);

		// Accumulate into temporary acceleration
		for (int dim = 0; dim < DIM; dim++) {
			a_tmp[dim] += a_cm[dim];
			adot_tmp[dim] += adot_cm[dim];
		}

		cm_pairs = cm_batch.count;

		// Check for new members from CM particles
		for (int i = 0; i < cm_batch.count; i++) {
			double dx = cm_batch.pos_x[i] - pos[0];
			double dy = cm_batch.pos_y[i] - pos[1];
			double dz = cm_batch.pos_z[i] - pos[2];
			double dvx = cm_batch.vel_x[i] - vel[0];
			double dvy = cm_batch.vel_y[i] - vel[1];
			double dvz = cm_batch.vel_z[i] - vel[2];
			double r2_check = dx*dx + dy*dy + dz*dz;
			double vx_check = dx*dvx + dy*dvy + dz*dvz;
			if (sqrt(r2_check) < r_search && vx_check < 0) {
				this->new_members[this->new_num_members++] = cm_batch.indices[i];
			}
		}
	}

	PROFILE_STOP(TimerID::IrregularCMLoop);
	PROFILE_WORK(TimerID::IrregularPairsEvaluated, cm_pairs);

	// Phase 19: Stop cache counters after neighbor and CM loops
	auto cache_end = std::chrono::high_resolution_clock::now();
	double cache_elapsed_s = std::chrono::duration<double>(cache_end - cache_start).count();
	PROFILE_CACHE_STOP(cache_elapsed_s);

	double a2, a3, da_dt2, adot_dt, dt2, dt3, dt4, dt5;
	double dt_ex = (new_time - this->current_time_reg)*enzo_time_step;

	double A, B, C;

	dt2 = dt*dt;
	dt3 = dt2*dt;
	dt4 = dt3*dt;
	dt5 = dt4*dt; // VERY IMPORTANT BUG FIXED by EW 2025.7.18

	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	PROFILE_START(TimerID::IrregularCorrection);
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
		this->new_position[dim] = pos[dim] + a2*dt + a3*dt;
		this->new_velocity[dim] = vel[dim] + 4*a2  + 5*a3;


		// note that these higher order terms and lowers have different neighbors
		this->acc_irregular[dim][0] = a_tmp[dim];
		this->acc_irregular[dim][1] = adot_tmp[dim];
		this->acc_irregular[dim][2] = a2*24/dt3;
		this->acc_irregular[dim][3] = a3*120/dt4;
#else
		// do the higher order correcteion
		da_dt2  = (this->acc_irregular[dim][0] - a_tmp[dim]) / dt2;
		adot_dt = (this->acc_irregular[dim][1] + adot_tmp[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*this->acc_irregular[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		// 4th order correction
		// save the values in the temporary variables
		this->new_position[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->new_velocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		// note that these higher order terms and lowers have different neighbors
		this->acc_irregular[dim][0] = a_tmp[dim];
		this->acc_irregular[dim][1] = adot_tmp[dim];
		this->acc_irregular[dim][2] = a2;
		this->acc_irregular[dim][3] = a3;
#endif
	}
	PROFILE_STOP(TimerID::IrregularCorrection);

	// Phase 15: Record total neighbor count (including CM particles) with timing
	auto particle_end_time = std::chrono::high_resolution_clock::now();
	long long particle_compute_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
	    particle_end_time - particle_start_time).count();

	// Total neighbors = original + CM particles processed
	int total_neighbors = neighbor_pairs + cm_pairs;
	PROFILE_NEIGHBOR_TIME(total_neighbors, particle_compute_ns);

	// Phase 17: Record worker compute time and check for heavy particles
	// Phase 18: Record particle type (CM vs regular)
#ifdef PERFORMANCETRACE
	{
		auto& profiler = Profiler::instance();
		// Phase 21: Use global MPI rank directly - workers are ranks 1-N, root is 0
		// getCurrentParticleWorkerRank() returns 0 on workers because recordWorkerAssignment()
		// is only called on root. Workers need to use my_rank for their own identity.
		int worker_rank = my_rank;
		if (worker_rank > 0) {
			profiler.recordWorkerComputeTime(worker_rank, particle_compute_ns);

			// Check if this is a heavy particle (>2σ compute time)
			const auto& neighbor_stats = profiler.getIntervalNeighborStats();
			if (neighbor_stats.count() > 10 &&
			    neighbor_stats.isOutlier(particle_compute_ns, 2.0)) {
				profiler.recordHeavyParticle(this->particle_index, worker_rank,
				                             particle_compute_ns, total_neighbors);
			}
		}
		// Phase 18: Track CM vs regular particle
		profiler.recordParticleType(this->is_cm_particle, particle_compute_ns);
	}
#endif

	for (int dim=0; dim<DIM; dim++) {
		this->acc_total[dim][0] = this->acc_regular[dim][0] + this->acc_irregular[dim][0] + this->acc_regular[dim][1]*dt_ex; // affect the next
		this->acc_total[dim][1] = this->acc_regular[dim][1] + this->acc_irregular[dim][1];
		this->acc_total[dim][2] = this->acc_regular[dim][2] + this->acc_irregular[dim][2];
		this->acc_total[dim][3] = this->acc_regular[dim][3] + this->acc_irregular[dim][3];
	}
}




// Modified by EW 2025.3.3 for no CUDA version

void Particle::compute_acceleration_reg() {

	double dt, mdot, epsilon=1e-6;
	double new_time; // 0 for current and 1 for advanced times

	double x[DIM], v[DIM]; // 0 for current and 1 for predicted positions and velocities
	double r2, vx; // 0 for current and 1 for predicted values
	double a[DIM], adot[DIM], a_new[DIM], adot_new[DIM]; // 0 for current and 1 for predicted accelerations
	double pos[DIM], vel[DIM];
	double pos_neighbor[DIM], vel_neighbor[DIM];
	double m_r3;
	int    j=0;
	Particle* ptcl;
	new_time = this->current_time_reg+this->time_step_reg; // the time to be advanced to
	dt       = this->time_step_reg*enzo_time_step; // interval of time step
	if (dt == 0.0) {
		dt  = this->time_block_reg*time_step*enzo_time_step; // interval of time step
		assert(dt != 0.0);
	}
	this->new_num_neighbors = 0;

	std::unordered_set<int> RealNeighbors; // neighbors containing CM ptcls, not members
	RealNeighbors.reserve(this->num_neighbors);

	// initialize irregular force terms for ith particle just in case
	for (int dim=0; dim<DIM; dim++){
		a[dim]        = 0.0;
		adot[dim]     = 0.0;
		a_new[dim]    = 0.0;
		adot_new[dim] = 0.0;
		this->acc_irregular[dim][0] = 0.;
		this->acc_irregular[dim][1] = 0.;
	}

	for (int i=0; i<this->num_neighbors; i++) {
		ptcl = &particles[neighbors[this->neighbors_offset + i]];
		if (ptcl->is_active)
			RealNeighbors.insert(neighbors[this->neighbors_offset + i]);
		else if (ptcl->cm_particle_index != -1)
			RealNeighbors.insert(ptcl->cm_particle_index);
	}
	if (this->is_cm_particle)
		RealNeighbors.erase(this->particle_index);


	/*******************************************************
	 * Regular Acceleartion Calculation
	 ********************************************************/
	if (this->num_neighbors == 0)
		this->predict_particle_second_order(this->time_step_reg, pos, vel);
	else
		this->predict_particle_second_order(0, pos, vel);


	for (int i=0; i<=g_state->last_particle_index; i++) {
		ptcl = &particles[i];

		if (!ptcl->is_active || this->pid == ptcl->pid)
			continue;

		// reset temporary variables at the start of a new calculation
		r2 = 0.0;
		vx = 0.0;

		if (ptcl->num_neighbors == 0)
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_reg, pos_neighbor, vel_neighbor);
		else
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		for (int dim=0; dim<DIM; dim++) {
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

		m_r3 = ptcl->mass/r2/sqrt(r2);


		//std::cout << "PIDs=" <<  this->neighbors[j] << ', ' << ptcl->pid << NumberOfNeighbor<< std::endl;
		//if (this->neighbors[j] == ptcl->pid) {
		if (RealNeighbors.find(ptcl->particle_index) != RealNeighbors.end()) {
			//std::cout << this->pid << ", PIDs=" <<  this->neighbors[j] << ", " << ptcl->pid << std::endl;
			j++;
		} 
		else {
			for (int dim=0; dim<DIM; dim++){
				a[dim]    += m_r3*x[dim];
				adot[dim] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}


		if (r2 < this->neighbor_radius_sq) {
			if (!ptcl->is_cm_particle) {
				new_neighbors[this->neighbors_offset + this->new_num_neighbors] = ptcl->particle_index;
				this->new_num_neighbors++;
			}
			else {
				for (int k=0; k<ptcl->num_members; k++) {
					new_neighbors[this->neighbors_offset + this->new_num_neighbors] = ptcl->members[k];
					this->new_num_neighbors++;
				}
			}
			for (int dim=0; dim<DIM; dim++){
				this->acc_irregular[dim][0] += m_r3*x[dim];
				this->acc_irregular[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
		else {
			for (int dim=0; dim<DIM; dim++){
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
	for (int dim=0; dim<DIM; dim++) {

		// do the higher order correcteion
		da_dt2  = (this->acc_regular[dim][0] - a[dim]) / dt2;
		adot_dt = (this->acc_regular[dim][1] + adot[dim]) / dt;
		a2 =  -6*da_dt2  - 2*adot_dt - 2*this->acc_regular[dim][1]/dt;
		a3 =  (12*da_dt2 + 6*adot_dt)/dt;

		//fprintf(stderr, "DIM=%d, pid=%d, a=%.2e, da_dt2=%.2e, adot_dt=%.2e, a2=%.2e, a3=%.2e\n", dim,this->pid,a[dim], da_dt2, adot_dt, a2, a3);

		// 4th order correction
		// save the values in the temporary variables
		this->new_position[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->new_velocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;


		this->acc_regular[dim][2] = a2;
		this->acc_regular[dim][3] = a3;
	}


	for (int dim=0; dim<DIM; dim++) {
		this->acc_regular[dim][0] = a_new[dim];
		this->acc_regular[dim][1] = adot_new[dim];
		this->acc_total[dim][0] = this->acc_regular[dim][0] + this->acc_irregular[dim][0]; 
		this->acc_total[dim][1] = this->acc_regular[dim][1] + this->acc_irregular[dim][1];
		if (this->new_num_neighbors == 0) {
			this->acc_total[dim][2] = this->acc_regular[dim][2];
			this->acc_total[dim][3] = this->acc_regular[dim][3];
		}
	}
}




// Modified by EW 2025.1.30

void Particle::update_regular_particle_cuda() {

	double new_a[DIM], new_adot[DIM];
	double new_time = this->current_time_reg+this->time_step_reg;
	double pos[DIM], vel[DIM];
	if (this->num_neighbors == 0)
		this->predict_particle_second_order(this->time_step_reg, pos, vel);
	else
		this->predict_particle_second_order(0, pos, vel);

	double a_tmp[DIM], adot_tmp[DIM];

	for (int dim=0; dim<DIM; dim++) {
		new_a[dim]			= this->acc_irregular[dim][0];
		new_adot[dim]		= this->acc_irregular[dim][1];
		this->acc_irregular[dim][0] = 0.;
		this->acc_irregular[dim][1] = 0.;

		a_tmp[dim]          = 0.;
		adot_tmp[dim]       = 0.;
	}

	std::unordered_set<int> hashTableOld;
	hashTableOld.reserve(this->num_neighbors);
	std::unordered_set<int> hashTableNew;
	hashTableNew.reserve(this->new_num_neighbors); // We are including myself in neighbor from GPU kernel by EW 2025.8.26


	hashTableNew.insert(new_neighbors + this->neighbors_offset, new_neighbors + this->neighbors_offset + this->new_num_neighbors);
	hashTableNew.erase(this->particle_index);
	this->new_num_neighbors--;

	for (int i = 0; i < this->num_neighbors; i++) {
		if (particles[neighbors[this->neighbors_offset + i]].is_active)
			hashTableOld.insert(neighbors[this->neighbors_offset + i]);
		else if (particles[neighbors[this->neighbors_offset + i]].cm_particle_index != -1)
			hashTableOld.insert(particles[neighbors[this->neighbors_offset + i]].cm_particle_index);
	}

	Particle* ptcl;
	double pos_neighbor[DIM], vel_neighbor[DIM];
	double dx[DIM], dv[DIM];
	double dr2;
	double dxdv;
	double m_r3;

	for (int _OldNeighborIndex: hashTableOld) {
		ptcl = &particles[_OldNeighborIndex];

		// neighbor in old but not in new
		if ( hashTableNew.find(ptcl->particle_index) == hashTableNew.end() ) {

			if (ptcl->num_neighbors == 0)
				ptcl->predict_particle_second_order(new_time-ptcl->current_time_reg, pos_neighbor, vel_neighbor);
			else
				ptcl->predict_particle_second_order(new_time-ptcl->current_time_irr, pos_neighbor, vel_neighbor);

			dr2  = 0.0;
			dxdv = 0.0;
			for (int dim=0; dim<DIM; dim++) {
				dx[dim] = pos_neighbor[dim] - pos[dim];
				dv[dim] = vel_neighbor[dim] - vel[dim];
				dr2    += dx[dim]*dx[dim];
				dxdv   += dx[dim]*dv[dim];
			}

			m_r3 = ptcl->mass/dr2/sqrt(dr2);

			for (int dim=0; dim<DIM; dim++){
				a_tmp[dim]    -= m_r3*dx[dim];
				adot_tmp[dim] -= m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
			}
		}
	}

	for (int _NewNeighborIndex: hashTableNew) {
		ptcl = &particles[_NewNeighborIndex];

		if (ptcl->num_neighbors == 0)
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_reg, pos_neighbor, vel_neighbor);
		else
			ptcl->predict_particle_second_order(new_time-ptcl->current_time_irr, pos_neighbor, vel_neighbor);

		dr2  = 0.0;
		dxdv = 0.0;
		for (int dim=0; dim<DIM; dim++) {
			dx[dim] = pos_neighbor[dim] - pos[dim];
			dv[dim] = vel_neighbor[dim] - vel[dim];
			dr2    += dx[dim]*dx[dim];
			dxdv   += dx[dim]*dv[dim];
		}

		m_r3 = ptcl->mass/dr2/sqrt(dr2);

		for (int dim=0; dim<DIM; dim++){
			this->acc_irregular[dim][0] += m_r3*dx[dim];
			this->acc_irregular[dim][1] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
		}

		// neighbor in new but not in old
		if ( hashTableOld.find(ptcl->particle_index) == hashTableOld.end() ) {
			for (int dim=0; dim<DIM; dim++){
				a_tmp[dim]    += m_r3*dx[dim];
				adot_tmp[dim] += m_r3*(dv[dim] - 3*dx[dim]*dxdv/dr2);
			}
		}
	}


	/*******************************************************
	 * Position and velocity correction due to 4th order correction
	 ********************************************************/
	double dt  = this->time_step_reg*enzo_time_step;  // unit conversion
	if (dt == 0.0) {
		dt  = this->time_block_reg*time_step*enzo_time_step; // unit conversion
		assert(dt != 0.0);
	}
	double dt2 = dt*dt;
	double dt3 = dt2*dt;
	double dt4 = dt3*dt;
	double dt5 = dt4*dt;
	double da_dt2, adot_dt, a2, a3;

	//fprintf(stdout, "PID=%d\n", ptcl->pid);
	for (int dim=0; dim<DIM; dim++) {
		da_dt2  = (this->acc_regular[dim][0] - new_a[dim]   - a_tmp[dim]    ) / dt2;
		adot_dt = (this->acc_regular[dim][1] + new_adot[dim] + adot_tmp[dim]) / dt;


		a2 =  -6*da_dt2 - 2*adot_dt - 2*this->acc_regular[dim][1]/dt;
		a3 = (12*da_dt2 + 6*adot_dt)/dt;

		// note that these higher order terms and lowers have different neighbors

		// 4th order correction
		// save the values in the temporary variables
		this->new_position[dim] = pos[dim] + a2*dt4/24 + a3*dt5/120;
		this->new_velocity[dim] = vel[dim] + a2*dt3/6  + a3*dt4/24;

		this->acc_regular[dim][2] = a2;
		this->acc_regular[dim][3] = a3;

		// reset for future use
		a_tmp[dim]    = 0.;
		adot_tmp[dim] = 0.;
	}

	int _NewNumberOfNeighbor = 0;
	for (int _NewNeighborIndex: hashTableNew) {
		ptcl = &particles[_NewNeighborIndex];
		if (ptcl->is_cm_particle) {
			for (int j=0; j<ptcl->num_members; j++) {
				new_neighbors[this->neighbors_offset + _NewNumberOfNeighbor] = ptcl->members[j];
				_NewNumberOfNeighbor++;
			}
		}
		else {
			new_neighbors[this->neighbors_offset + _NewNumberOfNeighbor] = ptcl->particle_index;
			_NewNumberOfNeighbor++;
		}
	}
	this->new_num_neighbors = _NewNumberOfNeighbor;


	for (int dim=0; dim<DIM; dim++) {
		this->acc_regular[dim][0] = new_a[dim];
		this->acc_regular[dim][1] = new_adot[dim];
		this->acc_total[dim][0] = this->acc_regular[dim][0] + this->acc_irregular[dim][0];
		this->acc_total[dim][1] = this->acc_regular[dim][1] + this->acc_irregular[dim][1];
		if (this->new_num_neighbors == 0) {
			this->acc_total[dim][2] = this->acc_regular[dim][2];
			this->acc_total[dim][3] = this->acc_regular[dim][3];
		}
	}
}

