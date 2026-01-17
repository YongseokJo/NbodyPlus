#ifdef FEWBODY
#include <map>
#include <unordered_map>
#include "../global.h"

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void mergeGroupCandidates(int OriginalLastParticleIndex);

void formPrimordialBinaries(int OriginalLastParticleIndex) {

	Particle* ptcl;
	Particle* NewCM;

	for (int i=0; i<=OriginalLastParticleIndex; i++) {
		ptcl = &particles[i];
		if (ptcl->new_num_members > 0) {
			NewCM = &particles[last_particle_index+1];
			NewCM->clear();
			NewCM->copy_new_members(ptcl);
			NewCM->new_members[NewCM->new_num_members++] = ptcl->particle_index;

			last_particle_index++;
		}
	}

	if (OriginalLastParticleIndex == last_particle_index) return;

	mergeGroupCandidates(OriginalLastParticleIndex);	// Merge group candidates
					// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!
	
	for (int i=OriginalLastParticleIndex+1; i<=last_particle_index; i++) {
		// deleteNeighbors(i);
		NewCM = &particles[i];
		NewCM->particle_index = i;
		NewCM->pid = new_cm_pid;
		new_cm_pid++;
		std::cout << "New Primordial CM ParticleIndex: " << i << std::endl;
		std::cout << "New Primordial CM PID: " << NewCM->pid << std::endl;
		NewCM->set_binary_interrupt_state(BinaryInterruptState::none);

		num_particles += 1 - NewCM->new_num_members;
	}
	g_state->last_particle_index = last_particle_index;
}

void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls,
 std::unordered_map<int,int>& existing, std::unordered_map<int,int>& terminated) {

	Particle* ptcl;
	Particle* NewCM;

	int OriginalLastParticleIndex = last_particle_index;

	for (int i=0; i<ParticleList.size(); i++) {
		ptcl = &particles[ParticleList[i]];
		if (ptcl->new_num_members > 0) {
			// /* // test by EW 2025.9.18 // It seems to work well by EW 2025.9.20
			int OriginalNewNumberOfMember = ptcl->new_num_members;
			for (int j = 0; j < OriginalNewNumberOfMember; j++) {
				Particle* member = &particles[ptcl->new_members[j]];
				if (!member->is_active) {
					fprintf(stdout, "Inactive member (PID: %d) found in formBinaries! Remove it from the group candidate of PID %d\n", member->pid, ptcl->pid);
					if (j < OriginalNewNumberOfMember - 1) {
						ptcl->new_members[j] = ptcl->new_members[OriginalNewNumberOfMember - 1];
						j--;
					}
					else
						ptcl->new_num_members--;
				}
			}
			if (ptcl->new_num_members <= 0) {
				fprintf(stdout, "Skipping forming group for PID %d since there are no more NewMembers!\n", ptcl->pid);
				continue;
			}
			// */
			// fprintf(stdout, "GAR. Num: %d\n", ptcl->new_num_members + 1); // for debugging by EW 2025.1.23
			// fprintf(stdout, "GAR. PID: %d\n", ptcl->pid); // for debugging by EW 2025.1.23
			NewCM = &particles[last_particle_index+1];
			NewCM->clear();
			NewCM->copy_new_members(ptcl);
			/* // for debugging by EW 2025.1.23
			for (int j = 0; j < ptcl->new_num_members; j++) {
				fprintf(stdout, "GAR. PID: %d\n", particles[ptcl->new_members[j]].pid);
			}
			*/
			NewCM->new_members[NewCM->new_num_members++] = ptcl->particle_index;

			last_particle_index++;
		}
	}

	if (OriginalLastParticleIndex == last_particle_index) return;

	// fprintf(stdout, "A. g_state->NOP: %d, NOP: %d\n", g_state->last_particle_index, last_particle_index);

	mergeGroupCandidates(OriginalLastParticleIndex);	// Merge group candidates
					// ex) A & B form a group and B & C form a group --> Merge so that A & B & C become one group!
	
	// fprintf(stdout, "B. g_state->NOP: %d, NOP: %d\n", g_state->last_particle_index, last_particle_index);

	assert(last_particle_index > OriginalLastParticleIndex); 

	while (terminated.size() != 0) {

		Particle* ptcl = &particles[last_particle_index];
		auto it = terminated.begin();
		particles[it->first].copy_new_members(ptcl);

		existing.insert({it->first, it->second});
		newCMptcls.push_back(it->first);
		terminated.erase(it);
		last_particle_index--;
		if (OriginalLastParticleIndex == last_particle_index)
			break;
	}
	if (OriginalLastParticleIndex != last_particle_index) {
		for (int i = OriginalLastParticleIndex+1; i <= last_particle_index; i++) {
			existing.insert({i, existing.size() % num_workers + 1});
			newCMptcls.push_back(i);
		}
	}

	for (int i: newCMptcls) {
		// deleteNeighbors(i);
		NewCM = &particles[i];
		if (NewCM->new_num_members < 2) {
			fprintf(stderr, "TASK_ERROR in GroupAcceleratonRoutine.cpp: NewCM->new_num_members < 2\n");
			for (int j=0; j<NewCM->new_num_members; j++) {
				ptcl = &particles[NewCM->new_members[j]];
				fprintf(stderr, "NewCM->new_members[%d]: %d\n", j, NewCM->new_members[j]);
				ptcl->new_num_members = 0;
				assert(ptcl->is_active);
			}
			if (i != last_particle_index)
				terminated.insert({i, existing[i]});
			else
				last_particle_index--;
			existing.erase(i);
			NewCM->clear();
			continue;
		}
		NewCM->particle_index = i;
		NewCM->neighbors_offset = NewCM->particle_index * MAX_NUM_NEIGHBOR;
		NewCM->pid = new_cm_pid;
		new_cm_pid++;
#ifdef DEBUG
		std::cout << "New CM ParticleIndex: " << i << std::endl;
		std::cout << "New CM PID: " << NewCM->pid << std::endl;
#endif
		NewCM->set_binary_interrupt_state(BinaryInterruptState::none);

		num_particles += 1 - NewCM->new_num_members;
	}
	g_state->last_particle_index = last_particle_index;

	ParticleList.insert(ParticleList.end(), newCMptcls.begin(), newCMptcls.end());
}


void mergeGroupCandidates(int OriginalLastParticleIndex) {

    bool merged = true;

    while (merged) {
        merged = false;
		for (int i = OriginalLastParticleIndex+1; i <= last_particle_index; i++) {

			Particle* currentCM = &particles[i];
			if (currentCM->new_num_members == 0) continue; // Skip already deleted groups

			for (int j = i + 1; j <= last_particle_index; j++) {

				Particle* otherCM = &particles[j];
				if (otherCM->new_num_members == 0) continue; // Skip already deleted groups

                // Check if there's any common member between group1 and group2
                bool commonFound = false;
				for (int k=0; k < currentCM->new_num_members; k++) {
					int member1 = currentCM->new_members[k];
					if (std::find(otherCM->new_members, otherCM->new_members + otherCM->new_num_members, member1) != otherCM->new_members + otherCM->new_num_members) {
						commonFound = true;
						break;
					}
                }

                // If common members are found, merge group2 into group1
                if (commonFound) {
                    // Merge group2 into group1, avoiding duplicates
					for (int l=0; l < otherCM->new_num_members; l++) {
						int member2 = otherCM->new_members[l];
						if (std::find(currentCM->new_members, currentCM->new_members + currentCM->new_num_members, member2) == currentCM->new_members + currentCM->new_num_members) {
							currentCM->new_members[currentCM->new_num_members] = member2;
							currentCM->new_num_members++;
						}
					}
					merged = true;

                    // Mark otherGroup for deletion after the loop
					if (j != last_particle_index)
						particles[j].copy_new_members(&particles[last_particle_index]);
					particles[last_particle_index].clear();
					last_particle_index--;
                }
            }
        }
    }
}


void makePrimordialGroup(Particle* ptclCM) {

	ptclCM->is_active = true;
	ptclCM->is_cm_particle = true;

	Group* ptclGroup = new Group();

	ptclCM->group_info = ptclGroup;
	ptclGroup->groupCM = ptclCM;

	for (int i = 0; i < ptclCM->new_num_members; ++i) {
		Particle* members = &particles[ptclCM->new_members[i]];
		members->is_active = false;
    }

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->new_num_members); // Binary tree is made and CM particle is made automatically.

	// ptclCM = &ptclGroup->sym_int.particles.cm;
	for (int dim=0; dim<DIM; dim++) {
		ptclCM->position[dim] = ptclGroup->sym_int.particles.cm.position[dim];
		ptclCM->velocity[dim] = ptclGroup->sym_int.particles.cm.velocity[dim];
		ptclCM->mass = ptclGroup->sym_int.particles.cm.mass;
	}

	// ptclCM->neighbor_radius_sq = initial_neighbor_radius*initial_neighbor_radius;
	ptclCM->neighbor_radius_sq = particles[ptclCM->new_members[0]].neighbor_radius_sq;

	fprintf(worker_output_file, "The ID of CM is %d.\n",ptclCM->pid);

	fprintf(worker_output_file, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		members->print_particle_info(worker_output_file);
    }

	ptclGroup->sym_int.initialIntegration(0); // This is primordial binary!
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);


// /* // Eunwoo test
	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->neighbor_radius_sq));
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
		fprintf(worker_output_file, "period: %e Myr\n", bin_root.period*1e4);
		fprintf(worker_output_file, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(worker_output_file, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
// */ // Eunwoo test

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	CalculateAcceleration01(ptclCM);
	CalculateAcceleration23(ptclCM);

	fprintf(worker_output_file, "\nResult of CM particle value calculation from function NewPrimordialBinaries\n");

	fprintf(worker_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", ptclCM->position[0]*position_unit, ptclCM->position[1]*position_unit, ptclCM->position[2]*position_unit);
	fprintf(worker_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", ptclCM->velocity[0]*velocity_unit/yr*pc/1e5, ptclCM->velocity[1]*velocity_unit/yr*pc/1e5, ptclCM->velocity[2]*velocity_unit/yr*pc/1e5);
	fprintf(worker_output_file, "Mass (Msol) - %e, \n", ptclCM->mass*mass_unit);
	fprintf(worker_output_file, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->acc_total[0][0], ptclCM->acc_total[1][0], ptclCM->acc_total[2][0]);
	// fprintf(worker_output_file, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->acc_total[0][1], ptclCM->acc_total[1][1], ptclCM->acc_total[2][1]);
	// fprintf(worker_output_file, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->acc_total[0][2], ptclCM->acc_total[1][2], ptclCM->acc_total[2][2]);
	// fprintf(worker_output_file, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->acc_total[0][3], ptclCM->acc_total[1][3], ptclCM->acc_total[2][3]);
	fprintf(worker_output_file, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->acc_regular[0][0], ptclCM->acc_regular[1][0], ptclCM->acc_regular[2][0]);
	// fprintf(worker_output_file, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->acc_regular[0][1], ptclCM->acc_regular[1][1], ptclCM->acc_regular[2][1]);
	// fprintf(worker_output_file, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->acc_regular[0][2], ptclCM->acc_regular[1][2], ptclCM->acc_regular[2][2]);
	// fprintf(worker_output_file, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->acc_regular[0][3], ptclCM->acc_regular[1][3], ptclCM->acc_regular[2][3]);
	fprintf(worker_output_file, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->acc_irregular[0][0], ptclCM->acc_irregular[1][0], ptclCM->acc_irregular[2][0]);
	// fprintf(worker_output_file, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->acc_irregular[0][1], ptclCM->acc_irregular[1][1], ptclCM->acc_irregular[2][1]);
	// fprintf(worker_output_file, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->acc_irregular[0][2], ptclCM->acc_irregular[1][2], ptclCM->acc_irregular[2][2]);
	// fprintf(worker_output_file, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->acc_irregular[0][3], ptclCM->acc_irregular[1][3], ptclCM->acc_irregular[2][3]);

	fprintf(worker_output_file, "------------------END-OF-NEW-PRIMORDIAL-BINARIES------------------\n\n");
	fflush(worker_output_file);
}

#endif

