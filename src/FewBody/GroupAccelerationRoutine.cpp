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
		if (ptcl->NewNumberOfMember > 0) {
			NewCM = &particles[LastParticleIndex+1];
			NewCM->clear();
			NewCM->copyNewMembers(ptcl);
			NewCM->NewMembers[NewCM->NewNumberOfMember++] = ptcl->ParticleIndex;

			LastParticleIndex++;
		}
	}

	if (OriginalLastParticleIndex == LastParticleIndex) return;

	mergeGroupCandidates(OriginalLastParticleIndex);	// Merge group candidates
					// ex) A & B are a group and B & C are a group --> Merge so that A & B & C become one group!
	
	for (int i=OriginalLastParticleIndex+1; i<=LastParticleIndex; i++) {
		// deleteNeighbors(i);
		NewCM = &particles[i];
		NewCM->ParticleIndex = i;
		NewCM->PID = NewCMPID;
		NewCMPID++;
		std::cout << "New Primordial CM ParticleIndex: " << i << std::endl;
		std::cout << "New Primordial CM PID: " << NewCM->PID << std::endl;
		NewCM->setBinaryInterruptState(BinaryInterruptState::none);

		NumberOfParticle += 1 - NewCM->NewNumberOfMember;
	}
	global_variable->LastParticleIndex = LastParticleIndex;
}

void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls,
 std::unordered_map<int,int>& existing, std::unordered_map<int,int>& terminated) {

	Particle* ptcl;
	Particle* NewCM;

	int OriginalLastParticleIndex = LastParticleIndex;

	for (int i=0; i<ParticleList.size(); i++) {
		ptcl = &particles[ParticleList[i]];
		if (ptcl->NewNumberOfMember > 0) {
			// /* // test by EW 2025.9.18 // It seems to work well by EW 2025.9.20
			int OriginalNewNumberOfMember = ptcl->NewNumberOfMember;
			for (int j = 0; j < OriginalNewNumberOfMember; j++) {
				Particle* member = &particles[ptcl->NewMembers[j]];
				if (!member->isActive) {
					fprintf(stdout, "Inactive member (PID: %d) found in formBinaries! Remove it from the group candidate of PID %d\n", member->PID, ptcl->PID);
					if (j < OriginalNewNumberOfMember - 1) {
						ptcl->NewMembers[j] = ptcl->NewMembers[OriginalNewNumberOfMember - 1];
						j--;
					}
					else
						ptcl->NewNumberOfMember--;
				}
			}
			if (ptcl->NewNumberOfMember <= 0) {
				fprintf(stdout, "Skipping forming group for PID %d since there are no more NewMembers!\n", ptcl->PID);
				continue;
			}
			// */
			// fprintf(stdout, "GAR. Num: %d\n", ptcl->NewNumberOfMember + 1); // for debugging by EW 2025.1.23
			// fprintf(stdout, "GAR. PID: %d\n", ptcl->PID); // for debugging by EW 2025.1.23
			NewCM = &particles[LastParticleIndex+1];
			NewCM->clear();
			NewCM->copyNewMembers(ptcl);
			/* // for debugging by EW 2025.1.23
			for (int j = 0; j < ptcl->NewNumberOfMember; j++) {
				fprintf(stdout, "GAR. PID: %d\n", particles[ptcl->NewMembers[j]].PID);
			}
			*/
			NewCM->NewMembers[NewCM->NewNumberOfMember++] = ptcl->ParticleIndex;

			LastParticleIndex++;
		}
	}

	if (OriginalLastParticleIndex == LastParticleIndex) return;

	// fprintf(stdout, "A. global_variable->NOP: %d, NOP: %d\n", global_variable->LastParticleIndex, LastParticleIndex);

	mergeGroupCandidates(OriginalLastParticleIndex);	// Merge group candidates
					// ex) A & B form a group and B & C form a group --> Merge so that A & B & C become one group!
	
	// fprintf(stdout, "B. global_variable->NOP: %d, NOP: %d\n", global_variable->LastParticleIndex, LastParticleIndex);

	assert(LastParticleIndex > OriginalLastParticleIndex); 

	while (terminated.size() != 0) {

		Particle* ptcl = &particles[LastParticleIndex];
		auto it = terminated.begin();
		particles[it->first].copyNewMembers(ptcl);

		existing.insert({it->first, it->second});
		newCMptcls.push_back(it->first);
		terminated.erase(it);
		LastParticleIndex--;
		if (OriginalLastParticleIndex == LastParticleIndex)
			break;
	}
	if (OriginalLastParticleIndex != LastParticleIndex) {
		for (int i = OriginalLastParticleIndex+1; i <= LastParticleIndex; i++) {
			existing.insert({i, existing.size() % NumberOfWorker + 1});
			newCMptcls.push_back(i);
		}
	}

	for (int i: newCMptcls) {
		// deleteNeighbors(i);
		NewCM = &particles[i];
		if (NewCM->NewNumberOfMember < 2) {
			fprintf(stderr, "Error in GroupAcceleratonRoutine.cpp: NewCM->NewNumberOfMember < 2\n");
			for (int j=0; j<NewCM->NewNumberOfMember; j++) {
				ptcl = &particles[NewCM->NewMembers[j]];
				fprintf(stderr, "NewCM->NewMembers[%d]: %d\n", j, NewCM->NewMembers[j]);
				ptcl->NewNumberOfMember = 0;
				assert(ptcl->isActive);
			}
			if (i != LastParticleIndex)
				terminated.insert({i, existing[i]});
			else
				LastParticleIndex--;
			existing.erase(i);
			NewCM->clear();
			continue;
		}
		NewCM->ParticleIndex = i;
		NewCM->NeighborsOffset = NewCM->ParticleIndex * MaxNumNeighbor;
		NewCM->PID = NewCMPID;
		NewCMPID++;
#ifdef DEBUG
		std::cout << "New CM ParticleIndex: " << i << std::endl;
		std::cout << "New CM PID: " << NewCM->PID << std::endl;
#endif
		NewCM->setBinaryInterruptState(BinaryInterruptState::none);

		NumberOfParticle += 1 - NewCM->NewNumberOfMember;
	}
	global_variable->LastParticleIndex = LastParticleIndex;

	ParticleList.insert(ParticleList.end(), newCMptcls.begin(), newCMptcls.end());
}


void mergeGroupCandidates(int OriginalLastParticleIndex) {

    bool merged = true;

    while (merged) {
        merged = false;
		for (int i = OriginalLastParticleIndex+1; i <= LastParticleIndex; i++) {

			Particle* currentCM = &particles[i];
			if (currentCM->NewNumberOfMember == 0) continue; // Skip already deleted groups

			for (int j = i + 1; j <= LastParticleIndex; j++) {

				Particle* otherCM = &particles[j];
				if (otherCM->NewNumberOfMember == 0) continue; // Skip already deleted groups

                // Check if there's any common member between group1 and group2
                bool commonFound = false;
				for (int k=0; k < currentCM->NewNumberOfMember; k++) {
					int member1 = currentCM->NewMembers[k];
					if (std::find(otherCM->NewMembers, otherCM->NewMembers + otherCM->NewNumberOfMember, member1) != otherCM->NewMembers + otherCM->NewNumberOfMember) {
						commonFound = true;
						break;
					}
                }

                // If common members are found, merge group2 into group1
                if (commonFound) {
                    // Merge group2 into group1, avoiding duplicates
					for (int l=0; l < otherCM->NewNumberOfMember; l++) {
						int member2 = otherCM->NewMembers[l];
						if (std::find(currentCM->NewMembers, currentCM->NewMembers + currentCM->NewNumberOfMember, member2) == currentCM->NewMembers + currentCM->NewNumberOfMember) {
							currentCM->NewMembers[currentCM->NewNumberOfMember] = member2;
							currentCM->NewNumberOfMember++;
						}
					}
					merged = true;

                    // Mark otherGroup for deletion after the loop
					if (j != LastParticleIndex)
						particles[j].copyNewMembers(&particles[LastParticleIndex]);
					particles[LastParticleIndex].clear();
					LastParticleIndex--;
                }
            }
        }
    }
}


void makePrimordialGroup(Particle* ptclCM) {

	ptclCM->isActive = true;
	ptclCM->isCMptcl = true;

	Group* ptclGroup = new Group();

	ptclCM->GroupInfo = ptclGroup;
	ptclGroup->groupCM = ptclCM;

	for (int i = 0; i < ptclCM->NewNumberOfMember; ++i) {
		Particle* members = &particles[ptclCM->NewMembers[i]];
		members->isActive = false;
    }

	ptclGroup->initialManager();
	ptclGroup->initialIntegrator(ptclCM->NewNumberOfMember); // Binary tree is made and CM particle is made automatically.

	// ptclCM = &ptclGroup->sym_int.particles.cm;
	for (int dim=0; dim<Dim; dim++) {
		ptclCM->Position[dim] = ptclGroup->sym_int.particles.cm.Position[dim];
		ptclCM->Velocity[dim] = ptclGroup->sym_int.particles.cm.Velocity[dim];
		ptclCM->Mass = ptclGroup->sym_int.particles.cm.Mass;
	}

	// ptclCM->RadiusOfNeighbor = InitialNeighborRadius*InitialNeighborRadius;
	ptclCM->RadiusOfNeighbor = particles[ptclCM->NewMembers[0]].RadiusOfNeighbor;

	fprintf(workerout, "The ID of CM is %d.\n",ptclCM->PID);

	fprintf(workerout, "------------------NEW-GROUP-MEMBER-INFORMATION------------------\n");
	for (int i=0; i < ptclGroup->sym_int.particles.getSize(); i++) {
		Particle* members = &ptclGroup->sym_int.particles[i];
		members->printParticleInfo(workerout);
    }

	ptclGroup->sym_int.initialIntegration(0); // This is primordial binary!
    ptclGroup->sym_int.info.calcDsAndStepOption(ptclGroup->manager.step.getOrder(), ptclGroup->manager.interaction.gravitational_constant, ptclGroup->manager.ds_scale);


// /* // Eunwoo test
	auto& bin_root = ptclGroup->sym_int.info.getBinaryTreeRoot();
	if (bin_root.semi>0.0) {
		ptclGroup->sym_int.info.r_break_crit = fmin(2*bin_root.semi, sqrt(ptclCM->RadiusOfNeighbor));
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
		fprintf(workerout, "period: %e Myr\n", bin_root.period*1e4);
		fprintf(workerout, "t_peri: %e Myr\n\t", abs(bin_root.t_peri*1e4));
		fprintf(workerout, "r_break_crit: %e pc\n", ptclGroup->sym_int.info.r_break_crit*position_unit);
	}
// */ // Eunwoo test

	// Find neighbors for CM particle and calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately 
	CalculateAcceleration01(ptclCM);
	CalculateAcceleration23(ptclCM);

	fprintf(workerout, "\nResult of CM particle value calculation from function NewPrimordialBinaries\n");

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

	fprintf(workerout, "------------------END-OF-NEW-PRIMORDIAL-BINARIES------------------\n\n");
	fflush(workerout);
}

#endif

