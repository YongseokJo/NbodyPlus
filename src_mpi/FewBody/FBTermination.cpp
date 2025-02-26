#ifdef FEWBODY
#include "../global.h"

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);

void insertNeighbors(Particle* ptclCM) {

	// std::cout << "CM ptcl PID: " << ptclCM->PID << ", NewNumberOfNeighbor: " << ptclCM->NewNumberOfNeighbor << std::endl; // for debugging by EW 2025.1.22

	for (int i=0; i<=LastParticleIndex; i++) {
		Particle* ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;
		
		for (int j = 0; j < ptcl->NumberOfNeighbor; ++j) {
			if (ptcl->Neighbors[j] == ptclCM->ParticleIndex) {
				// std::cout << "Terminated neighbor PID: " << particles[ptcl->Neighbors[j]].PID << " of particle PID: " << ptcl->PID << std::endl; // for debugging by EW 2025.1.22
				// std::cout << "Original NumberOfNeighbor: " << ptcl->NumberOfNeighbor << std::endl; // for debugging by EW 2025.1.22
				ptcl->Neighbors[j] = ptcl->Neighbors[ptcl->NumberOfNeighbor - 1];
				ptcl->NumberOfNeighbor--;

				for (int k=0; k<ptclCM->NewNumberOfNeighbor; k++) {
					if (particles[ptclCM->NewNeighbors[k]].Mass == 0.0)
						continue;
					ptcl->Neighbors[ptcl->NumberOfNeighbor] = ptclCM->NewNeighbors[k];
					// std::cout << "Newly added neighbor PID: " << particles[ptclCM->NewNeighbors[k]].PID << std::endl; // for debugging by EW 2025.1.22
					ptcl->NumberOfNeighbor++;
				}
				// std::cout << "After NumberOfNeighbor: " << ptcl->NumberOfNeighbor << std::endl; // for debugging by EW 2025.1.22
				break;
			}
		}
	}
}

void FBTermination(Particle* ptclCM) {

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination.cpp... (CM PID: %d)\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);
	fprintf(binout, "CurrentTimeIrr of the first member (Myr): %e\n", particles[ptclCM->Members[0]].CurrentTimeIrr*EnzoTimeStep*1e4);
	fprintf(binout, "N_member: %d\n", ptclCM->NumberOfMember);

	NumberOfParticle--; // CM particle should be inactive by EW 2025.1.20
	
	Particle* members;
	for (int i=0; i<ptclCM->NumberOfMember; i++) {
		members = &particles[ptclCM->Members[i]];

		if (members->Mass == 0.0)
			continue;

		NumberOfParticle++; // by EW 2025.1.20

		if (ptclCM->NumberOfMember == 2)
			members->setBinaryInterruptState(BinaryInterruptState::none);
		else if (ptclCM->NumberOfMember > 3)
			members->setBinaryInterruptState(BinaryInterruptState::manybody);

		members->CurrentBlockIrr	= ptclCM->CurrentBlockIrr;
		members->CurrentBlockReg	= ptclCM->CurrentBlockReg;
		// members->CurrentTimeIrr		= ptclCM->CurrentTimeIrr;
		members->CurrentTimeReg		= ptclCM->CurrentTimeReg;
		members->NewCurrentBlockIrr	= ptclCM->NewCurrentBlockIrr;

		members->TimeLevelIrr		= ptclCM->TimeLevelIrr; // test by EW 2025.1.29
		members->TimeLevelReg		= ptclCM->TimeLevelReg;

		members->RadiusOfNeighbor = ACRadius*ACRadius; // added by EW 2025.1.16
		// members->RadiusOfNeighbor = ptclCM->RadiusOfNeighbor; // modified by EW 2025.1.30

		CalculateAcceleration01(members);
		CalculateAcceleration23(members);

		if (members->CurrentTimeIrr != ptclCM->CurrentTimeIrr) {
			assert(members->CurrentTimeIrr < ptclCM->CurrentTimeIrr); // for debugging by EW
			fprintf(stdout, "SDAR: binary merger happened!\n");

			double pos[Dim], vel[Dim];
			members->predictParticleSecondOrder(ptclCM->CurrentTimeIrr - members->CurrentTimeIrr, pos, vel);
			// /*
			for (int dim=0; dim<Dim; dim++) {
				members->Position[dim] =  pos[dim];
				members->Velocity[dim] =  vel[dim];
			}
			// */
			// members->correctParticleFourthOrder(ptclCM->CurrentTimeIrr - members->CurrentTimeIrr, pos, vel, members->a_tot);
			// members->updateParticle();
			members->CurrentTimeIrr = ptclCM->CurrentTimeIrr;
		}

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

		if (members->NumberOfNeighbor != 0) {
			members->calculateTimeStepIrr2();

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

			members->correctParticleFourthOrder(members->CurrentTimeIrr - members->CurrentTimeReg, members->Position, members->Velocity, members->a_tot);
			members->updateParticle();
			members->CurrentTimeReg = members->CurrentTimeIrr;
		}
/* // test11
		if (ptclCM->NumberOfMember == 2) {
			members->calculateTimeStepReg();

			if (members->TimeLevelReg <= ptclCM->TimeLevelReg-1 
					&& members->TimeBlockReg/2+members->CurrentBlockReg > ptclCM->CurrentBlockIrr)  { // this ensures that irr time of any particles is smaller than adjusted new reg time.
				members->TimeLevelReg = ptclCM->TimeLevelReg-1;
			}
			else if  (members->TimeLevelReg >= ptclCM->TimeLevelReg+1) {
				members->TimeLevelReg = ptclCM->TimeLevelReg+1;
			}
			else 
				members->TimeLevelReg = ptclCM->TimeLevelReg;
			members->TimeStepReg  = static_cast<double>(pow(2, members->TimeLevelReg));
			members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block));

			// members->calculateTimeStepIrr();
			members->calculateTimeStepIrr2();
		}
		else {
			members->calculateTimeStepReg();

			while (members->TimeStepReg*EnzoTimeStep*1e4 > 2e-7) {
				members->TimeLevelReg--;
				members->TimeStepReg  = static_cast<double>(pow(2, members->TimeLevelReg));
				members->TimeBlockReg = static_cast<ULL>(pow(2, members->TimeLevelReg-time_block));
			}
			if (members->TimeBlockReg + members->CurrentBlockReg <= ptclCM->CurrentBlockIrr) {
				members->CurrentBlockReg = NextRegTimeBlock - members->TimeBlockReg;
				members->CurrentTimeReg = members->CurrentBlockReg * time_step;
			}

			// members->calculateTimeStepIrr();
			members->calculateTimeStepIrr2();
			while (members->TimeStepIrr*EnzoTimeStep*1e4 > 1e-10) {
				members->TimeLevelIrr--;
				members->TimeStepIrr = static_cast<double>(pow(2, members->TimeLevelIrr));
				members->TimeBlockIrr = static_cast<ULL>(pow(2, members->TimeLevelIrr-time_block));
			}

			members->NewCurrentBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
			members->NextBlockIrr = members->CurrentBlockIrr + members->TimeBlockIrr;
		}
*/
		fprintf(binout,"PID: %d\n", members->PID);
		fprintf(binout, "Position (pc) - x:%e, y:%e, z:%e, \n", members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
		fprintf(binout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(binout, "Mass (Msol) - %e, \n", members->Mass*mass_unit);

		fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
		// fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
		// fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
		// fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
		fprintf(binout, "Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
		// fprintf(binout, "Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
		// fprintf(binout, "Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
		// fprintf(binout, "Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
		fprintf(binout, "Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
		// fprintf(binout, "Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
		// fprintf(binout, "Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
		// fprintf(binout, "Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
		fprintf(binout, "Time Steps (Myr) - irregular:%e, regular:%e \n", members->TimeStepIrr*EnzoTimeStep*1e4, members->TimeStepReg*EnzoTimeStep*1e4);
		// fprintf(binout, "Time Blocks - irregular:%llu, regular:%llu \n", members->TimeBlockIrr, members->TimeBlockReg);
		// fprintf(binout, "Current Blocks - irregular: %llu, regular:%llu \n", members->CurrentBlockIrr, members->CurrentBlockReg);
	}
	// insertNeighbors(ptclCM);
	ptclCM->clear();

	fflush(binout);
}
#endif