#ifdef FEWBODY
#include "../global.h"
#include "Queue.h"
#include <unordered_set>

void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void CalculateAcceleration01MPI(Particle* ptcl1, std::vector<double>& accIrrReg_send);
void CalculateAcceleration23MPI(Particle* ptcl1, std::vector<double>& accIrrReg_send);

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

		if (members->Mass < 0.0)
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

		// members->RadiusOfNeighbor = ACRadius*ACRadius; // added by EW 2025.1.16
		members->RadiusOfNeighbor = ptclCM->RadiusOfNeighbor; // modified by EW 2025.1.30

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
			// members->calculateTimeStepIrr2();
			members->calculateTimeStepIrr();

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

void FBTerminationRoot(Particle* ptclCM) {

	fprintf(binout,"--------------------------------------\n");
	fprintf(binout,"In FBTermination.cpp... (CM PID: %d)\n", ptclCM->PID);
	fprintf(binout, "CurrentTimeIrr of ptclCM (Myr): %e\n", ptclCM->CurrentTimeIrr*EnzoTimeStep*1e4);
	fprintf(binout, "CurrentTimeIrr of the first member (Myr): %e\n", particles[ptclCM->Members[0]].CurrentTimeIrr*EnzoTimeStep*1e4);
	fprintf(binout, "N_member: %d\n", ptclCM->NumberOfMember);

	NumberOfParticle--; // CM particle should be inactive by EW 2025.1.20

	Queue queue{FBTerminationMPI, ptclCM->ParticleIndex, -1.0};
	std::vector<MPI_Request> requests;

	for (int rank = 1; rank < NumberOfProcessor; rank++) {
		MPI_Request req;
		MPI_Isend(&queue, 1, QueueType, rank, QUEUE_TAG, MPI_COMM_WORLD, &req);
		requests.push_back(req);
	}
	MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
	
	Particle* members;
	for (int i=0; i<ptclCM->NumberOfMember; i++) {
		members = &particles[ptclCM->Members[i]];

		if (members->Mass < 0.0)
			continue;

		std::vector<double> accIrrReg_send(Dim * HERMITE_ORDER * 2); // all the members are initialized as 0.0 by EW 2025/7/16
		std::vector<double> accIrrReg_recv(Dim * HERMITE_ORDER * 2); // all the members are initialized as 0.0 by EW 2025/7/16
		MPI_Request request;
		MPI_Ireduce(accIrrReg_send.data(), accIrrReg_recv.data(), Dim * HERMITE_ORDER * 2, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD, &request);

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

		// members->RadiusOfNeighbor = ptclCM->RadiusOfNeighbor; // modified by EW 2025.1.30 // This was done in RootRoutines.cpp by EW 2025.7.16

		members->NumberOfNeighbor = ptclCM->NumberOfNeighbor;
		std::memcpy(members->Neighbors, ptclCM->Neighbors, sizeof(int) * ptclCM->NumberOfNeighbor);
		for (int j = 0; j < ptclCM->NumberOfMember; j++) {
			Particle* members_members = &particles[ptclCM->Members[j]];
			if (members_members->Mass > 0.0 && members_members->PID != members->PID)
				members->Neighbors[members->NumberOfNeighbor++] = members_members->ParticleIndex;
		}

		MPI_Wait(&request, MPI_STATUS_IGNORE);
		for (int dim = 0; dim < Dim; dim++) {
			for (int order = 0; order < HERMITE_ORDER; order++) {
				members->a_irr[dim][order] = accIrrReg_recv[dim * HERMITE_ORDER + order];
				members->a_reg[dim][order] = accIrrReg_recv[dim * HERMITE_ORDER + order + Dim * HERMITE_ORDER];
				members->a_tot[dim][order] = members->a_irr[dim][order] + members->a_reg[dim][order];
			}
		}

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
			// members->calculateTimeStepIrr2();
			members->calculateTimeStepIrr();

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

void FBTerminationWorker(Particle* ptclCM) {
	
	Particle* members;
	for (int i=0; i<ptclCM->NumberOfMember; i++) {
		members = &particles[ptclCM->Members[i]];

		if (members->Mass < 0.0)
			continue;

		std::vector<double> accIrrReg_send(Dim * HERMITE_ORDER * 2); // all the members are initialized as 0.0 by EW 2025/7/16
		std::vector<double> accIrrReg_recv(Dim * HERMITE_ORDER * 2); // all the members are initialized as 0.0 by EW 2025/7/16

		CalculateAcceleration01MPI(members, accIrrReg_send);
		CalculateAcceleration23MPI(members, accIrrReg_send);

		MPI_Request request;
		MPI_Ireduce(accIrrReg_send.data(), accIrrReg_recv.data(), Dim * HERMITE_ORDER * 2, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD, &request);
		
		MPI_Wait(&request, MPI_STATUS_IGNORE);
	}
}

void CalculateAcceleration01MPI(Particle* ptcl1, std::vector<double>& accIrrReg_send) {

	double x[Dim], v[Dim];
	double m_r3;
	double v2;
	double r2 = 0;
	double vx = 0;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]    = 0.;
		v[dim]    = 0.;
	}

	int size = global_variable->LastParticleIndex + 1;

	int start = (size * (MyRank - 1)) / NumberOfWorker;
	int end   = (size * MyRank) / NumberOfWorker;

	Particle *ptcl2;
	for (int i = start; i < end; i++) {
		ptcl2 = &particles[i];

		if (!ptcl2->isActive || ptcl1->PID == ptcl2->PID) {
			continue;
		}

		r2 = 0;
		vx = 0;
		v2 = 0;

		// updated the predicted positions and velocities just in case
		// if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
		//ptcl2->predictParticleSecondOrder(newTime);
		for (int dim=0; dim<Dim; dim++) {
			x[dim] = ptcl2->Position[dim] - ptcl1->Position[dim];
			v[dim] = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
			r2    += x[dim]*x[dim];
			vx    += v[dim]*x[dim];
			v2    += v[dim]*v[dim];
		}

		m_r3 = ptcl2->Mass/r2/sqrt(r2); 

		if (r2 < ptcl1->RadiusOfNeighbor) {
			for (int dim=0; dim<Dim; dim++) {
				accIrrReg_send[dim * HERMITE_ORDER + 0] += m_r3*x[dim];
				accIrrReg_send[dim * HERMITE_ORDER + 1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				// Calculate 0th and 1st derivatives of acceleration
				accIrrReg_send[dim * HERMITE_ORDER + 0 + Dim * HERMITE_ORDER] += m_r3*x[dim];
				accIrrReg_send[dim * HERMITE_ORDER + 1 + Dim * HERMITE_ORDER] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
			}
		} // endfor dim
	} // endfor ptcl2
	return;
}

void CalculateAcceleration23MPI(Particle* ptcl1, std::vector<double>& accIrrReg_send) {

	double x[Dim], v[Dim], a21[Dim], a21dot[Dim], a1[Dim], a2[Dim], a1dot[Dim], a2dot[Dim];
	double a, b, c;
	double rdf_r2, vdf_r2, rdfdot_r2, v2, r2, r3, vr, m_r3;
	double adot2, adot3;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]      = 0.;
		v[dim]      = 0.;
		a21[dim]    = 0.;
		a21dot[dim] = 0.;
		a1[dim]     = ptcl1->a_tot[dim][0];
		a1dot[dim]  = ptcl1->a_tot[dim][1];
	}

	int size = global_variable->LastParticleIndex + 1;

	int start = (size * (MyRank - 1)) / NumberOfWorker;
	int end   = (size * MyRank) / NumberOfWorker;

	Particle *ptcl2;
	for (int i = start; i <= end; i++) {
		ptcl2 = &particles[i];

		if (!ptcl2->isActive || ptcl1->PID == ptcl2->PID) {
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
		for (int dim=0; dim<Dim; dim++) {
			a2[dim]	   = accIrrReg_send[dim * HERMITE_ORDER + 0 + dim * HERMITE_ORDER + 0 + Dim * HERMITE_ORDER];
			a2dot[dim] = accIrrReg_send[dim * HERMITE_ORDER + 1 + dim * HERMITE_ORDER + 1 + Dim * HERMITE_ORDER];
			x[dim]     = ptcl2->Position[dim] - ptcl1->Position[dim];
			v[dim]     = ptcl2->Velocity[dim] - ptcl1->Velocity[dim];
			r2        += x[dim]*x[dim];
			vr        += v[dim]*x[dim];
			v2        += v[dim]*v[dim];
		}

		r3   = r2*sqrt(r2);
		m_r3 = ptcl2->Mass/r3; 

		for (int dim=0; dim<Dim; dim++) {
			a21[dim]    = m_r3*x[dim];
			a21dot[dim] = m_r3*(v[dim] - 3*x[dim]*vr/r2);
			rdf_r2     += x[dim]*(a1[dim]-a2[dim])/r2;
			vdf_r2     += v[dim]*(a1[dim]-a2[dim])/r2;
			rdfdot_r2  += x[dim]*(a1dot[dim]-a2dot[dim])/r2;
		}

		a = vr/r2;
		b = v2/r2 + rdf_r2 + a*a;
		c = 3*vdf_r2 + rdfdot_r2 + a*(3*b-4*a*a);

		if (r2 < ptcl1->RadiusOfNeighbor) {
			for (int dim=0; dim<Dim; dim++) {
				adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				accIrrReg_send[dim * HERMITE_ORDER + 2] += adot2;
				accIrrReg_send[dim * HERMITE_ORDER + 3] += adot3;
			}
		}
		else {
			for (int dim=0; dim<Dim; dim++) {
				adot2 = -ptcl2->Mass*(a1[dim]-a2[dim])/r3-6*a*a21dot[dim]-3*b*a21[dim];
				adot3 = -ptcl2->Mass*(a1dot[dim]-a2dot[dim])/r3-9*a*adot2-9*b*a21dot[dim]-3*c*a21[dim];
				accIrrReg_send[dim * HERMITE_ORDER + 2 + Dim * HERMITE_ORDER] += adot2;
				accIrrReg_send[dim * HERMITE_ORDER + 3 + Dim * HERMITE_ORDER] += adot3;
			}
		} // endfor if
	} //endfor ptcl2
}
#endif