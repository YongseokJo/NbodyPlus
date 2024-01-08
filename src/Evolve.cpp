#include "Particle/Particle.h"


void SortComputationChain(std::vector<Particle*> &particle);

void Evolve(std::vector<Particle*> &particle) {

	std::cout << "Evolve starting ..." << std::endl;
	Particle* ptcl;
	double ComputingTime, EvolvedTime; // time where computing is actaully happening

	std::cout << particle.size() << std::endl;
	std::cout << LevelList[0]	<< std::endl;
	ptcl = particle[LevelList[0]];

	// This part can be parallelized.
	ComputingTime = ptcl->CurrentTimeIrr + ptcl->TimeStepIrr;

	while (true) {
		// check if the particle is subject to calculation this level.
		fprintf(stderr, "in loop\n");
		for (int i=0; i<NNB; i++) {
			ptcl = particle[i];
			EvolvedTime = ptcl->CurrentTimeIrr+ptcl->TimeStepIrr; 
			fprintf(stderr, "PID=%d, com t= %e, evl t = %e\n",
					ptcl->getPID(),ComputingTime, EvolvedTime);
			//if (ComputingTime >= EvolvedTime) { // Condition should be added more.
			if (ptcl->checkNeighborForEvolution())
				break;
			if (i == NNB-1) {
				//std::cout << "Warning::Nothing to evolve!!!\n" <<std::flush;
				std::cout << "Error::Nothing to evolve!!!\n" <<std::flush;
				exit(EXIT_FAILURE); 
			}
		}

		fprintf(stderr, "Particle ID=%d\n",ptcl->getPID());
		fprintf(stderr, "Current Time = %e, time step = %e\n",
			 	ptcl->CurrentTimeIrr, ptcl->TimeStepIrr);

		if (ptcl->CurrentTimeIrr == ptcl->CurrentTimeReg + ptcl->TimeStepReg) {
			std::cout <<  "Regular force calculating...\n" << std::flush;
				ptcl->calculateRegForce(particle); // this only does acceleration computation.
		}
		std::cout << "Irregular force calculating...\n" << std::flush;
		ptcl->calculateIrrForce(); // this includes particle and time advances.
		//ptcl = ptcl->NextParticle;
	}
	//SortComputationChain(particle);
}
