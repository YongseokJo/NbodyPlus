#include <iostream>
#include "global.h"
#include "defs.h"

void NewKSInitialization(Particle* ptclI, std::vector<Particle*> &particle, double current_time);
void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle, double current_time);

void AddNewBinariesToList(std::vector<Particle*> &particle) {

    // add new binaries
    //
    for (Particle *ptcl : particle) {
    // if the irregular time step is too short, check if it is binary
    	if ((ptcl->TimeStepIrr<KSTime)&&(ptcl->isBinary == false)) {
    		ptcl->isKSCandidate(ptcl->CurrentTimeIrr + ptcl->TimeStepIrr);
    		if (ptcl->isBinary) {
		        std::cout << "AddNewBinaries ... new binary pair found" << std::endl;
			fprintf(binout, "BinaryAccelerationRoutine.cpp: new binary particle found!\n");
  			NewKSInitialization(ptcl,particle,ptcl->CurrentTimeIrr);
			std::cout << "New initialization finished ..." << std::endl;

		}
	}
    }

   fprintf(binout,"\n After binary addition, the number of particles are... %d \n",int(particle.size()));

}

void BinaryAccelerationRoutine(double next_time, std::vector<Particle*> &particle) {

	int count;
        int bincount = 0;

	count = 0;

	if (next_time == 0) {
		return;
	}

	for (Binary* ptclBin: BinaryList) {

	    	ptclBin->KSIntegration(next_time, bincount);

		count += 1;

		fprintf(binout, "\nBinaryAccelerationRoutine.cpp: After KS Integration of %dth binary....\n", count);
		fprintf(binout, "The ID of ith particle is %d \n",ptclBin->ptclCM->BinaryParticleI->PID);
		fprintf(binout, "The ID of ith particle is %d \n",ptclBin->ptclCM->BinaryParticleJ->PID);
		fflush(binout);	

		if (bincount>0) {
			std::cout << "Integrating Binary ..." << std::endl;

			fprintf(binout, "KS coordinates - u1:%e, u2:%e, u3:%e, u4:%e\n", ptclBin->u[0], ptclBin->u[1], ptclBin->u[2], ptclBin->u[3]);
			fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e\n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
			fprintf(binout, "KS coordinates - udot1:%e, udot2:%e, udot3:%e, udot4:%e\n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
			fprintf(binout, "Other important KS variables - r:%e, h:%e, gamma: %e, tau:%e, step:%e, currentTime: %e \n", ptclBin->r, ptclBin->h, ptclBin->gamma, ptclBin->dTau, ptclBin->TimeStep, ptclBin->CurrentTime);
			fprintf(binout, "loop number = %d \n", bincount);
			fflush(binout);
		}

		// check for termination
		
		if ((ptclBin->r>(ptclBin->r0*2.0)) || (ptclBin->TimeStep > 2.0*KSTime)) {
	        	std::cout << "Terminating Binary ..." << std::endl;
			fprintf(binout, "Terminating Binary at time : %e \n", next_time);
			KSTermination(ptclBin->ptclCM, particle, next_time);
		}

    }
}
