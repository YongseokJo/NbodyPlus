#include <iostream>
#include "global.h"
#include "defs.h"

void NewKSInitialization(Particle* ptclI, std::vector<Particle*> &particle, double current_time);
void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle);

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

}

void BinaryAccelerationRoutine(double next_time, std::vector<Particle*> &particle) {

	int count;

	count = 0;

	for (Binary* ptclBin: BinaryList) {

	        ptclBin->IntegrateBinary(next_time);

		count += 1;
		std::cout << "Integrating Binary ..." << std::endl;

		fprintf(binout, "BinaryAccelerationRoutine.cpp: After KS Integration of %dth binary....\n", count);
		fprintf(binout, "KS coordinates - u1:%f, u2:%f, u3:%f, u4:%f/n", ptclBin->u[0], ptclBin->u[1], ptclBin->u[2], ptclBin->u[3]);
		fprintf(binout, "KS coordinates - udot1:%f, udot2:%f, udot3:%f, udot4:%f/n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
		fprintf(binout, "KS coordinates - udot1:%f, udot2:%f, udot3:%f, udot4:%f/n", ptclBin->udot[0], ptclBin->udot[1], ptclBin->udot[2], ptclBin->udot[3]);
		fprintf(binout, "Other important KS variables - r:%f, h:%f, gamma: %f, tau:%f, step:%f/n", ptclBin->r, ptclBin->h, ptclBin->gamma, ptclBin->dTau, ptclBin->TimeStep);


		// check for termination
		
		if ((ptclBin->r>(ptclBin->r0*2)) || (ptclBin->TimeStep > 2*KSTime)) {
	                std::cout << "Terminating Binary ..." << std::endl;
			KSTermination(ptclBin->ptclCM, particle);
		}

    	}
}
