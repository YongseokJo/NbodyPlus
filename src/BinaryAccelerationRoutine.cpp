#include <iostream>
#include "global.h"

void NewKSInitialization(Particle* ptclI, std::vector<Particle*> &particle, double current_time);

void AddNewBinariesToList(std::vector<Particle*> &particle) {

    // add new binaries

	for (Particle *ptcl : particle) {
        // if the irregular time step is too short, check if it is binary
        if ((ptcl->TimeStepIrr<KSTime)&(ptcl->isBinary = false)) {
            ptcl->isKSCandidate(ptcl->CurrentTimeIrr + ptcl->TimeStepIrr);
            if (ptcl->isBinary) {
                NewKSInitialization(ptcl,particle,ptcl->CurrentTimeIrr);
            }
        }
    }

}

void BinaryAccelerationRoutine(double next_time) {

    for (Binary* ptclBin: BinaryList) {
        ptclBin->IntegrateBinary(next_time);
    }

}