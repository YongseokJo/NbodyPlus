#include "Particle.h"


void Particle::updateEvolveParticle(std::vector<Particle*> &list, double MinRegTime) {
		if ((MinRegTime >= this->CurrentTimeIrr+this->TimeStepIrr)
				&& (this->checkNeighborForEvolution())) {
			return;
		} else {
			int i=0;
			for (Particle* ptcl:list) {
				if (ptcl == this)
					break;
				i++;
			}
			list.erase(list.begin()+i);
		}
}


