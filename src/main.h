#include "Particle/Particle.h"

int readData(std::vector<Particle*> &particle);
int createComputationChain(std::vector<Particle*> &particle);
int initializeParticle(std::vector<Particle*> &particle);
void Evolve(std::vector<Particle*> &particle);
int Parser(int argc, char* argv[]);

