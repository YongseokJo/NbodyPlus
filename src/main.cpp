#include <iostream>
//#include "defs.h"
#include "global.h"
#include "main.h"


using namespace std;

//Global Variables
int NNB; double global_time; bool debug;
std::vector<int> LevelList;

int main() {
	cout << "This Nbody+." << endl;
	std::vector<Particle*> particle{};
	//particle = new std::vector<Particle*>();
	global_time = 0.;
	debug = true;

	if (readData(particle) == FAIL) 
		fprintf(stderr, "Read Data Failed!\n");

	initializeParticle(particle);
	//createParticleChain(particle);


	// Particle should be deleted at some point

	return 0;
}
