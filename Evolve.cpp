#include <cmath>
#include <iostream>
#include "defs.h"
#include "global.h"


int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
// int ReceiveFromEzno(std::vector<Particle*> &particle);
// int SendToEzno(std::vector<Particle*> &particle);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool RegularAccelerationRoutine(std::vector<Particle*> &particle);
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle);

bool IsOutput         = false;
double outputTime = 0.;
double NextRegTime    = 0.;
std::vector<Particle*> ComputationChain{};


void Evolve(std::vector<Particle*> &particle) {

	std::cout << "Evolve starting ..." << std::endl;
	int outNum = 0;
	int freq   = 0;

//	if (NNB == 0) { 
//		std::cout << "No particle to be calculated ..." << std::endl;
//		goto Communication;
//	}

	//CreateComputationChain(particle);

	while (true)
	{
		while (global_time < 1)
		{
			// It's time to compute regular force.
			RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
			IrregularAccelerationRoutine(particle);
			global_time = NextRegTime;

			// create output at appropriate time intervals
			
			if (outputTime > global_time) {
			writeParticle(particle, global_time, ++outNum);
			outputTime += outputTimeStep;
			}

		}

	//Communication:
	//	do
	//	{
	//		SendToEzno(particle);
	//		ReceiveFromEzno(particle);
	//	} while (NNB == 0);
	//	global_time = 0.;
	//	NextRegTime = 0.;
	}
}



