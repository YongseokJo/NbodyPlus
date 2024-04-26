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
void AddNewBinariesToList(std::vector<Particle*> &particle);
void BinaryAccelerationRoutine(double next_time);

double binary_time = 0;
bool IsOutput         = false;
double outputTime = 0;
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
	writeParticle(particle, global_time, outNum++);
	outputTime = outputTimeStep;

	while (true) {

		// It's time to compute regular force.
		AddNewBinariesToList(particle);

		if (BinaryList.size()>0) {
			BinaryAccelerationRoutine(binary_time);
		}

		RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
		IrregularAccelerationRoutine(particle);
		global_time = NextRegTime;
		binary_time = FirstComputation->CurrentTimeIrr + FirstComputation->TimeStepIrr;

		// create output at appropriate time intervals
		if (outputTime <= global_time ) {
			writeParticle(particle, global_time, outNum++);
			outputTime += outputTimeStep;
		}

		// end if the global time exceeds the end time
		if (global_time >= 1) {
			std::cout << EnzoTimeStep << std::endl;
			exit(EXIT_FAILURE);
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



