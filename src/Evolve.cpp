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
void BinaryAccelerationRoutine(double next_time, ULL next_block, std::vector<Particle*> &particle);

bool IsOutput         = false;
double binary_time = 0;
double binary_time_prev = 0;
ULL binary_block = 0;
double outputTime = 0;
ULL NextRegTimeBlock    = 0;
std::vector<Particle*> ComputationChain{};
TimeTracer _time;
int outNum = 0;

void Evolve(std::vector<Particle*> &particle) {

	std::cout << "Evolve starting ..." << std::endl;
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
#ifdef time_trace
		_time.reg.markStart();
#endif

		RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
		std::cout << "Adding new binaries to list ..." << std::endl;
		AddNewBinariesToList(particle);

		if ((BinaryList.size()>0)&(binary_time_prev != binary_time)) {
	                fprintf(binout, "-------------------------------------\n");                                                                        
                	fprintf(binout, "global_time = %e \n",global_time);                                                         
                 	fprintf(binout, "binary_time = %e \n",binary_time);

			std::cout << "Integrating Binaries ..." << std::endl;
			fprintf(binout, "Evolve.cpp: integrating binaries\n");
			fprintf(binout, "# of binaries = %d \n",int(BinaryList.size()));
			BinaryAccelerationRoutine(binary_time, binary_block, particle);
		}
																						//
#ifdef time_trace
		_time.reg.markEnd();
		_time.irr.markStart();
#endif

		IrregularAccelerationRoutine(particle);

#ifdef time_trace
		_time.irr.markEnd();

		_time.reg.getDuration();
		_time.irr.getDuration();
		_time.output();
#endif

		global_time = NextRegTimeBlock*time_step;
		binary_time_prev = binary_time;
		binary_time = FirstComputation->CurrentTimeIrr + FirstComputation->TimeStepIrr;
		binary_block = FirstComputation->CurrentBlockIrr + FirstComputation->TimeBlockIrr;

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



