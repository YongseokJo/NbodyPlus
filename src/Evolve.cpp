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
void BinaryAccelerationRoutine(double next_time, std::vector<Particle*> &particle);

bool IsOutput         = false;
double binary_time = 0;
double outputTime = 0;
double NextRegTime    = 0.;
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

                  fprintf(binout, "-------------------------------------\n");                                                
                  fprintf(binout, "Evolve.cpp: starting a new while loop\n");                                                
                  fprintf(binout, "global_time = %f\n",global_time);                                                         
                  fprintf(binout, "binary_time = %f\n",binary_time);

		// It's time to compute regular force.
#ifdef time_trace
		_time.reg.markStart();
#endif

		RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
		std::cout << "Adding new binaries to list ..." << std::endl;
		AddNewBinariesToList(particle);

		if (BinaryList.size()>0) {
			std::cout << "Integrating Binaries ..." << std::endl;
			fprintf(binout, "Evolve.cpp: integrating binaries\n");
			//fprintf(binout, "# of binaries = %d\n",BinaryList.size());
			BinaryAccelerationRoutine(binary_time, particle);
		}
		
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



