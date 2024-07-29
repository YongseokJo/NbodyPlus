#include <cmath>
#include <iostream>
#include "defs.h"
#include "global.h"
#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

int writeParticle(std::vector<Particle*> &particle, REAL MinRegTime, int outputNum);
// int ReceiveFromEzno(std::vector<Particle*> &particle);
// int SendToEzno(std::vector<Particle*> &particle);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool RegularAccelerationRoutine(std::vector<Particle*> &particle);
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle);
void AddNewBinariesToList(std::vector<Particle*> &particle);                                                               
void BinaryAccelerationRoutine(REAL next_time, ULL next_block, std::vector<Particle*> &particle);

bool IsOutput           = false;
REAL binary_time      = 0;
REAL binary_time_prev = 0;
ULL binary_block        = 0;
REAL outputTime       = 0;
ULL NextRegTimeBlock    = 0;
int outNum              = 0;
REAL global_time_irr  = 0;
std::vector<Particle*> ComputationChain{};
#ifdef time_trace
TimeTracer _time;
#endif

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
#ifdef NSIGHT
nvtxRangePushA("RegularAccelerationRoutine");
#endif
		RegularAccelerationRoutine(particle); // do not update particles unless NNB=0
		
#ifdef NSIGHT
nvtxRangePop();
#endif


#ifdef time_trace
		_time.reg.markEnd();
		_time.reg.getDuration();
		_time.irr.markStart();
#endif
		
#ifdef NSIGHT
nvtxRangePushA("IrregularAccelerationRoutine");
#endif

		IrregularAccelerationRoutine(particle);

#ifdef NSIGHT
nvtxRangePop();
#endif

#ifdef time_trace
		_time.irr.markEnd();

		_time.irr.getDuration();
		_time.output();
#endif

		global_time = NextRegTimeBlock*time_step;

		// create output at appropriate time intervals
		if (outputTime <= global_time ) {
			writeParticle(particle, global_time, outNum++);
			outputTime += outputTimeStep;
		}

		// end if the global time exceeds the end time
		if (global_time >= 1) {
			std::cout << EnzoTimeStep << std::endl;
			return;
			//exit(EXIT_FAILURE);
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



