#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#include <iostream>
#include <chrono>


struct Performance {

	long WholeRoutine = 0;

	long IrregularForce = 0;
	long IrregularUpdate = 0;

	long FewBodyTermination = 0;
	long FewBodySearch = 0;
	long FewBodyInitialization = 0;

	long RegularForce = 0; // use this only if CUDA is not defined
	long RegularSendAllParticlesToGPU = 0;
	long RegularGPU = 0;
	long RegularAdjust = 0;
	long RegularUpdate = 0;

	long SkipListCreate = 0;
	long SkipListUpdate = 0;
	long UpdateNextRegTime = 0;
	long UpdateNextRegTimeFB = 0;

#ifdef SEVN
	long StellarEvolution = 0;
#endif

};

#endif

