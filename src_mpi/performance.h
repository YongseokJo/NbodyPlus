#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#include <iostream>
#include <chrono>


struct Performance {
	//std::chrono::milliseconds IrregularForce{0};
	long IrregularForce = 0;
	long IrregularRoutine = 0;
	long RegularRoutine = 0;
	// std::chrono::milliseconds RegularForce{0};
};

#endif

