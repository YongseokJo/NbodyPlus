#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cassert>
#include <mpi.h>
#include "global.h"
#include "QueueScheduler.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif


void InitialAssignmentOfTasks(Queue queue, int NumTask, int TAG);
void ParticleSynchronization();
void InitializationRoutines(QueueScheduler &queue_scheduler, Worker *workers);
bool IrregularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList);
void RegularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList);

int writeParticle(double current_time, int outputNum);

#ifdef MULTIMAP
void createRegularMap(std::multimap<ULL,int>& RegularMap);
#else // no multimap
void updateNextRegTime(std::unordered_set<int>& RegularList);
#endif // multimap

#ifdef SEVN
void StellarEvolution();
#endif

Worker* workers;

void RootRoutines() {

	std::cout << "Root processor is ready." << std::endl;

	Particle* ptcl;
	TaskName task;
	int total_tasks;
	int completed_tasks=0;

	std::vector<int> newCMptcls; // by EW 2025.1.6 // unordered_set? by EW 2025.1.11
	// std::vector<int> EmptyIndex; // by EW 2025.1.7  empty slots in particles e.g., due to mergers
	// unordered_set? by EW 2025.1.11
	// merged particles & PISN will be contained here
	// new single Particle formed in Enzo can be formed in ParticleIndex of these ptcls
	// if empty, LastParticleIndex++

	bool bin_termination = false;
	bool new_binaries = false;
#ifdef MULTIMAP
	std::multimap<ULL,int> RegularMap;
#endif

	std::unordered_set<int> RegularList;

	QueueScheduler queue_scheduler;
	Queue queue;
	workers = new Worker[NumberOfWorker+1];

	for (int i=0; i<=NumberOfWorker; i++) {
		workers[i].initialize(i);
	}

#ifdef PERFORMANCETRACE
	std::chrono::high_resolution_clock::time_point start_point_whole;
	std::chrono::high_resolution_clock::time_point end_point_whole;

	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
#endif // performance

	/* Initialization */
	InitializationRoutines(queue_scheduler, workers);

	/* Actual Loop */
	{
#ifdef MULTIMAP
		createRegularMap(RegularMap);
#endif

		//ParticleSynchronization();
		while (1) {

			// create output at appropriate time intervals
			if (global_time >= outputTime) {
				writeParticle(global_time, outNum++);
				outputTime += outputTimeStep;
			}

			// end if the global time exceeds the end time
			if (global_time >= 1) {
				task=Ends;
				queue = {task, -1, -1.0};
				InitialAssignmentOfTasks(queue, NumberOfWorker, QUEUE_TAG);
				fprintf(stdout, "Simulation Done! Current Time: %e Myr\n", global_time*EnzoTimeStep*1e4);
				return;
			}

#ifdef PERFORMANCETRACE
			start_point_whole = std::chrono::high_resolution_clock::now();
#endif

#ifndef MULTIMAP

#ifdef PERFORMANCETRACE
			start_point_routine = std::chrono::high_resolution_clock::now();
#endif
#ifdef NSIGHT
			nvtxRangePushA("updateNextRegTime");
#endif
			updateNextRegTime(RegularList);
#ifdef NSIGHT
			nvtxRangePop();
#endif
#ifdef PERFORMANCETRACE
			end_point_routine = std::chrono::high_resolution_clock::now();
			performance.UpdateNextRegTime +=
				std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#endif // no multimap

			if (!IrregularRoutines(queue_scheduler, workers, RegularList))
				continue;

			RegularRoutines(queue_scheduler, workers, RegularList);

			global_time = NextRegTimeBlock*time_step;

#ifdef SEVN

#ifdef PERFORMANCETRACE
			start_point_routine = std::chrono::high_resolution_clock::now();
#endif

			if (!SEVNList.empty() && SEVNList.begin()->first <= global_time*EnzoTimeStep*1e4)
				StellarEvolution(); // Currently, evolving all the particles upto global_time

#ifdef PERFORMANCETRACE
			end_point_routine = std::chrono::high_resolution_clock::now();
			performance.StellarEvolution +=
				std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

#endif

#ifdef PERFORMANCETRACE
			end_point_whole = std::chrono::high_resolution_clock::now();
			performance.WholeRoutine +=
				std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_whole - start_point_whole).count();
#endif
			//exit(SUCCESS);
		} // While(1)
	} // Actual Loop
}

#ifdef MULTIMAP
void createRegularMap(std::multimap<ULL,int>& RegularMap) {

	assert(RegularMap.empty());

	Particle* ptcl;
	for (int i=0; i<=LastParticleIndex; i++)
	{
		ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;

		RegularMap.insert({ptcl->CurrentBlockReg + ptcl->TimeBlockReg, ptcl->ParticleIndex});
	}
	NextRegTimeBlock = RegularMap.begin()->first;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
#else // no multimap
void updateNextRegTime(std::unordered_set<int>& RegularList) {

	ULL time_tmp=0, time=block_max;
	Particle *ptcl;

	RegularList.clear();

	for (int i=0; i<=LastParticleIndex; i++)
	{
		//std::cout << i << std::endl;
		ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;
		// Next regular time step
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;

		// Find the minum regular time step
		if (time_tmp <= time) {
			//fprintf(stderr, "PID=%d, time_tme=%llu\n", ptcl->PID, time_tmp);
			if (time_tmp < time) {
				RegularList.clear();
				time = time_tmp;
			}
			//RegularList.push_back(ptcl->ParticleIndex);
			RegularList.insert(ptcl->ParticleIndex);
		}
	}
	NextRegTimeBlock = time;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
#endif // multimap