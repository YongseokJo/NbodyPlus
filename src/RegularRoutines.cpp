#include "QueueScheduler.h"
#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

void calculateRegAccelerationOnGPU(std::unordered_set<int>& RegularList, QueueScheduler &queue_scheduler);
#ifdef MULTIMAP
void getRegularList(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList);
void updateRegularMap(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList);
#endif

void RegularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList) {

    double next_time = NextRegTimeBlock*time_step;

// Performance tracing variables (kept for backward compatibility)
#ifdef PERFORMANCETRACE
	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
#endif // performance

#ifdef MULTIMAP
    PROFILE_START(TimerID::RegularMap);
#ifdef PERFORMANCETRACE
    start_point_routine = std::chrono::high_resolution_clock::now();
#endif
    getRegularList(RegularMap, RegularList);

#ifdef PERFORMANCETRACE
    end_point_routine = std::chrono::high_resolution_clock::now();
    performance.RegularMap +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
    PROFILE_STOP(TimerID::RegularMap);
#endif // multimap

#ifdef CUDA

#ifdef NSIGHT
    nvtxRangePushA("calculateRegAccelerationOnGPU");
#endif

#ifdef DEBUG
    std::cout << "calculateRegAccelerationOnGPU starts" << std::endl;
    std::cout << "RegularList size: " << RegularList.size() << std::endl;
#endif

    calculateRegAccelerationOnGPU(RegularList, queue_scheduler);

#ifdef DEBUG
    std::cout << "calculateRegAccelerationOnGPU ended" << std::endl;
#endif

#ifdef NSIGHT
    nvtxRangePop();
#endif

#else // cuda regular routine ends // no cuda regular routine starts

    PROFILE_START(TimerID::RegularForce);
#ifdef PERFORMANCETRACE
    start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
    nvtxRangePushA("RegForce");
#endif

#ifdef DEBUG
    std::cout << "Regular force starts" << std::endl;
#endif
    // Regular force
    queue_scheduler.initialize(RegForce);
    queue_scheduler.takeQueueRegularList(RegularList);
    do
    {
        queue_scheduler.assignQueueAutoRegularList();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());
#ifdef DEBUG
    std::cout << "Regular force ended" << std::endl;
#endif

#ifdef NSIGHT
    nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
    end_point_routine = std::chrono::high_resolution_clock::now();
    performance.RegularForce +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
    PROFILE_STOP(TimerID::RegularForce);

#endif // no cuda regular routine ends

    PROFILE_START(TimerID::RegularUpdate);
#ifdef PERFORMANCETRACE
    start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
    nvtxRangePushA("RegUpdate");
#endif

#ifdef DEBUG
    std::cout << "update regular starts" << std::endl;
#endif
    // Update Regular
    queue_scheduler.initialize(RegUpdate);
    queue_scheduler.takeQueueRegularList(RegularList);
    do
    {
        queue_scheduler.assignQueueAutoRegularList();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());
#ifdef DEBUG
    std::cout << "update regular ended" << std::endl;
#endif

#ifdef NSIGHT
    nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
    end_point_routine = std::chrono::high_resolution_clock::now();
    performance.RegularUpdate +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
    PROFILE_STOP(TimerID::RegularUpdate);

#ifdef MULTIMAP
    PROFILE_START(TimerID::RegularMap);
#ifdef PERFORMANCETRACE
    start_point_routine = std::chrono::high_resolution_clock::now();
#endif
    updateRegularMap(RegularMap, RegularList);

#ifdef PERFORMANCETRACE
    end_point_routine = std::chrono::high_resolution_clock::now();
    performance.RegularMap +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
    PROFILE_STOP(TimerID::RegularMap);
#endif // multimap

}

#ifdef MULTIMAP
void getRegularList(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList) {

	// assert(RegularMap.size() == NumberOfParticle);
	if (RegularMap.size() != NumberOfParticle) { // PISN case

		int num_erased = RegularMap.size() - NumberOfParticle;
		int num = 0;

		fprintf(stdout, "PISN search (number: %d) in RegularMap...\n", num_erased);

		auto it = RegularMap.begin();
		while (it != RegularMap.end()) {
			if (particles[it->second].Mass > 0) {
				it++;
			}
			else {
				fprintf(stdout, "PISN (PID: %d) is erased in RegularMap\n", particles[it->second].PID);
				it = RegularMap.erase(it);
				num++;
				if (num == num_erased)
					break;
			}
		}
	}
	assert(RegularMap.begin()->first == NextRegTimeBlock);
	assert(RegularList.empty());

	auto it = RegularMap.begin();
	while (it->first == NextRegTimeBlock) {
		RegularList.insert(it->second);
		it = RegularMap.erase(it);
	}
	assert(RegularMap.size() + RegularList.size() == NumberOfParticle);
}

void updateRegularMap(std::multimap<ULL,int>& RegularMap, std::unordered_set<int>& RegularList) {

	assert(RegularMap.size() + RegularList.size() == NumberOfParticle);

	Particle* ptcl;
	while (!RegularList.empty()) {
		ptcl = &particles[*RegularList.begin()];
		RegularMap.insert({ptcl->CurrentBlockReg + ptcl->TimeBlockReg, ptcl->ParticleIndex});
		RegularList.erase(RegularList.begin());
	}
	assert(RegularList.empty());
	assert(RegularMap.size() == NumberOfParticle);

	NextRegTimeBlock = RegularMap.begin()->first;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
#endif