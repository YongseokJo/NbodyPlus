#include "queue_scheduler.h"
#include "particle_data.h"
#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

// ============================================================================
// Regular force routines - SoA Integration Notes
// ============================================================================
// The force calculation path uses:
// 1. calculateRegAccelerationOnGPU() - GPU path (Phase 3 SoA integration)
// 2. Particle::compute_acceleration_reg() - CPU path (SoA helpers in Phase 4)
// 3. Particle::update_regular_particle_cuda() - Post-GPU neighbor update
//
// Current data flow:
// - particles[] (AoS) remains source of truth for routine orchestration
// - particle_data (SoA) is synced for GPU transfers (Phase 3)
// - Force calculation uses SoA helpers when beneficial
//
// Future optimization: Add sync_from_particle()/sync_to_particle() at
// RegularRoutines entry/exit to enable pure SoA force calculation.
// ============================================================================

void calculateRegAccelerationOnGPU(std::unordered_set<int>& RegularList, QueueScheduler &queue_scheduler);
#ifdef MULTIMAP
void getRegularList(std::multimap<ull_t,int>& RegularMap, std::unordered_set<int>& RegularList);
void updateRegularMap(std::multimap<ull_t,int>& RegularMap, std::unordered_set<int>& RegularList);
#endif

void RegularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList) {

    double next_time = next_reg_time_block*time_step;

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
    nvtxRangePushA("TASK_REG_FORCE");
#endif

#ifdef DEBUG
    std::cout << "Regular force starts" << std::endl;
#endif
    // Regular force
    queue_scheduler.initialize(TASK_REG_FORCE);
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
    nvtxRangePushA("TASK_REG_UPDATE");
#endif

#ifdef DEBUG
    std::cout << "update regular starts" << std::endl;
#endif
    // Update Regular
    queue_scheduler.initialize(TASK_REG_UPDATE);
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
void getRegularList(std::multimap<ull_t,int>& RegularMap, std::unordered_set<int>& RegularList) {

	// assert(RegularMap.size() == num_particles);
	if (RegularMap.size() != num_particles) { // PISN case

		int num_erased = RegularMap.size() - num_particles;
		int num = 0;

		fprintf(stdout, "PISN search (number: %d) in RegularMap...\n", num_erased);

		auto it = RegularMap.begin();
		while (it != RegularMap.end()) {
			if (particles[it->second].mass > 0) {
				it++;
			}
			else {
				fprintf(stdout, "PISN (PID: %d) is erased in RegularMap\n", particles[it->second].pid);
				it = RegularMap.erase(it);
				num++;
				if (num == num_erased)
					break;
			}
		}
	}
	assert(RegularMap.begin()->first == next_reg_time_block);
	assert(RegularList.empty());

	auto it = RegularMap.begin();
	while (it->first == next_reg_time_block) {
		RegularList.insert(it->second);
		it = RegularMap.erase(it);
	}
	assert(RegularMap.size() + RegularList.size() == num_particles);
}

void updateRegularMap(std::multimap<ull_t,int>& RegularMap, std::unordered_set<int>& RegularList) {

	assert(RegularMap.size() + RegularList.size() == num_particles);

	Particle* ptcl;
	while (!RegularList.empty()) {
		ptcl = &particles[*RegularList.begin()];
		RegularMap.insert({ptcl->current_block_reg + ptcl->time_block_reg, ptcl->particle_index});
		RegularList.erase(RegularList.begin());
	}
	assert(RegularList.empty());
	assert(RegularMap.size() == num_particles);

	next_reg_time_block = RegularMap.begin()->first;
	g_state->next_reg_time_block = next_reg_time_block;
}
#endif