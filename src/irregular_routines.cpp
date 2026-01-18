#include "queue_scheduler.h"
#include "skip_list.h"
#include "particle_data.h"
#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

// ============================================================================
// Irregular force routines - SoA Integration Notes
// ============================================================================
// The irregular force loop uses:
// 1. Particle::compute_acceleration_irr() - Uses SoA helpers from Phase 4
// 2. Particle::update_particle() - State update (inline in particle.h)
// 3. Skip list for time stepping - Uses particles[] indices
//
// Current data flow:
// - particles[] (AoS) remains source of truth for:
//   - Skip list operations (next_block_irr, current_block_irr)
//   - FewBody group management (is_cm_particle, members[], etc.)
//   - Active particle filtering
// - Force calculation (compute_acceleration_irr) uses SoA helpers internally
//
// The SoA infrastructure enables future optimization where the inner force
// loop can operate on contiguous SoA arrays for better cache efficiency.
// ============================================================================


bool createSkipList(SkipList *skiplist);
bool updateSkipList(SkipList *skiplist, int ptcl_id);

void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls, std::unordered_map<int, int>& existing, std::unordered_map<int, int>& terminated);
void FBTermination(Particle* ptclCM);
void Merge(Particle* p1, Particle* p2);

#ifdef SEVN_BINARY
bool makeSEVNBinary(Particle* ptclCM);
void deleteSEVNBinary(Particle* ptclCM);
void BinaryEvolution(Particle* ptclCM);
#endif

bool IrregularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList) {

    int max_level = 20;
    double prob = 0.5;
    SkipList *skiplist;
    Node* ThisLevelNode;
    double current_time_irr=0;
    double next_time=0;
    Worker* worker;

    Particle* ptcl;
    Queue queue;
#ifdef FEWBODY
    bool bin_termination = false;
    bool new_binaries = false;
    std::vector<int> newCMptcls; // by EW 2025.1.6 // unordered_set? by EW 2025.1.11
#endif

// Performance tracing variables (kept for backward compatibility)
#ifdef PERFORMANCETRACE
	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
#ifdef SEVN_BINARY
	std::chrono::high_resolution_clock::time_point start_point_BSE;
	std::chrono::high_resolution_clock::time_point end_point_BSE;
#endif
#ifdef MULTIMAP
	std::chrono::high_resolution_clock::time_point start_point_map;
	std::chrono::high_resolution_clock::time_point end_point_map;
#endif
#endif // performance
        
    PROFILE_START(TimerID::SkipListCreate);
#ifdef PERFORMANCETRACE
    start_point_routine = std::chrono::high_resolution_clock::now();
#endif
#ifdef NSIGHT
    nvtxRangePushA("createSkipList");
#endif
    skiplist = new SkipList(max_level, prob);
    if (createSkipList(skiplist) == false)
        fprintf(stderr, "There are no irregular particles!\nBut is it really happening? check skiplist->display()\n");
#ifdef NSIGHT
    nvtxRangePop();
#endif
#ifdef PERFORMANCETRACE
    end_point_routine = std::chrono::high_resolution_clock::now();
    performance.SkipListCreate +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
    PROFILE_STOP(TimerID::SkipListCreate);

    // Irregular
    while ( skiplist->getFirstNode() != nullptr) {

        ThisLevelNode = skiplist->getFirstNode();
        ThisLevelNode->ParticleList.erase(
            std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
                [](int i) {
                return !particles[i].is_active;
                }
            ),
            ThisLevelNode->ParticleList.end()
        );
        if (ThisLevelNode->ParticleList.size() == 0) {
            skiplist->deleteFirstNode();
            continue;
        }

        next_time     = particles[ThisLevelNode->ParticleList[0]].current_time_irr\
                                    + particles[ThisLevelNode->ParticleList[0]].time_step_irr;

        PROFILE_START(TimerID::IrregularForce);
#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef DEBUG
        // print out particlelist
        fprintf(stdout, "(IRR_FORCE) next_time: %e Myr\n", next_time*enzo_time_step*1e4);
        // /*
        fprintf(stdout, "PID: %d. CurrentTimeIrr: %e Myr, TimeStepIrr: %e Myr\n", 
                    particles[ThisLevelNode->ParticleList[0]].pid, 
                    particles[ThisLevelNode->ParticleList[0]].current_time_irr*enzo_time_step*1e4, 
                    particles[ThisLevelNode->ParticleList[0]].time_step_irr*enzo_time_step*1e4);

        // fprintf(stdout, "PID (%d) = ", ThisLevelNode->ParticleList.size());
        for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            // fprintf(stdout, "%d, ", ptcl->pid);
            fprintf(stdout, "PID: %d. %e Myr, %e Myr\n", 
                    ptcl->pid,
                    ptcl->current_time_irr*enzo_time_step*1e4,
                    ptcl->time_step_irr*enzo_time_step*1e4);
        }
        fprintf(stdout, "\n");
        // fflush(stdout);
        // */
#endif

        // Irregular Force
#ifdef FEWBODY
#ifdef NSIGHT
        nvtxRangePushA("IrregularForce");
#endif
#ifdef DEBUG
        std::cout << "Irr force starts" << std::endl;
#endif
// /*
        int cm_pid;
        queue_scheduler.initializeIrr(TASK_IRR_FORCE, next_time, ThisLevelNode->ParticleList);
        auto iter = queue_scheduler.CMPtcls.begin();
#ifdef SEVN_BINARY
        std::vector<int> CMPtclsForSEVN;
#endif
        // Pre-post receives for async pattern (Phase 12)
        queue_scheduler.postAllReceives();

        do {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAsync();  // Async sends (Phase 12)

            // Start async window timer (Phase 12)
            PROFILE_START(TimerID::AsyncWindow);
#ifdef DEBUG
            int cm_tasks_dispatched = 0;  // Track local work during async window (Phase 12)
#endif
            do {
                worker = queue_scheduler.testQueueAsync();  // Non-blocking test (Phase 12)
                // if there's any CMPtcl - this is our local work during async window
                if (queue_scheduler.CMPtcls.size() > 0)
                {
                    PROFILE_START(TimerID::OverlapWork);
                    // check if there's any CM ptcl ready to go for SDAR
                    if (iter == queue_scheduler.CMPtcls.end())
                        iter = queue_scheduler.CMPtcls.begin();
                    cm_pid = *(iter);
                    ptcl = &particles[cm_pid];
#ifdef DEBUG
                    fprintf(stdout, "(ASYNC_CM) Dispatching CM particle PID: %d to worker %d\n",
                            ptcl->pid, cm_particle_worker_map[cm_pid]);
                    cm_tasks_dispatched++;
#endif
                    /*
                    for (int j = 0; j < ptcl->num_neighbors; j++)
                    {
                        // if (particles[ptcl->neighbors[j]].is_up_to_date == false) // original code
                        if (particles[ptcl->neighbors[j]].is_active && !particles[ptcl->neighbors[j]].is_up_to_date) // modified by EW 2025.2.26
                        {
                            iter++;
                            goto skip_to_next;
                        }
                    }
                    */
                    queue.task = TASK_AR_INTEGRATION;
                    queue.pid = cm_pid;
                    queue.next_time = next_time;
                    workers[cm_particle_worker_map[cm_pid]].add_queue(queue);
                    queue_scheduler.assignWorker(&workers[cm_particle_worker_map[cm_pid]]);
                    iter = queue_scheduler.CMPtcls.erase(iter);
#ifdef SEVN_BINARY
                    if (ptcl->binary_evolution != nullptr)
                        CMPtclsForSEVN.push_back(cm_pid);
#endif
                    PROFILE_STOP(TimerID::OverlapWork);
                // skip_to_next:;
                }
                if (worker != nullptr)
                {
#ifdef DEBUG
                    fprintf(stdout, "(ASYNC_COMPLETE) Worker %d completed task for PID: %d\n",
                            worker->rank, worker->getCurrentQueue()->pid);
#endif
                    queue_scheduler.callbackAsync(worker);  // Async callback (Phase 12)
                }
            } while (worker == nullptr);

            // Stop async window timer (Phase 12)
            PROFILE_STOP(TimerID::AsyncWindow);
#ifdef DEBUG
            if (cm_tasks_dispatched > 0) {
                fprintf(stdout, "(ASYNC_OVERLAP) Dispatched %d CM tasks during async window\n",
                        cm_tasks_dispatched);
            }
#endif
        } while (queue_scheduler.isComplete());

#ifdef SEVN_BINARY
#ifdef PERFORMANCETRACE
        start_point_BSE = std::chrono::high_resolution_clock::now();
#endif
        for (int cm_index: CMPtclsForSEVN) {
            ptcl = &particles[cm_index];
            BinaryEvolution(ptcl);
        }
#ifdef PERFORMANCETRACE
        end_point_BSE = std::chrono::high_resolution_clock::now();
        performance.BinaryStellarEvolution +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_BSE - start_point_BSE).count();
#endif
#endif
// */
/*
        queue_scheduler.initialize(TASK_IRR_FORCE, next_time);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());

        for (int ptcl_id : ThisLevelNode->ParticleList)
        {
            ptcl = &particles[ptcl_id];

            if (ptcl->is_cm_particle) {
                int rank = cm_particle_worker_map[ptcl->particle_index];
                queue.task = TASK_AR_INTEGRATION;
                queue.pid = ptcl->particle_index;
                queue.next_time = next_time;
                workers[rank].add_queue(queue);
                workers[rank].run_queue();
                workers[rank].callback();
            }
        }
*/
#ifdef DEBUG
            std::cout << "Irregular Force done" << std::endl;
#endif
#ifdef NSIGHT
        nvtxRangePop();
#endif
#else
        queue_scheduler.initialize(TASK_IRR_FORCE, next_time);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        queue_scheduler.postAllReceives();  // Pre-post receives (Phase 12)
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAsync();  // Async sends (Phase 12)
            Worker* completed_worker = queue_scheduler.waitQueueAsync();  // Async wait (Phase 12)
            if (completed_worker != nullptr) {
                queue_scheduler.callbackAsync(completed_worker);  // Async callback (Phase 12)
            }
        } while (queue_scheduler.isComplete());
#endif

#ifdef PERFORMANCETRACE
        end_point_routine = std::chrono::high_resolution_clock::now();
        performance.IrregularForce +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
        PROFILE_STOP(TimerID::IrregularForce);

        //ParticleSynchronization();


        PROFILE_START(TimerID::IrregularUpdate);
#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
        nvtxRangePushA("IrregularUpdate");
#endif
        /* // For this kind of simple work, using queue_scheduler is slower
        queue_scheduler.initialize(TASK_IRR_UPDATE);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());
        */
        for (int ptcl_id : ThisLevelNode->ParticleList) {
            ptcl = &particles[ptcl_id];

            if (ptcl->num_neighbors != 0) // IAR modified
                ptcl->update_particle();
            ptcl->current_block_irr = ptcl->new_current_block_irr;
            ptcl->current_time_irr  = ptcl->current_block_irr*time_step;
        }
#ifdef DEBUG
        for (int i: ThisLevelNode->ParticleList) {
            ptcl = &particles[i];
            if (ptcl->current_time_irr != next_time) {
                fprintf(stdout, "TASK_ERROR! PID: %d, CurrentTimeIrr: %e Myr, next_time: %e Myr\n", ptcl->pid, ptcl->current_time_irr*enzo_time_step*1e4, next_time*enzo_time_step*1e4);
                assert(ptcl->current_time_irr == next_time);
            }
        }
            std::cout << "Irregular update done" << std::endl;
#endif
#ifdef NSIGHT
        nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
        end_point_routine = std::chrono::high_resolution_clock::now();
        performance.IrregularUpdate +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
        PROFILE_STOP(TimerID::IrregularUpdate);
        int OriginalParticleListSize;
#ifdef FEWBODY

        PROFILE_START(TimerID::FewBodyTermination);
#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
        nvtxRangePushA("FewBodyTermination");
#endif
        OriginalParticleListSize = ThisLevelNode->ParticleList.size();
        for (int i=0; i<OriginalParticleListSize; i++) {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            if (ptcl->get_binary_interrupt_state() == BinaryInterruptState::merger ||
                ptcl->get_binary_interrupt_state() == BinaryInterruptState::terminated) {

                assert(ptcl->is_cm_particle); // for debugging by EW 2025.1.20
            
                if (ptcl->get_binary_interrupt_state() == BinaryInterruptState::merger) {
                    if (ptcl->num_members == 2) { // binary merger

                        Particle* donor = &particles[ptcl->members[0]];
                        Particle* accretor = &particles[ptcl->members[1]];

                        Merge(donor, accretor);
#ifdef SEVN
                        if (donor->stellar_evolution != nullptr && accretor->stellar_evolution != nullptr) {
                            fprintf(stdout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->pid, donor->mass*mass_unit, donor->stellar_evolution->get_zams());
                            fprintf(stdout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->pid, accretor->mass*mass_unit, accretor->stellar_evolution->get_zams());
                            fprintf(stdout, "Donor: amiempty(): %d\n", donor->stellar_evolution->amiempty());
                            fprintf(stdout, "Accretor: amiempty(): %d\n", accretor->stellar_evolution->amiempty());
                            fflush(stdout);
                        }
#endif
                    }
                    else { // from NewFBInitialization3

                        assert(ptcl->num_members > 2); // for debugging by EW 2025.1.20

                        Particle* donor;
                        Particle* accretor;

                        for (int j=0; j<ptcl->num_members; j++) {
                            if (particles[ptcl->members[j]].get_binary_interrupt_state() == BinaryInterruptState::collision) {
                                donor = &particles[ptcl->members[j]];
                                accretor = &particles[donor->get_binary_pair_id()];

                                assert(accretor->get_binary_interrupt_state() == BinaryInterruptState::collision);
                                assert(accretor->get_binary_pair_id() == donor->particle_index);
                                break;
                            }
                        }

                        Merge(donor, accretor);
#ifdef SEVN
                        if (donor->stellar_evolution != nullptr && accretor->stellar_evolution != nullptr) {
                            fprintf(stdout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->pid, donor->mass*mass_unit, donor->stellar_evolution->get_zams());
                            fprintf(stdout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->pid, accretor->mass*mass_unit, accretor->stellar_evolution->get_zams());
                            fprintf(stdout, "Donor: amiempty(): %d\n", donor->stellar_evolution->amiempty());
                            fprintf(stdout, "Accretor: amiempty(): %d\n", accretor->stellar_evolution->amiempty());
                            fflush(stdout);
                        }
#endif
                        int rank = cm_particle_worker_map[ptcl->particle_index];
                        queue.task = TASK_MERGE_MANYBODY;
                        queue.pid = ptcl->particle_index;
                        workers[rank].add_queue(queue);
                        workers[rank].run_queue();
                        workers[rank].callback();
#ifdef SEVN // This code is updated first in Enzo-Abyss by EW 2025.5.23
                        Particle* ptcl_erased = donor->mass < 0.0 ? donor : accretor;
                        fprintf(stdout, "ptcl_erased... PID: %d\n", ptcl_erased->pid);
                        if (ptcl_erased->stellar_evolution != nullptr) {

                            auto it = SEVNList.begin();
                            while (it != SEVNList.end()) {
                                if (it->second == ptcl_erased->particle_index) {
                                    it = SEVNList.erase(it);
                                    fprintf(stdout, "Merger induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl_erased->pid);
                                    break;
                                }
                                else
                                    it++;
                            }

                            delete ptcl_erased->stellar_evolution;
                            ptcl_erased->stellar_evolution = nullptr;
                            fprintf(stdout, "Merger induced zero mass particle (PID: %d) SEVN memory is free now\n", ptcl_erased->pid);
                        }
                        fflush(stdout);
#endif
                        continue;
                    }
                }

                bin_termination = true;

                if (ptcl->particle_index == last_particle_index) {
                    last_particle_index--;
                    g_state->last_particle_index = last_particle_index;
                }
                else
                    prev_cm_particle_worker_map.insert({ptcl->particle_index, cm_particle_worker_map[ptcl->particle_index]});
                cm_particle_worker_map.erase(ptcl->particle_index);

                for (int j=0; j < ptcl->num_members; j++) {
                    if (particles[ptcl->members[j]].mass < 0.0) {
#ifdef SEVN
                        Particle* ptcl_erased = &particles[ptcl->members[j]];
                        fprintf(stdout, "ptcl_erased... PID: %d\n", ptcl_erased->pid);
                        if (ptcl_erased->stellar_evolution != nullptr) {

                            auto it = SEVNList.begin();
                            while (it != SEVNList.end()) {
                                if (it->second == ptcl_erased->particle_index) {
                                    it = SEVNList.erase(it);
                                    fprintf(stdout, "Merger induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl_erased->pid);
                                    break;
                                }
                                else
                                    it++;
                            }

                            delete ptcl_erased->stellar_evolution;
                            ptcl_erased->stellar_evolution = nullptr;
                            fprintf(stdout, "Merger induced zero mass particle (PID: %d) SEVN memory is free now\n", ptcl_erased->pid);
                        }
                        fflush(stdout);
#endif
                        continue;
                    }
                    ThisLevelNode->ParticleList.push_back(ptcl->members[j]);
                }
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_map = std::chrono::high_resolution_clock::now();
#endif
                auto range = RegularMap.equal_range(ptcl->current_block_reg + ptcl->time_block_reg);
                for (auto it = range.first; it != range.second; ++it) {
                    if (ptcl->pid == particles[it->second].pid) {
                        RegularMap.erase(it);
                        break;
                    }
                }
#ifdef PERFORMANCETRACE
                end_point_map = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#endif // multimap

#ifdef SEVN_BINARY
                if (ptcl->binary_evolution != nullptr) {
#ifdef PERFORMANCETRACE
                    start_point_BSE = std::chrono::high_resolution_clock::now();
#endif
                    deleteSEVNBinary(ptcl);
#ifdef PERFORMANCETRACE
                    end_point_BSE = std::chrono::high_resolution_clock::now();
                    performance.BinaryStellarEvolution +=
                        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_BSE - start_point_BSE).count();
#endif
                }
#endif
                // Temporary test... If this works, RegularMap version should be updated too by EW 2025.08.18 // It works well!!
                RegularList.erase(ptcl->particle_index);
                FBTermination(ptcl);
            }
        }

        if (bin_termination) {
            for (int i=OriginalParticleListSize; i<ThisLevelNode->ParticleList.size(); i++) {
                ptcl = &particles[ThisLevelNode->ParticleList[i]];
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_map = std::chrono::high_resolution_clock::now();
#endif

                RegularMap.insert({ptcl->current_block_reg + ptcl->time_block_reg, ptcl->particle_index});

#ifdef PERFORMANCETRACE
                end_point_map = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#else // no multimap
                if (ptcl->current_block_reg + ptcl->time_block_reg == next_reg_time_block)
                    RegularList.insert(ptcl->particle_index);
#endif // multimap
                /*
                ptcl->new_num_members = 0;
                for (int j=OriginalParticleListSize; j<ThisLevelNode->ParticleList.size(); j++) {
                    if (i == j) continue;
                    ptcl->new_members[ptcl->new_num_members++] = ThisLevelNode->ParticleList[j];
                }
                */
                /*
                if (ptcl->time_step_irr < t_search)
                    ptcl->check_new_group_v4();
                */
            }

            // Erase terminated CM particles by EW 2025.1.6
            ThisLevelNode->ParticleList.erase(
                std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
                    [](int i) {
                    return !particles[i].is_active;
                    }
                ),
                ThisLevelNode->ParticleList.end()
            );
        }

#ifdef NSIGHT
        nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
        end_point_routine = std::chrono::high_resolution_clock::now();
        performance.FewBodyTermination +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
        PROFILE_STOP(TimerID::FewBodyTermination);
#ifdef unused
#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
        nvtxRangePushA("FewBodySearch");
#endif

#ifdef DEBUG
        std::cout << "FB search starts" << std::endl;
#endif
        /* // FB search is united with TASK_IRR_FORCE, so we don't need to do this again.
        // std::cerr << "FB search starts" << std::endl;
        // Few-body group search
        queue_scheduler.initialize(TASK_SEARCH_GROUP);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            // queue_scheduler.printStatus();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());

        // std::cerr << "FB search ended" << std::endl;
        */
        for (int ptcl_id : ThisLevelNode->ParticleList)
        {
            ptcl = &particles[ptcl_id];
            
            if (ptcl->new_num_members != 0)
                ptcl->check_new_group_v3();
            /*
            if (ptcl->get_binary_interrupt_state() == BinaryInterruptState::threebody) {
                ptcl->set_binary_interrupt_state(BinaryInterruptState::none);
            }
            else if (ptcl->get_binary_interrupt_state()==BinaryInterruptState::manybody)
            {
                ptcl->new_num_members = 0;
                ptcl->check_new_group_v2();
                ptcl->set_binary_interrupt_state(BinaryInterruptState::none);
            }
            else if (ptcl->new_num_members != 0)
            {	
                ptcl->check_new_group_v3();
            }
            */
        }
#ifdef DEBUG
        std::cout << "FB search ended" << std::endl;
#endif

#ifdef NSIGHT
        nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
        end_point_routine = std::chrono::high_resolution_clock::now();
        performance.FewBodySearch +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
#endif // unused
        PROFILE_START(TimerID::FewBodyInitialization);
#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
        nvtxRangePushA("FormBinaries");
#endif
        OriginalParticleListSize = ThisLevelNode->ParticleList.size();
        int rank_delete, rank_new;
#ifdef DEBUG
        std::cout << "formBinaries starts" << std::endl;
#endif
        formBinaries(ThisLevelNode->ParticleList, newCMptcls, cm_particle_worker_map, prev_cm_particle_worker_map);
#ifdef DEBUG
        std::cout << "formBinaries ended" << std::endl;
#endif
        if (OriginalParticleListSize != ThisLevelNode->ParticleList.size()) {
#ifdef DEBUG
            std::cout << "New Binary!" << std::endl;
#endif
            new_binaries = true;

            // new code by EW 2025.1.26
            Particle* ptclCM;
            Particle* mem_ptclCM;
            for (int i=0; i<newCMptcls.size(); i++) {
                ptclCM = &particles[newCMptcls[i]]; // 2025.01.10 edited to newCMptcls[i] by YS

                for (int j=0; j<ptclCM->new_num_members; j++) {
                    mem_ptclCM = &particles[ptclCM->new_members[j]];
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                    start_point_map = std::chrono::high_resolution_clock::now();
#endif

                    auto range = RegularMap.equal_range(mem_ptclCM->current_block_reg + mem_ptclCM->time_block_reg);
                    for (auto it = range.first; it != range.second; ++it) {
                        if (mem_ptclCM->pid == particles[it->second].pid) {
                            RegularMap.erase(it);
                            break;
                        }
                    }
#ifdef PERFORMANCETRACE
                    end_point_map = std::chrono::high_resolution_clock::now();
                    performance.RegularMap +=
                        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#endif // multimap
                    if (mem_ptclCM->is_cm_particle) {
                        fprintf(stdout, "manybody group detected; PID %d should be deleted first\n", mem_ptclCM->pid);
                        rank_delete = cm_particle_worker_map[mem_ptclCM->particle_index];
                        fprintf(stdout, "Rank of CM ptcl %d: %d\n", mem_ptclCM->pid, rank_delete);
                        queue.task = TASK_DELETE_GROUP;
                        queue.pid = mem_ptclCM->particle_index;
                        workers[rank_delete].add_queue(queue);
                        workers[rank_delete].run_queue();
                        workers[rank_delete].callback();

                        if (mem_ptclCM->particle_index == last_particle_index) {
                            last_particle_index--;
                            g_state->last_particle_index = last_particle_index;
                        }
                        else
                            prev_cm_particle_worker_map.insert({mem_ptclCM->particle_index, cm_particle_worker_map[mem_ptclCM->particle_index]});
                        cm_particle_worker_map.erase(mem_ptclCM->particle_index);
                    }
                }

                rank_new = cm_particle_worker_map[ptclCM->particle_index];
#ifdef DEBUG
                fprintf(stdout, "Rank of CM ptcl %d: %d\n", ptclCM->pid, rank_new);
#endif

#ifdef SEVN_BINARY
#ifdef PERFORMANCETRACE
                start_point_BSE = std::chrono::high_resolution_clock::now();
#endif
                if (!makeSEVNBinary(ptclCM)) {
                    if (ptclCM->particle_index == last_particle_index) {
                        last_particle_index--;
                        g_state->last_particle_index = last_particle_index;
                    }
                    else
                        prev_cm_particle_worker_map.insert({ptclCM->particle_index, cm_particle_worker_map[ptclCM->particle_index]});
                    cm_particle_worker_map.erase(ptclCM->particle_index);
                    continue;
                }
#ifdef PERFORMANCETRACE
                end_point_BSE = std::chrono::high_resolution_clock::now();
                performance.BinaryStellarEvolution +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_BSE - start_point_BSE).count();
#endif
#endif
                queue.task = TASK_MAKE_GROUP;
                queue.pid = ptclCM->particle_index;
                workers[rank_new].add_queue(queue);
                workers[rank_new].run_queue();
                workers[rank_new].callback();
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_map = std::chrono::high_resolution_clock::now();
#endif
                RegularMap.insert({ptclCM->current_block_reg + ptclCM->time_block_reg, ptclCM->particle_index});
#ifdef PERFORMANCETRACE
                end_point_map = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#else // no multimap
                // Temporary test... If this works, RegularMap version should be updated too by EW 2025.08.18 // It works well!!
                for (int i = 0; i < ptclCM->new_num_members; i++) {
                    RegularList.erase(ptclCM->new_members[i]);
                    particles[ptclCM->new_members[i]].new_num_members = 0;
                }
                if (ptclCM->current_block_reg + ptclCM->time_block_reg == next_reg_time_block)
                    RegularList.insert(ptclCM->particle_index); // VERY IMPORTANT BUG FIXED by EW 2025.7.18
#endif
            }
#ifdef DEBUG
            std::cout << "All new fewbody objects are initialized." << std::endl;
#endif

            ThisLevelNode->ParticleList.erase(
                std::remove_if(
                        ThisLevelNode->ParticleList.begin(), 
                        ThisLevelNode->ParticleList.end(),
                        [](int i) { return !particles[i].is_active; }
                ),
                ThisLevelNode->ParticleList.end()
            );
        }
        newCMptcls.clear();

#ifdef NSIGHT
        nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
        end_point_routine = std::chrono::high_resolution_clock::now();
        performance.FewBodyInitialization +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
        PROFILE_STOP(TimerID::FewBodyInitialization);

#endif

        PROFILE_START(TimerID::SkipListUpdate);
#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif
#ifdef DEBUG
        std::cout << "updateSkipList starts" << std::endl;
#endif
        OriginalParticleListSize = ThisLevelNode->ParticleList.size();
        for (int i=0; i<OriginalParticleListSize; i++)
            updateSkipList(skiplist, ThisLevelNode->ParticleList[i]);
#ifdef DEBUG
        std::cout << "updateSkipList ended" << std::endl;
#endif
#ifdef PERFORMANCETRACE
        end_point_routine = std::chrono::high_resolution_clock::now();
        performance.SkipListUpdate +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif
        PROFILE_STOP(TimerID::SkipListUpdate);

        current_time_irr = particles[ThisLevelNode->ParticleList[0]].current_block_irr*time_step;
#ifdef DEBUG
        std::cout << "skiplist->deleteFirstNode() starts" << std::endl;
#endif
        skiplist->deleteFirstNode();
#ifdef DEBUG
        std::cout << "skiplist->deleteFirstNode() ended" << std::endl;
#endif

#ifdef IRR_TEST
        std::cout << "current_time_irr=" << current_time_irr<< std::endl;
        // create output at appropriate time intervals
        if (current_time_irr >= output_time) {
            writeParticle(current_time_irr, output_num++);
            output_time += output_time_step;
        }

        // end if the global time exceeds the end time
        if (current_time_irr >= 1) {
            task=TASK_END;
            queue = {task, -1, -1.0};
            InitialAssignmentOfTasks(Queue, num_workers, QUEUE_TAG);
            fprintf(stdout, "Simulation Done! Current Time: %e Myr\n", global_time*enzo_time_step*1e4);
            return;
        }
#endif
    } // Irr

#ifdef DEBUG
    std::cout << "delete skiplist" << std::endl;
#endif
    delete skiplist;
    skiplist = nullptr;
    //exit(SUCCESS);

#ifdef FEWBODY
    if (bin_termination || new_binaries) {
#ifdef MULTIMAP
        if (next_reg_time_block != RegularMap.begin()->first) {
            next_reg_time_block = RegularMap.begin()->first;
            g_state->next_reg_time_block = next_reg_time_block;
            continue;
        }
#else // no multimap
        /* // I didn't erase this yet because MultiMap should be fixed too by EW 2025.8.18
        for (auto it = RegularList.begin(); it != RegularList.end(); ) {
            if (!particles[*it].is_active) {
                it = RegularList.erase(it);
                fprintf(stderr, "bin_term or new_binaries... Why PID: %d is inactive? (NRTB: %llu)\n", particles[*it].pid, next_reg_time_block);
                assert(particles[*it].is_active);
            }
            else
                ++it;
        }
        */
        if (RegularList.empty())
            return false;
#endif
    }
#endif
    return true;
}

bool createSkipList(SkipList *skiplist) {

	bool debug = false;
	//fprintf(stdout, "create level starts!\n");
	//fflush(stdout);

	if (debug) {
		fprintf(stdout, "create level starts!\n");
		fflush(stdout);
	}

	Particle* ptcl;

	for (int i=0; i<=last_particle_index; i++) {
		ptcl =  &particles[i];

		// if ((ptcl->num_neighbors != 0) && (ptcl->next_block_irr <= next_reg_time_block)) { // IAR original
		if (ptcl->is_active && ptcl->next_block_irr <= next_reg_time_block) {	// IAR modified
			//fprintf(stdout, "PID=%d, NBI=%llu\n", ptcl->pid, ptcl->next_block_irr);
			if (!skiplist->search(ptcl->next_block_irr, ptcl->particle_index))
				skiplist->insert(ptcl->next_block_irr, ptcl->particle_index);
		}
	}

	/*
	if (debug) {
		skiplist->display();
	}
	*/

	//fprintf(stdout, "create level ends!\n");
	//fflush(stdout);

	if (skiplist->getFirstNode() == nullptr)
		return false;
	else
		return true;
}



bool updateSkipList(SkipList *skiplist, int ptcl_id) {
	bool debug = false;
	if (debug) {
		fprintf(stdout, "update level starts!\n");
		fflush(stdout);
	}

	/* Update New Time Steps */
	//Node* ThisLevelNode = skiplist->getFirstNode();

	/*
	if (this->debug) {
	fprintf(stdout, "PID=%d, NBI=%llu, size=%lu\n", ptcl->pid, ptcl->next_block_irr, ThisLevelNode->particle_list.size());
	fprintf(stdout, "NextBlockIrr=%llu\n",ptcl->next_block_irr);
	fflush(stdout);
	}
	*/

	Particle * ptcl = &particles[ptcl_id];

	//std::cout << "NextBlockIrr of "<< ptcl_id<<" = " << ptcl->next_block_irr << std::endl;
	if (ptcl->next_block_irr > next_reg_time_block)
		return true;

	if (!skiplist->search(ptcl->next_block_irr, ptcl->particle_index))
		skiplist->insert(ptcl->next_block_irr, ptcl->particle_index);

	if (debug) {
	}


	if (debug) {
		//skiplist->display();
		//fprintf(stderr, "This is it.\n\n\n");
	}

	return true;
}