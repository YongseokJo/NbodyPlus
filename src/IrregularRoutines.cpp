#include "QueueScheduler.h"
#include "SkipList.h"
#include <algorithm>
#ifdef NSIGHT
#include <nvToolsExt.h>
#endif


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

    // Irregular
    while ( skiplist->getFirstNode() != nullptr) {

        ThisLevelNode = skiplist->getFirstNode();
        ThisLevelNode->ParticleList.erase(
            std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
                [](int i) {
                return !particles[i].isActive;
                }
            ),
            ThisLevelNode->ParticleList.end()
        );
        if (ThisLevelNode->ParticleList.size() == 0) {
            skiplist->deleteFirstNode();
            continue;
        }

        next_time     = particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr\
                                    + particles[ThisLevelNode->ParticleList[0]].TimeStepIrr;

#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif
    
#ifdef DEBUG
        // print out particlelist
        fprintf(stdout, "(IRR_FORCE) next_time: %e Myr\n", next_time*EnzoTimeStep*1e4);
        // /*
        fprintf(stdout, "PID: %d. CurrentTimeIrr: %e Myr, TimeStepIrr: %e Myr\n", 
                    particles[ThisLevelNode->ParticleList[0]].PID, 
                    particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr*EnzoTimeStep*1e4, 
                    particles[ThisLevelNode->ParticleList[0]].TimeStepIrr*EnzoTimeStep*1e4);

        // fprintf(stdout, "PID (%d) = ", ThisLevelNode->ParticleList.size());
        for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            // fprintf(stdout, "%d, ", ptcl->PID);
            fprintf(stdout, "PID: %d. %e Myr, %e Myr\n", 
                    ptcl->PID,
                    ptcl->CurrentTimeIrr*EnzoTimeStep*1e4,
                    ptcl->TimeStepIrr*EnzoTimeStep*1e4);
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
        queue_scheduler.initializeIrr(IrrForce, next_time, ThisLevelNode->ParticleList);
        auto iter = queue_scheduler.CMPtcls.begin();
#ifdef SEVN_BINARY
        std::vector<int> CMPtclsForSEVN;
#endif
        do {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            do {
                worker = queue_scheduler.waitQueue(1); // non-blocking wait
                // if there's any CMPtcl
                if (queue_scheduler.CMPtcls.size() > 0)
                {
                    // check if there's any CM ptcl ready to go for SDAR
                    if (iter == queue_scheduler.CMPtcls.end())
                        iter = queue_scheduler.CMPtcls.begin();
                    cm_pid = *(iter);
                    ptcl = &particles[cm_pid];
                    /*
                    for (int j = 0; j < ptcl->NumberOfNeighbor; j++)
                    {
                        // if (particles[ptcl->Neighbors[j]].isUpdateToDate == false) // original code
                        if (particles[ptcl->Neighbors[j]].isActive && !particles[ptcl->Neighbors[j]].isUpdateToDate) // modified by EW 2025.2.26
                        {
                            iter++;
                            goto skip_to_next;
                        }
                    }
                    */
                    queue.task = ARIntegration;
                    queue.pid = cm_pid;
                    queue.next_time = next_time;
                    workers[CMPtclWorker[cm_pid]].addQueue(queue);
                    queue_scheduler.assignWorker(&workers[CMPtclWorker[cm_pid]]);
                    iter = queue_scheduler.CMPtcls.erase(iter);
#ifdef SEVN_BINARY
                    if (ptcl->BinaryEvolution != nullptr)
                        CMPtclsForSEVN.push_back(cm_pid);
#endif
                // skip_to_next:;
                }
                if (worker != nullptr) 
                {

                }
            } while (worker == nullptr);
            queue_scheduler.callback(worker);
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
        queue_scheduler.initialize(IrrForce, next_time);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());

        for (int ptcl_id : ThisLevelNode->ParticleList)
        {
            ptcl = &particles[ptcl_id];

            if (ptcl->isCMptcl) {
                int rank = CMPtclWorker[ptcl->ParticleIndex];
                queue.task = ARIntegration;
                queue.pid = ptcl->ParticleIndex;
                queue.next_time = next_time;
                workers[rank].addQueue(queue);
                workers[rank].runQueue();
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
        queue_scheduler.initialize(IrrForce, next_time);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());
#endif

#ifdef PERFORMANCETRACE
        end_point_routine = std::chrono::high_resolution_clock::now();
        performance.IrregularForce +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif

        //ParticleSynchronization();


#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
        nvtxRangePushA("IrregularUpdate");
#endif
        /* // For this kind of simple work, using queue_scheduler is slower
        queue_scheduler.initialize(IrrUpdate);
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

            if (ptcl->NumberOfNeighbor != 0) // IAR modified
                ptcl->updateParticle();
            ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
            ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*time_step;
        }
#ifdef DEBUG
        for (int i: ThisLevelNode->ParticleList) {
            ptcl = &particles[i];
            if (ptcl->CurrentTimeIrr != next_time) {
                fprintf(stdout, "Error! PID: %d, CurrentTimeIrr: %e Myr, next_time: %e Myr\n", ptcl->PID, ptcl->CurrentTimeIrr*EnzoTimeStep*1e4, next_time*EnzoTimeStep*1e4);
                assert(ptcl->CurrentTimeIrr == next_time);
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
        int OriginalParticleListSize;
#ifdef FEWBODY

#ifdef PERFORMANCETRACE
        start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
        nvtxRangePushA("FewBodyTermination");
#endif
        OriginalParticleListSize = ThisLevelNode->ParticleList.size();
        for (int i=0; i<OriginalParticleListSize; i++) {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            if (ptcl->getBinaryInterruptState() == BinaryInterruptState::merger ||
                ptcl->getBinaryInterruptState() == BinaryInterruptState::terminated) {

                assert(ptcl->isCMptcl); // for debugging by EW 2025.1.20
            
                if (ptcl->getBinaryInterruptState() == BinaryInterruptState::merger) {
                    if (ptcl->NumberOfMember == 2) { // binary merger

                        Particle* donor = &particles[ptcl->Members[0]];
                        Particle* accretor = &particles[ptcl->Members[1]];

                        Merge(donor, accretor);
#ifdef SEVN
                        if (donor->StellarEvolution != nullptr && accretor->StellarEvolution != nullptr) {
                            fprintf(stdout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->PID, donor->Mass*mass_unit, donor->StellarEvolution->get_zams());
                            fprintf(stdout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->PID, accretor->Mass*mass_unit, accretor->StellarEvolution->get_zams());
                            fprintf(stdout, "Donor: amiempty(): %d\n", donor->StellarEvolution->amiempty());
                            fprintf(stdout, "Accretor: amiempty(): %d\n", accretor->StellarEvolution->amiempty());
                            fflush(stdout);
                        }
#endif
                    }
                    else { // from NewFBInitialization3

                        assert(ptcl->NumberOfMember > 2); // for debugging by EW 2025.1.20

                        Particle* donor;
                        Particle* accretor;

                        for (int j=0; j<ptcl->NumberOfMember; j++) {
                            if (particles[ptcl->Members[j]].getBinaryInterruptState() == BinaryInterruptState::collision) {
                                donor = &particles[ptcl->Members[j]];
                                accretor = &particles[donor->getBinaryPairID()];

                                assert(accretor->getBinaryInterruptState() == BinaryInterruptState::collision);
                                assert(accretor->getBinaryPairID() == donor->ParticleIndex);
                                break;
                            }
                        }

                        Merge(donor, accretor);
#ifdef SEVN
                        if (donor->StellarEvolution != nullptr && accretor->StellarEvolution != nullptr) {
                            fprintf(stdout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->PID, donor->Mass*mass_unit, donor->StellarEvolution->get_zams());
                            fprintf(stdout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->PID, accretor->Mass*mass_unit, accretor->StellarEvolution->get_zams());
                            fprintf(stdout, "Donor: amiempty(): %d\n", donor->StellarEvolution->amiempty());
                            fprintf(stdout, "Accretor: amiempty(): %d\n", accretor->StellarEvolution->amiempty());
                            fflush(stdout);
                        }
#endif
                        int rank = CMPtclWorker[ptcl->ParticleIndex];
                        queue.task = MergeManyBody;
                        queue.pid = ptcl->ParticleIndex;
                        workers[rank].addQueue(queue);
                        workers[rank].runQueue();
                        workers[rank].callback();
#ifdef SEVN // This code is updated first in Enzo-Abyss by EW 2025.5.23
                        Particle* ptcl_erased = donor->Mass < 0.0 ? donor : accretor;
                        fprintf(stdout, "ptcl_erased... PID: %d\n", ptcl_erased->PID);
                        if (ptcl_erased->StellarEvolution != nullptr) {

                            auto it = SEVNList.begin();
                            while (it != SEVNList.end()) {
                                if (it->second == ptcl_erased->ParticleIndex) {
                                    it = SEVNList.erase(it);
                                    fprintf(stdout, "Merger induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl_erased->PID);
                                    break;
                                }
                                else
                                    it++;
                            }

                            delete ptcl_erased->StellarEvolution;
                            ptcl_erased->StellarEvolution = nullptr;
                            fprintf(stdout, "Merger induced zero mass particle (PID: %d) SEVN memory is free now\n", ptcl_erased->PID);
                        }
                        fflush(stdout);
#endif
                        continue;
                    }
                }

                bin_termination = true;

                if (ptcl->ParticleIndex == LastParticleIndex) {
                    LastParticleIndex--;
                    global_variable->LastParticleIndex = LastParticleIndex;
                }
                else
                    PrevCMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker[ptcl->ParticleIndex]});
                CMPtclWorker.erase(ptcl->ParticleIndex);

                for (int j=0; j < ptcl->NumberOfMember; j++) {
                    if (particles[ptcl->Members[j]].Mass < 0.0) {
#ifdef SEVN
                        Particle* ptcl_erased = &particles[ptcl->Members[j]];
                        fprintf(stdout, "ptcl_erased... PID: %d\n", ptcl_erased->PID);
                        if (ptcl_erased->StellarEvolution != nullptr) {

                            auto it = SEVNList.begin();
                            while (it != SEVNList.end()) {
                                if (it->second == ptcl_erased->ParticleIndex) {
                                    it = SEVNList.erase(it);
                                    fprintf(stdout, "Merger induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl_erased->PID);
                                    break;
                                }
                                else
                                    it++;
                            }

                            delete ptcl_erased->StellarEvolution;
                            ptcl_erased->StellarEvolution = nullptr;
                            fprintf(stdout, "Merger induced zero mass particle (PID: %d) SEVN memory is free now\n", ptcl_erased->PID);
                        }
                        fflush(stdout);
#endif
                        continue;
                    }
                    ThisLevelNode->ParticleList.push_back(ptcl->Members[j]);
                }
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_map = std::chrono::high_resolution_clock::now();
#endif
                auto range = RegularMap.equal_range(ptcl->CurrentBlockReg + ptcl->TimeBlockReg);
                for (auto it = range.first; it != range.second; ++it) {
                    if (ptcl->PID == particles[it->second].PID) {
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
                if (ptcl->BinaryEvolution != nullptr) {
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
                RegularList.erase(ptcl->ParticleIndex);
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

                RegularMap.insert({ptcl->CurrentBlockReg + ptcl->TimeBlockReg, ptcl->ParticleIndex});

#ifdef PERFORMANCETRACE
                end_point_map = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#else // no multimap
                if (ptcl->CurrentBlockReg + ptcl->TimeBlockReg == NextRegTimeBlock)
                    RegularList.insert(ptcl->ParticleIndex);
#endif // multimap
                /*
                ptcl->NewNumberOfMember = 0;
                for (int j=OriginalParticleListSize; j<ThisLevelNode->ParticleList.size(); j++) {
                    if (i == j) continue;
                    ptcl->NewMembers[ptcl->NewNumberOfMember++] = ThisLevelNode->ParticleList[j];
                }
                */
                /*
                if (ptcl->TimeStepIrr * EnzoTimeStep * 1e4 < TSEARCH)
                    ptcl->checkNewGroup4();
                */
            }

            // Erase terminated CM particles by EW 2025.1.6
            ThisLevelNode->ParticleList.erase(
                std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
                    [](int i) {
                    return !particles[i].isActive;
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
        /* // FB search is united with IrrForce, so we don't need to do this again.
        // std::cerr << "FB search starts" << std::endl;
        // Few-body group search
        queue_scheduler.initialize(SearchGroup);
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
            
            if (ptcl->NewNumberOfMember != 0)
                ptcl->checkNewGroup3();
            /*
            if (ptcl->getBinaryInterruptState() == BinaryInterruptState::threebody) {
                ptcl->setBinaryInterruptState(BinaryInterruptState::none);
            }
            else if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody)
            {
                ptcl->NewNumberOfMember = 0;
                ptcl->checkNewGroup2();
                ptcl->setBinaryInterruptState(BinaryInterruptState::none);
            }
            else if (ptcl->NewNumberOfMember != 0)
            {	
                ptcl->checkNewGroup3();
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
        formBinaries(ThisLevelNode->ParticleList, newCMptcls, CMPtclWorker, PrevCMPtclWorker);
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

                for (int j=0; j<ptclCM->NewNumberOfMember; j++) {
                    mem_ptclCM = &particles[ptclCM->NewMembers[j]];
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                    start_point_map = std::chrono::high_resolution_clock::now();
#endif

                    auto range = RegularMap.equal_range(mem_ptclCM->CurrentBlockReg + mem_ptclCM->TimeBlockReg);
                    for (auto it = range.first; it != range.second; ++it) {
                        if (mem_ptclCM->PID == particles[it->second].PID) {
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
                    if (mem_ptclCM->isCMptcl) {
                        fprintf(stdout, "manybody group detected; PID %d should be deleted first\n", mem_ptclCM->PID);
                        rank_delete = CMPtclWorker[mem_ptclCM->ParticleIndex];
                        fprintf(stdout, "Rank of CM ptcl %d: %d\n", mem_ptclCM->PID, rank_delete);
                        queue.task = DeleteGroup;
                        queue.pid = mem_ptclCM->ParticleIndex;
                        workers[rank_delete].addQueue(queue);
                        workers[rank_delete].runQueue();
                        workers[rank_delete].callback();

                        if (mem_ptclCM->ParticleIndex == LastParticleIndex) {
                            LastParticleIndex--;
                            global_variable->LastParticleIndex = LastParticleIndex;
                        }
                        else
                            PrevCMPtclWorker.insert({mem_ptclCM->ParticleIndex, CMPtclWorker[mem_ptclCM->ParticleIndex]});
                        CMPtclWorker.erase(mem_ptclCM->ParticleIndex);
                    }
                }

                rank_new = CMPtclWorker[ptclCM->ParticleIndex];
#ifdef DEBUG
                fprintf(stdout, "Rank of CM ptcl %d: %d\n", ptclCM->PID, rank_new);
#endif

#ifdef SEVN_BINARY
#ifdef PERFORMANCETRACE
                start_point_BSE = std::chrono::high_resolution_clock::now();
#endif
                if (!makeSEVNBinary(ptclCM)) {
                    if (ptclCM->ParticleIndex == LastParticleIndex) {
                        LastParticleIndex--;
                        global_variable->LastParticleIndex = LastParticleIndex;
                    }
                    else
                        PrevCMPtclWorker.insert({ptclCM->ParticleIndex, CMPtclWorker[ptclCM->ParticleIndex]});
                    CMPtclWorker.erase(ptclCM->ParticleIndex);
                    continue;
                }
#ifdef PERFORMANCETRACE
                end_point_BSE = std::chrono::high_resolution_clock::now();
                performance.BinaryStellarEvolution +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_BSE - start_point_BSE).count();
#endif
#endif
                queue.task = MakeGroup;
                queue.pid = ptclCM->ParticleIndex;
                workers[rank_new].addQueue(queue);
                workers[rank_new].runQueue();
                workers[rank_new].callback();
#ifdef MULTIMAP
#ifdef PERFORMANCETRACE
                start_point_map = std::chrono::high_resolution_clock::now();
#endif
                RegularMap.insert({ptclCM->CurrentBlockReg + ptclCM->TimeBlockReg, ptclCM->ParticleIndex});
#ifdef PERFORMANCETRACE
                end_point_map = std::chrono::high_resolution_clock::now();
                performance.RegularMap +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_map - start_point_map).count();
#endif
#else // no multimap
                // Temporary test... If this works, RegularMap version should be updated too by EW 2025.08.18 // It works well!!
                for (int i = 0; i < ptclCM->NewNumberOfMember; i++) {
                    RegularList.erase(ptclCM->NewMembers[i]);
                    particles[ptclCM->NewMembers[i]].NewNumberOfMember = 0;
                }
                if (ptclCM->CurrentBlockReg + ptclCM->TimeBlockReg == NextRegTimeBlock)
                    RegularList.insert(ptclCM->ParticleIndex); // VERY IMPORTANT BUG FIXED by EW 2025.7.18
#endif
            }
#ifdef DEBUG
            std::cout << "All new fewbody objects are initialized." << std::endl;
#endif

            ThisLevelNode->ParticleList.erase(
                std::remove_if(
                        ThisLevelNode->ParticleList.begin(), 
                        ThisLevelNode->ParticleList.end(),
                        [](int i) { return !particles[i].isActive; }
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

#endif

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

        current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr*time_step;
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
        if (current_time_irr >= outputTime) {
            writeParticle(current_time_irr, outNum++);
            outputTime += outputTimeStep;
        }

        // end if the global time exceeds the end time
        if (current_time_irr >= 1) {
            task=Ends;
            queue = {task, -1, -1.0};
            InitialAssignmentOfTasks(Queue, NumberOfWorker, QUEUE_TAG);
            fprintf(stdout, "Simulation Done! Current Time: %e Myr\n", global_time*EnzoTimeStep*1e4);
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
        if (NextRegTimeBlock != RegularMap.begin()->first) {
            NextRegTimeBlock = RegularMap.begin()->first;
            global_variable->NextRegTimeBlock = NextRegTimeBlock;
            continue;
        }
#else // no multimap
        /* // I didn't erase this yet because MultiMap should be fixed too by EW 2025.8.18
        for (auto it = RegularList.begin(); it != RegularList.end(); ) {
            if (!particles[*it].isActive) {
                it = RegularList.erase(it);
                fprintf(stderr, "bin_term or new_binaries... Why PID: %d is inactive? (NRTB: %llu)\n", particles[*it].PID, NextRegTimeBlock);
                assert(particles[*it].isActive);
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

	for (int i=0; i<=LastParticleIndex; i++) {
		ptcl =  &particles[i];

		// if ((ptcl->NumberOfNeighbor != 0) && (ptcl->NextBlockIrr <= NextRegTimeBlock)) { // IAR original
		if (ptcl->isActive && ptcl->NextBlockIrr <= NextRegTimeBlock) {	// IAR modified
			//fprintf(stdout, "PID=%d, NBI=%llu\n", ptcl->PID, ptcl->NextBlockIrr);
			if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
				skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);
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
	fprintf(stdout, "PID=%d, NBI=%llu, size=%lu\n", ptcl->PID, ptcl->NextBlockIrr, ThisLevelNode->particle_list.size());
	fprintf(stdout, "NextBlockIrr=%llu\n",ptcl->NextBlockIrr);
	fflush(stdout);
	}
	*/

	Particle * ptcl = &particles[ptcl_id];

	//std::cout << "NextBlockIrr of "<< ptcl_id<<" = " << ptcl->NextBlockIrr << std::endl;
	if (ptcl->NextBlockIrr > NextRegTimeBlock)
		return true;

	if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
		skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);

	if (debug) {
	}


	if (debug) {
		//skiplist->display();
		//fprintf(stderr, "This is it.\n\n\n");
	}

	return true;
}