#include <iostream>
#include <vector>
#include <errno.h>
#include <unordered_set>
#include "global.h"
#include "Queue.h"
#include "performance.h"

void ComputeAcceleration(int ptcl_id, double next_time);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void makePrimordialGroup(Particle* ptclCM);
void NewFBInitialization(Particle* ptclCM);
void deleteGroup(Particle* ptclCM);
void FBdeleteGroup(Group* group);
void NewFBInitialization3(Group* group);


#define noWORKER_DEBUG

void WorkerRoutines() {

	//std::cout << "Processor " << MyRank << " is ready." << std::endl;

	//TaskName task = Error;
	int task;
	MPI_Status status;
	MPI_Request request;
	int ptcl_id;
	double next_time;
	int NewNumberOfNeighbor;
	int NewNeighbors[NumNeighborMax];
	int size=0;
	double new_a[Dim];
	double new_adot[Dim];
	Particle *ptcl;
	Group *group;
	std::unordered_set<int> MyCMPtcls;
	//int *MyQueue = &queues[MyRank];
	int MyCurrentQueue;

    MPI_Request requests[MAX_COMMUNICATION];
    MPI_Status statuses[MAX_COMMUNICATION];
    int NumberOfCommunication=0;

	while (true) {
		MPI_Recv(NULL, 0, MPI_BYTE, ROOT, TASK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//MPI_Recv(&task, 1, MPI_INT, ROOT, TASK_TAG, MPI_COMM_WORLD, &statuses[0]);
		//MPI_Win_sync(win);
		//MPI_Win_sync(win4);
		//MPI_Win_sync(win5);
		//task = tasks[MyRank];
		task = tasks[0];
		MyCurrentQueue = MyRank-1;
		//std::cerr << "Processor " << MyRank << " received task " << task << std::endl;
		#ifdef WORKER_DEBUG
		fprintf(workerout, "Received task %d\n", task);
		#endif

		switch (task) {

/*=======================================*
 *           Irregular Force             *
 *=======================================*/
			case IrrForce: // Irregular Acceleration
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				MPI_Win_sync(win4);
				MPI_Win_sync(win5);
#ifdef IRR_CM_SEPARATE
				next_time = global_variable->next_time;

#ifdef WORKER_DEBUG
				fprintf(workerout, "IrrForce starts, next_time =  %e\n", global_variable->next_time);
				fprintf(workerout, "Queues = ");
				for (int i = 0; i < global_variable->QueueSize; i++)
				{
					fprintf(workerout, "%d, ", queues[i]);
				}
				fprintf(workerout, "\n");
#endif

				while (MyCurrentQueue < global_variable->QueueSize)
				{
#ifdef WORKER_DEBUG
					fprintf(workerout, "MyCurrentQueue =  %d\n", MyCurrentQueue);
#endif
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED)
					{
						std::cerr << "Queue skipped" << std::endl;
						exit(1);
#ifdef WORKER_DEBUG
					fprintf(workerout, "skipped\n");
#endif
					MyCurrentQueue++;
					continue;
				}
				queues[MyCurrentQueue] = ALREADY_COMPUTED;
				// std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
				ptcl = &particles[ptcl_id];
				ptcl->computeAccelerationIrr();
				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->calculateTimeStepIrr();
				ptcl->NextBlockIrr = ptcl->NewCurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->isUpdateToDate = true;
				// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
				// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
				MPI_Send(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD);
				// MyCurrentQueue++;
				MyCurrentQueue += NumberOfWorker;
#ifdef WORKER_DEBUG
				 fprintf(workerout, "PID= %d,CurrentTimeIrr =  %e\n", ptcl_id, ptcl->CurrentTimeIrr);
#endif
			}
#else
#endif
#ifdef WORKER_DEBUG
				fprintf(workerout, "IrrForce ends\n--------------------------\n");
				fflush(workerout);
				#endif
				break;


/*=======================================*
 *           Regular Update              *
 *=======================================*/
			case RegForce: // Regular Acceleration
		MPI_Win_fence(0, win4);
		MPI_Win_fence(0, win5);
				//std::cout << "RegCal start " << MyRank << std::endl;
				MPI_Recv(&ptcl_id,   1, MPI_INT,    ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, MPI_COMM_WORLD, &status);

				particles[ptcl_id].computeAccelerationReg();
				//ComputeAcceleration(ptcl_id, next_time);
				//std::cout << "RegCal end" << MyRank << std::endl;
				break;


/*=======================================*
 *           Irregular Update            *
 *=======================================*/
			case IrrUpdate: // Irregular Update Particle
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
        MPI_Win_sync(win4);
        MPI_Win_sync(win5);
				// std::cout << "IrrUp Processor " << MyRank << " " << global_variable->next_time <<std::endl;
#ifdef WORKER_DEBUG
				fprintf(workerout, "IrrUpdate starts, next_time =  %e\n", global_variable->next_time);
				fprintf(workerout, "Queues = ");
				for (int i = 0; i < global_variable->QueueSize; i++)
				{
					fprintf(workerout, "%d, ", queues[i]);
				}
				fprintf(workerout, "\n");
#endif
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					//MPI_Win_fence(0, win5); 
					ptcl_id = queues[MyCurrentQueue];
					//std::cout << "Processor " << MyRank << ": MyCurrentQueue= " << MyCurrentQueue <<
					// "  QueueSize= " << global_variable->QueueSize << std::endl;
					#ifdef WORKER_DEBUG
					fprintf(workerout, "MyCurrentQueue =  (%d/%d), PID = %d\n",
					 MyCurrentQueue, global_variable->QueueSize, ptcl_id);
					#endif
					if (ptcl_id == ALREADY_COMPUTED)
					{
						std::cerr << "Queue skipped" << std::endl;
						exit(1);
#ifdef WORKER_DEBUG
					fprintf(workerout, "skipped\n");
#endif
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					ptcl = &particles[ptcl_id];
					if (ptcl->NumberOfNeighbor != 0) // IAR modified
						ptcl->updateParticle();
					ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockIrr * time_step;
					// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					//MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					MPI_Send(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD);
					//MyCurrentQueue++;
					MyCurrentQueue += NumberOfWorker;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					//std::cout << "Processor " << MyRank << ", pid=" << ptcl_id << ", CurrentTimeIrr=" << ptcl->CurrentTimeIrr << std::endl;
					#ifdef WORKER_DEBUG
					fprintf(workerout, "PID= %d,CurrentTimeIrr =  %e\n", ptcl_id, ptcl->CurrentTimeIrr);
					#endif
				}
				#ifdef WORKER_DEBUG
				fprintf(workerout, "IrrUpdate ends\n=================================\n");
				fflush(workerout);
				#endif
				//std::cout  << "Processor " << MyRank << " IrrUp end!" << std::endl;
				break;


/*=======================================*
 *           Regular Update              *
 *=======================================*/
			case RegUpdate: // Regular Update Particle
		MPI_Win_fence(0, win4);
		MPI_Win_fence(0, win5);
				//std::cout << "RegUp start " << MyRank << std::endl;
				MPI_Recv(&ptcl_id, 1, MPI_INT, ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				//std::cout << "ptcl " << ptcl_id << std::endl;

				ptcl = &particles[ptcl_id];
				ptcl->updateParticle();

				for (int i=0; i<ptcl->NewNumberOfNeighbor; i++)
					ptcl->Neighbors[i] = ptcl->NewNeighbors[i];
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->CurrentBlockReg += ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg   = ptcl->CurrentBlockReg*time_step;
				ptcl->calculateTimeStepReg();
				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockReg;
				ptcl->calculateTimeStepIrr();
				ptcl->updateRadius();
				if (ptcl->NumberOfNeighbor == 0) {
					ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*time_step;
				}
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				break;


/*=======================================*
 *           Regular CUDA                *
 *=======================================*/
			case RegCuda: // Update Regular Particle CUDA
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED)
					{
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					ptcl->updateRegularParticleCuda();
					// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					//MyCurrentQueue++;
					MyCurrentQueue += NumberOfWorker;
				}
				break;


/*=======================================*
 *           Regular CUDA Update         *
 *=======================================*/
			case RegCudaUpdate: // Update Regular Particle CUDA II
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED)
					{
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					for (int j = 0; j < ptcl->NewNumberOfNeighbor; j++)
						ptcl->Neighbors[j] = ptcl->NewNeighbors[j];
					ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

					ptcl->updateParticle();
					ptcl->CurrentBlockReg = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;
					ptcl->CurrentTimeReg = ptcl->CurrentBlockReg * time_step;
					ptcl->calculateTimeStepReg();
					ptcl->calculateTimeStepIrr();
					if (ptcl->NumberOfNeighbor == 0)
					{
						ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
						ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg * time_step;
					}
					ptcl->updateRadius();
					ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of ptcl particle
					// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					//MyCurrentQueue++;
					MyCurrentQueue += NumberOfWorker;
				}
				break;


/*=======================================*
 *           Initialize Acceleration 1   *
 *=======================================*/
			case InitAcc1: // Initialize Acceleration(01)
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				//std::cout << "Processor " << MyRank << " initialization starts." << std::endl;
				//fprintf(workerout, "InitAcc1 starts, next_time =  %e\n", global_variable->next_time);
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					//fprintf(workerout, "MyCurrentQueue =  %d\n", MyCurrentQueue);
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED) {
						std::cerr << "Queue skipped" << std::endl;
						exit(1);
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					CalculateAcceleration01(ptcl);
					//MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					MyCurrentQueue += NumberOfWorker;
					//fprintf(workerout, "PID= %d,CurrentTimeIrr =  %e\n", ptcl_id, ptcl->CurrentTimeIrr);
				}
				//fprintf(workerout, "InitAcc1 ends\n--------------------------\n");
				fflush(workerout);
				//std::cout << "Processor " << MyRank<< " done." << std::endl;
				//MPI_Waitall(NumberOfCommunication, requests, statuses);
				//NumberOfCommunication = 0;
				break;


/*=======================================*
 *           Initialize Acceleration 2   *
 *=======================================*/
			case InitAcc2: // Initialize Acceleration(23)
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED) {
						std::cerr << "Queue skipped" << std::endl;
						exit(1);
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					CalculateAcceleration23(ptcl);
					//MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					MyCurrentQueue += NumberOfWorker;
				} while (MyCurrentQueue < global_variable->QueueSize);
				break;


/*=======================================*
 *           Initialize Time             *
 *=======================================*/
			case InitTime: // Initialize Time Step
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED) {
						std::cerr << "Queue skipped" << std::endl;
						exit(1);
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					if (ptcl->isActive)
						ptcl->initializeTimeStep();
					//MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					//MyCurrentQueue++;
					MyCurrentQueue += NumberOfWorker;
				}
				break;


/*=======================================*
 *           Time Synchronization        *
 *=======================================*/
			case TimeSync: // Initialize Timestep variables
				MPI_Win_fence(0, win);
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				broadcastFromRoot(time_block);
				broadcastFromRoot(block_max);
				broadcastFromRoot(time_step);
				//MPI_Win_sync(win);  // Synchronize memory
				//MPI_Barrier(MPI_COMM_WORLD);
				//MPI_Win_fence(0, win);
				fprintf(stderr, "(%d) nbody+:time_block = %d, EnzoTimeStep=%e\n", MyRank, time_block, EnzoTimeStep);
				fflush(stderr);
				break;

#ifdef FEWBODY
/*=======================================*
 *           Search Primordial FewBody   *
 *=======================================*/
			case SearchPrimordialGroup: // Primordial binary search
		MPI_Win_fence(0, win4);
		MPI_Win_fence(0, win5);
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED) {
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					if (ptcl->isActive)
						ptcl->initializeTimeStep();
					//MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					//MyCurrentQueue++;
					MyCurrentQueue += NumberOfWorker;
				}
				break;

			case SearchGroup: // Few-body group search
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				MPI_Win_sync(win4);
				MPI_Win_sync(win5);
				//fprintf(workerout, "SearchGroup starts, next_time =  %e\n", global_variable->next_time);
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					ptcl_id = queues[MyCurrentQueue];
					//fprintf(workerout, "MyCurrentQueue =  %d\n", MyCurrentQueue);
					if (ptcl_id == ALREADY_COMPUTED) {
						MyCurrentQueue++;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					//std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					ptcl->NewNumberOfNeighbor = 0;

					if (ptcl->getBinaryInterruptState() == BinaryInterruptState::manybody)
					{
						ptcl->checkNewGroup2();
						ptcl->setBinaryInterruptState(BinaryInterruptState::none);
					}
					else
					{
						if (ptcl->TimeStepIrr * EnzoTimeStep * 1e4 < TSEARCH)
							ptcl->checkNewGroup();
					}

					// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					//MPI_Send(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD);
					MyCurrentQueue += NumberOfWorker;
					//fprintf(workerout, "PID= %d,CurrentTimeIrr =  %e\n", ptcl_id, ptcl->CurrentTimeIrr);
				}
				//fprintf(workerout, "SearchGroup ends\n--------------------------\n");
				//fflush(workerout);
				// std::cerr << "FB search of particle  " << ptcl_id << " is successfully finished on rank " << MyRank << "." <<std::endl;
				break;

			case MakePrimordialGroup: // Make a primordial group
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
				MPI_Recv(&ptcl_id  , 1, MPI_INT   , ROOT, PTCL_TAG, MPI_COMM_WORLD, &status);
				ptcl = &particles[ptcl_id];
				std::cout << "Primordial FewBody object of particle " << ptcl->PID
						  << " is successfully initialized on rank " << MyRank << "." <<std::endl;
				makePrimordialGroup(ptcl);

				break;

			case MakeGroup: // Make a group
				MPI_Win_sync(win4);
				MPI_Win_sync(win5);
				ptcl_id = queues[MyCurrentQueue];
				ptcl = &particles[ptcl_id];
				NewFBInitialization(ptcl);
				std::cout << "FewBody object of particle " << ptcl->PID
						  << " is successfully initialized on rank " << MyRank << "." <<std::endl;
				ptcl->GroupInfo->MyRank = MyRank;
				// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
				MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
				break;

			case DeleteGroup: // Delete a Group struct
		MPI_Win_sync(win4);
		MPI_Win_sync(win5);
				ptcl_id = queues[MyCurrentQueue];
				ptcl = &particles[ptcl_id];
				deleteGroup(ptcl);
				std::cout << "deleting CM particle " << ptcl->PID
						  << " is successfully initialized on rank " << MyRank << "." <<std::endl;
				// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
				MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
				break;

			case ARIntegration: // SDAR for few body encounters
				MPI_Win_fence(0, win4);
				MPI_Win_fence(0, win5);
#ifdef IRR_CM_SEPARATE
				next_time = global_variable->next_time;
				while (MyCurrentQueue < global_variable->QueueSize)
				{
					ptcl_id = queues[MyCurrentQueue];
					if (ptcl_id == ALREADY_COMPUTED)
					{
						MyCurrentQueue++;
						continue;
					}
					std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << std::endl;
					ptcl = &particles[ptcl_id];
					if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr)
					{
						fprintf(stderr, "Rank = %d / CM PID = %d\n", MyRank, ptcl_id);
						fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
						exit(EXIT_FAILURE);
					}
					std::cout << "Processor " << MyRank << ": CM PID= " << ptcl->GroupInfo->MyRank << std::endl;
					if (ptcl->GroupInfo->MyRank != MyRank)
					{
						MyCurrentQueue++;
						//MyCurrentQueue += NumberOfWorker;
						continue;
					}
					queues[MyCurrentQueue] = ALREADY_COMPUTED;
					MPI_Win_sync(win4);

					std::cout << "Processor " << MyRank << ": PID= " << ptcl_id << ". SDAR starts." << std::endl;

					ptcl->GroupInfo->ARIntegration(next_time);
					if (!ptcl->GroupInfo->isMerger && !ptcl->GroupInfo->isTerminate)
						ptcl->GroupInfo->isTerminate = ptcl->GroupInfo->CheckBreak();

					if (ptcl->GroupInfo->isTerminate)
					{
						FBdeleteGroup(ptcl->GroupInfo);
						std::cout << "(SDAR) Processor " << MyRank << ": PID= " << ptcl->PID << " deleted!" << std::endl;
						std::cout << "isActive=" << ptcl->isActive << std::endl;
					}
					else
						std::cout << "(SDAR) Processor " << MyRank << ": PID= " << ptcl->PID << " done!" << std::endl;
					// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
					//MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
					MPI_Win_sync(win);
					MPI_Win_flush_all(win);
					MPI_Send(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD);
					MyCurrentQueue++;
					//MyCurrentQueue += NumberOfWorker;
				}
#else
#endif
				break;
			
			case MergeManyBody: // Merger insided many-body (>2) group
		MPI_Win_sync(win4);
		MPI_Win_sync(win5);
				ptcl_id = queues[MyCurrentQueue];
				ptcl = &particles[ptcl_id];
				std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << std::endl;

				if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
					exit(EXIT_FAILURE);
				}

				NewFBInitialization3(ptcl->GroupInfo);

				ptcl->GroupInfo->isMerger = false;
				ptcl->setBinaryInterruptState(BinaryInterruptState::none);

				// MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[NumberOfCommunication++]);
				//MPI_Isend(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD, &requests[0]);
				MPI_Send(NULL, 0, MPI_BYTE, ROOT, FINISH_TAG, MPI_COMM_WORLD);
				std::cout << "(SDAR) Processor " << MyRank<< ": PID= "<<ptcl->PID << " NewFBInitialization3 done!" <<std::endl;
				break;
#endif 

			case Ends: // Simualtion ends
				std::cout << "Processor " << MyRank<< " returns." << std::endl;
				return;
				break;

			case Error:
				perror("Error task assignments");
				exit(EXIT_FAILURE);
				break;

			default:
				break;
		}

		MPI_Win_fence(0, win);
        //MPI_Barrier(MPI_COMM_WORLD);
		//MPI_Win_fence(0, win);
		//MPI_Win_sync(win);
		//MPI_Win_flush_all(win);
		//std::cerr << "Processor " << MyRank << " done." << std::endl;
	}
}


void ComputeAcceleration(int ptcl_id, double next_time) {
	Particle *ptcl = &particles[ptcl_id];
	Particle *neighbor;
	double neighbor_pos[Dim], neighbor_vel[Dim];
	for (int i=0; i<ptcl->NumberOfNeighbor; i++) {
		neighbor = &particles[ptcl->Neighbors[i]];
		//neighbor->predictParticle(next_time, neighbor_pos, neighbor_vel);
		//ptcl->getAcceleration(neighbor_pos, neighbor_vel);
	}
}

