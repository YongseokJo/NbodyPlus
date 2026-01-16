#ifndef QUEUE_SCHEDULER_H
#define QUEUE_SCHEDULER_H
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iomanip>
#include "global.h"
#include "Worker.h"


class QueueScheduler
{
public:
    std::unordered_set<Worker*> WorkersToGo;

    QueueScheduler()
    {
        _FreeWorkers.reserve(NumberOfWorker);
        WorkersToGo.reserve(NumberOfWorker);
    }

    void initialize(TaskName task, double next_time)
    {
        _initialize();
        _task = task;
        _next_time = next_time;
    }

    void initialize(TaskName task)
    {
        _initialize();
        _task = task;
        _next_time = -1.;
    }

    void takeQueue(std::vector<int> &queue_list) {
        _queue_list = queue_list;
        _total_queues = _queue_list.size();
    }


    void assignQueueAuto() {
        if (_queue_list.size() == 0)
            return;
        for (auto worker = _FreeWorkers.begin(); worker != _FreeWorkers.end();)
        {
            if (_queue_list.size() > 0 && (*worker)->NumberOfQueues == 0)
            {
                _queue.task = _task;
                _queue.pid = _queue_list.back();
                _queue_list.pop_back();
                _queue.next_time = _next_time;
                (*worker)->addQueue(_queue);
                WorkersToGo.insert(*worker);
                _assigned_queues++;
                worker = _FreeWorkers.erase(worker);
                //_queue.print();
            }
            else {
               ++worker; 
            }
        }
    }






    void runQueueAuto() {
        for (auto worker = WorkersToGo.begin(); worker != WorkersToGo.end();)
        {
            if ((*worker)->NumberOfQueues > 0 && !(*worker)->onDuty)
            {
                (*worker)->runQueue();
                worker = WorkersToGo.erase(worker);
            }
            else
            {
                ++worker;
            }
        }
    }

    void sendQueueforRegCuda(Worker *worker) {
        worker->sendTask(_queue);
    }


    void takeQueueRegularList(std::unordered_set<int> &queue_list) {
        _queue_list_ = queue_list;
        _total_queues = _queue_list_.size();
    }

    void assignQueueAutoRegularList() {
        if (_queue_list_.size() == 0)
            return;
        for (auto worker = _FreeWorkers.begin(); worker != _FreeWorkers.end();)
        {
            if (_queue_list_.size() > 0 && (*worker)->NumberOfQueues == 0)
            {
                auto it = _queue_list_.begin();
                _queue.task = _task;
                _queue.pid = *it;
                _queue_list_.erase(it);
                _queue.next_time = _next_time;
                (*worker)->addQueue(_queue);
                WorkersToGo.insert(*worker);
                _assigned_queues++;
                worker = _FreeWorkers.erase(worker);
                //_queue.print();
            }
            else {
               ++worker; 
            }
        }
    }



    Worker* waitQueue(int type) {
        if (_total_queues == 0) // in case of no queue
            return nullptr;

        if (type == 0) // blocking
        {
            PROFILE_START(TimerID::QueueWait);
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &_status);
            PROFILE_STOP(TimerID::QueueWait);
            _completed_queues++;
            // Retrieve the rank of the source processor
            _rank = _status.MPI_SOURCE;
            //fprintf(stdout, "returned rank = %d\n",_rank);
            workers[_rank].callback();
            if (workers[_rank].NumberOfQueues > 0)
                WorkersToGo.insert(&workers[_rank]);
            else
                _FreeWorkers.insert(&workers[_rank]);
            return &workers[_rank];
        }
        else if (type == 1) // non-blocking
        {
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &_flag, &_status);
            if (_flag)
            {
                _rank = _status.MPI_SOURCE;
                //fprintf(stdout,"Iprobe probed someting on rank %d.\n", _rank);
                return &workers[_rank];
            }
            else
            {
                //fprintf(stdout,"Iprobe probed nothing yet.\n");
                return nullptr;
            }
        }
        return nullptr;  // Safety return for undefined type
    }


    // this is only for a non-blocking wait.
    void callback(Worker* worker) {
        if (worker->getCurrentQueue()->task == ARIntegration) 
            _completed_cm_queues++;
        worker->callback();
        _completed_queues++;
        if (worker->NumberOfQueues > 0)
            WorkersToGo.insert(worker);
        else
            _FreeWorkers.insert(worker);
    }




    bool isComplete() {
        if (_completed_queues == _total_queues && _assigned_queues == _total_queues)
            return _complete; // but returns false
        else
            return _not_complete; 
    }


    void printFreeWorker() {
        std::cout << std::left << "FreeWorker List: ";
        for (Worker* worker: _FreeWorkers) {
            std::cout << std::setw(3) << worker->MyRank;
        }
        std::cout << std::endl;
    }

    void printWorkerToGo() {
        std::cout << std::left << "WorkersToGo List: ";
        for (Worker* worker: WorkersToGo) {
            std::cout << std::setw(3) << worker->MyRank;
        }
        std::cout << std::endl;
    }

    void printStatus() {
        std::cout << "-----------Queue Status-----------" << std::endl;
        std::cout << "Total Queues = " << _total_queues << std::endl;
        std::cout << "Assigned Queues = " << _assigned_queues << std::endl;
        std::cout << "Completed Queues = " << _completed_queues << std::endl;
        std::cout << "Total CM Queues = " << _total_cm_queues << std::endl;
        std::cout << "Completed CM Queues = " << _completed_cm_queues << std::endl;

        std::cout << "-----------Worker Status-----------" << std::endl;
        std::cout << std::left << std::setw(10) << "MyRank";
        for (int i=1; i<=NumberOfWorker; i++) {
            std::cout << "|  " << std::setw(4) << workers[i].MyRank;
        }
        std::cout << std::endl;

        std::cout << std::left << std::setw(10) << "OnDuty";
        for (int i=1; i<=NumberOfWorker; i++) {
            std::cout << "|  " << std::setw(4) << workers[i].onDuty ? "T" : "F" ;
        }
        std::cout << std::endl;

        std::cout << std::left << std::setw(10) << "#ofQueue";
        for (int i=1; i<=NumberOfWorker; i++) {
            std::cout << "|  " << std::setw(4) << workers[i].NumberOfQueues;
        }
        std::cout << std::endl;
        std::cout << "-----------Details-----------" << std::endl;
        printFreeWorker();
        printWorkerToGo();
    }

    void setTotalQueue(int total_queues) {_total_queues = total_queues;}

    void assignWorker(Worker *worker) {
        if (!(worker->onDuty)) {
            WorkersToGo.insert(worker);
            _FreeWorkers.erase(worker);
        }
        _assigned_queues++;
    }

#ifdef FEWBODY
    std::unordered_set<int> CMPtcls;
    //std::unordered_map<int, int> CM_worker;
    //std::vector<int> UpdatedCMPtcl;
    //int fb_total_tasks;    // this is for the SDAR computations by YS 2025.01.06
    //int fb_assigned_tasks; // this is for the SDAR computations by YS 2025.01.06

    void initializeIrr(TaskName task, double next_time, std::vector<int> &queue_list)
    {
        _initialize();
        int pid;
        _task = task;
        _next_time = next_time;
        _total_queues = queue_list.size();
        _queue_list.reserve(_total_queues);
        _queue_list.resize(_total_queues);
        int single_ptcl=0, cm_ptcl=_total_queues-1;
        _total_cm_queues=0;
        _completed_cm_queues=0;

        for (int i=0; i<_total_queues; i++)
        {
            pid = queue_list[i];
            particles[pid].isUpdateToDate = false;
            if (!particles[pid].isCMptcl)
            {
                _queue_list[single_ptcl++] = pid;
            }
            else
            {
                _queue_list[cm_ptcl--] = pid;
                CMPtcls.insert(pid);
                _total_cm_queues++;
            }
        }
        _total_queues += _total_cm_queues;
    }

    /* check the neighbors of CM particles if they're up to date
    void checkCMPtclToGo() {
        Particle *ptcl;
        int i = 0;
        for (int pid:UpdatedCMPtcl) {
           ptcl = &particles[pid]; 
           if (!ptcl->isCMPtcl) {
            fprintf(stderr, "this is not CM ptcl!\n");
            exit(EXIT_FAILURE);
           }
           for (int j = 0; j < _ptcl->NumberOfNeighbor; j++)
           {
               if (particles[_ptcl->Neighbors[j]].isUpdateToDate == false)
                   break;
           }
           ReadyToGoCMPtcl.insert(pid); // (Query to myself) should I add worker rank?
        }
    }
    */



    ~QueueScheduler() {
    }
#endif

private:
    std::unordered_set<Worker*> _FreeWorkers;
    std::vector<int> _queue_list;
    std::unordered_set<int> _queue_list_;
    Queue _queue;
    TaskName _task;
    int _rank, _flag;
    double _next_time;
    int _total_queues, _assigned_queues, _completed_queues;
    int _total_cm_queues, _completed_cm_queues;
    const bool _complete = false;
    const bool _not_complete = true;
    MPI_Status _status;    // Pointer to the status object
    MPI_Request _request;    // Pointer to the status object



    void _initialize() {
        WorkersToGo.clear();
        _FreeWorkers.clear();
        for (int i = 1; i <= NumberOfWorker; i++) {
            workers[i].initialize();
            _FreeWorkers.insert(&workers[i]);
        }
        _total_queues=0;
        _assigned_queues=0;
        _completed_queues=0;
    }

#ifdef unuse
        void doSimpleTask(Task task, std::vector<int> &list)
    {
        int return_value, i = 0;
        _task = task;
        _total_tasks = list.size();
        _setTask(_task);
        //fprintf(stdout, "task=%d, the size of FreeWorker is %d\n", _task, _FreeWorkers.size());
        do {
            if (_FreeWorkers.size() == 0 || _assigned_tasks == _total_tasks) {
                // have to add check all the sends are recved.
                // MPI_Waitall(NumberOfCommunication, requests, statuses);
                // NumberOfCommunication = 0;
                _checkCompletion(return_value);
                /* we can do something here */
                _Callback();
                _completed_tasks++;
            }

            if (_FreeWorkers.size() != 0 && _assigned_tasks < _total_tasks) {
                _WorkerTmp = _FreeWorkers.back();
                _FreeWorkers.pop_back();
                _WorkerTmp->PID = list[i++];
                _assignJobs(_WorkerTmp);
                _assigned_tasks++;
                //fprintf(stdout, "assigned_tasks = %d, number of free worker = %d pid = %d rank = %d\n",
                //_assigned_tasks, _FreeWorkers.size(), _WorkerTmp->PID, _WorkerTmp->MyRank);
                fflush(stdout);
            }
        } while(_completed_tasks < _total_tasks);
    }

    void doIrregularForce(std::vector<int> &ParticleList,
                          double &next_time)
    {
        _task=IrrForce;
        _ParticleList = ParticleList;
        _total_tasks = _ParticleList.size();
        std::vector<int> SingleParticleList;
        std::unordered_set<int> CMPtcls;
        std::vector<int> UpdatedCMPtcl;
        std::unordered_set<int> ReadyToGoCMPtcl;
        int pid_tmp, return_value;
        SingleParticleList.reserve(_total_tasks);

        _setTask(_task);

        for (int i=0; i<_ParticleList.size(); i++) {
            pid_tmp = _ParticleList[i];
            particles[pid_tmp].isUpdateToDate = false;
            if (particles[pid_tmp].isCMptcl)
            {
                SingleParticleList.push_back(pid_tmp);
            }
            else
            {
                CMPtcls.insert(pid_tmp);
            }
        }

        int fb_total_tasks = CMPtcls.size(); // this is for the SDAR computations by YS 2025.01.06
        int fb_assigned_tasks = 0; // this is for the SDAR computations by YS 2025.01.06
        UpdatedCMPtcl.resize(fb_total_tasks);

        do {
            if (_FreeWorkers.size() == 0 || _assigned_tasks == _total_tasks) 
            {
                // have to add check all the sends are recved.
                // MPI_Waitall(NumberOfCommunication, requests, statuses);
                // NumberOfCommunication = 0;
                _checkCompletion(return_value);
                // this ensures that cmptcls are ready to go (neighbors are update to date)
                _checkCMPtclNeighbor(UpdatedCMPtcl, ReadyToGoCMPtcl); 
                _WorkerTmp = _Callback();
                _completed_tasks++;

            }

            // SDAR assignment
            if (_FreeWokrers.size() != 0 && ReadyToGoCMPtcl.size() != 0) {
                // I have to add that if the worker is a CMworker, I should run SDAR integration (Query to myself)
                if (particles[_WorkerTmp->PID].isCMptcl && _WorkerTmp->isCMWorker) 
                {
                    // all neighbors are up to date.
                
                
                        fprintf(stderr, "Worker does not match!\n");
                        exit(EXIT_FAILURE);
                    }
                    _FreeWorkers.pop_back();
                    _WorkerTmp->task = 26;
                    _assignJobs(_WorkerTmp);
                    fb_assigned_tasks++;
                }
            }

            if (_FreeWorkers.size() != 0 && _assigned_tasks < _total_tasks)
            {
                _WorkerTmp = _FreeWorkers.back();
                _FreeWorkers.pop_back();

                // CM ptcl assignment
                // if this worker contains CM ptcl, see if any relevant CM ptcls in ParticleList.
                if (_WorkerTmp->isCMWorker)
                {
                    for (int cm_pid : _WorkerTmp->CMPtclIDs)
                    {
                        _WorkerTmp->isCMWorker = false;
                        if (CMPtcls.find(cm_pid) != CMPtcls.end())
                        {
                            _WorkerTmp->PID = cm_pid;
                            _WorkerTmp->isCMWorker = true;
                            CMPtcls.erase(cm_pid);
                            UpdatedCMPtcl.push_back(cm_pid);
                            break;
                        }
                    }
                }

                // Single particle assignment
                if (!_WorkerTmp->isCMWorker)
                {
                    pid_tmp = SingleParticleList.back();
                    SingleParticleList.pop_back();
                    _WorkerTmp->PID = pid_tmp;
                }

                _WorkerTmp->task = 0;
                _WorkerTmp->next_time = next_time;
                //fprintf(stdout, "assigned_tasks = %d/%d, number of free worker = %d pid = %d rank = %d\n",
                //_assigned_tasks, _total_tasks, _FreeWorkers.size(), _WorkerTmp->PID, _WorkerTmp->MyRank);
                _assignJobs(_WorkerTmp);
                _assigned_tasks++;
            }
        } while (_completed_tasks < _total_tasks && fb_completed_tasks < fb_total_tasks);

        if (SingleParticleList.size() != 0 || CMPtcls.size() != 0) {
            fprintf(stderr, "Something's wrong. The particle list is not empty.\n");
            exit(1);
        }

        UpdatedCMPtcl.clear();
        UpdatedCMPtcl.shrink_to_fit();
        SingleParticleList.clear();
        SingleParticleList.shrink_to_fit();
    }

//private:
    std::vector<int> _ParticleList;
    int _total_tasks, _completed_tasks, _assigned_tasks, _remaining_tasks;
    int  _pid;
    Task _task;
    std::vector<Worker*> _FreeWorkers;
    MPI_Request _request;  // Pointer to the request handle
    MPI_Status _status;    // Pointer to the status object
    int _rank;
    Particle* _ptcl;
    Worker* _WorkerTmp, _CompletedWorker;



    void _assignJobs(Worker *worker) {
        if ((worker->task == 0) || (worker->task == 1) || (worker->task == 26))
        {
            MPI_Send(&worker->task,      1, MPI_INT,    worker->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->PID,       1, MPI_INT,    worker->MyRank, PTCL_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->next_time, 1, MPI_DOUBLE, worker->MyRank, TIME_TAG, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&worker->task, 1, MPI_INT, worker->MyRank, TASK_TAG, MPI_COMM_WORLD);
            MPI_Send(&worker->PID,  1, MPI_INT, worker->MyRank, PTCL_TAG, MPI_COMM_WORLD);
        }
        worker->onDuty = true;
    }


    void _checkCompletion(int &return_value) {
        MPI_Irecv(&return_value, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD, &_request);
    }


    Worker* _Callback() {
        MPI_Wait(&_request, &_status);
        int rank = _status.MPI_SOURCE;
        if (!workers[rank-1].onDuty) {
            fprintf(stderr, "Something's worng! the worker was not on duty.1\n");
            exit(1);
        }
        workers[rank-1].onDuty = false;
        _FreeWorkers.push_back(&workers[rank-1]);
        return &workers[rank-1];
    }


    void _setTask(Task task) {
        for (int i = 0; i < NumberOfWorker; i++) {
            workers[i].task = task;
        }
    }

    void _checkCMPtclNeighbor(const std::vector<int>& UpdatedCMPtcl, std::unordered_set<int>& ReadyToGoCMPtcl) {
        Particle *ptcl;
        int i = 0;
        for (int pid:UpdatedCMPtcl) {
           ptcl = &particles[pid]; 
           if (!ptcl->isCMPtcl) {
            fprintf(stderr, "this is not CM ptcl!\n");
            exit(EXIT_FAILURE);
           }
           for (int j = 0; j < _ptcl->NumberOfNeighbor; j++)
           {
               if (particles[_ptcl->Neighbors[j]].isUpdateToDate == false)
                   break;
           }
           ReadyToGoCMPtcl.insert(pid); // (Query to myself) should I add worker rank?
        }
    }
#endif
};




#endif
