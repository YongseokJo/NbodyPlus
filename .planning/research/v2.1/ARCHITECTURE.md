# Architecture Research: Async MPI Integration

## Current Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                        ROOT (Rank 0)                         │
│  ┌─────────────────┐  ┌─────────────────┐  ┌──────────────┐ │
│  │ QueueScheduler  │  │ IrregularRoutines│  │   Workers[]  │ │
│  │                 │  │                 │  │              │ │
│  │ - assignQueue() │  │ - main loop     │  │ - send_task()│ │
│  │ - runQueue()    │  │ - skiplist      │  │ - callback() │ │
│  │ - waitQueue()   │──│ - CM handling   │──│ - on_duty    │ │
│  └─────────────────┘  └─────────────────┘  └──────────────┘ │
│           │                    │                   │         │
│           └────────────────────┼───────────────────┘         │
│                                ▼                             │
│                    ┌─────────────────────┐                   │
│                    │  Blocking MPI       │                   │
│                    │  MPI_Send / MPI_Recv│                   │
│                    │  MPI_Probe          │                   │
│                    └─────────────────────┘                   │
└─────────────────────────────────────────────────────────────┘
                              │
              ┌───────────────┼───────────────┐
              ▼               ▼               ▼
        ┌──────────┐   ┌──────────┐   ┌──────────┐
        │ Worker 1 │   │ Worker 2 │   │ Worker N │
        │ MPI_Recv │   │ MPI_Recv │   │ MPI_Recv │
        │ compute  │   │ compute  │   │ compute  │
        │ MPI_Send │   │ MPI_Send │   │ MPI_Send │
        └──────────┘   └──────────┘   └──────────┘
```

### Current Flow (Blocking)

1. Root: `assignQueueAuto()` — picks particle, assigns to free worker
2. Root: `runQueueAuto()` → `Worker::send_task()` → `MPI_Send`
3. Root: `waitQueue()` → `MPI_Probe` (blocks until any worker done)
4. Root: `callback()` → `MPI_Recv` (blocks to get result)
5. Repeat for each particle

**Problem**: Root sends 1 task, waits for completion, sends next. No overlap.

## Target Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                        ROOT (Rank 0)                         │
│  ┌─────────────────┐  ┌─────────────────┐  ┌──────────────┐ │
│  │ QueueScheduler  │  │ AsyncManager    │  │   Workers[]  │ │
│  │                 │  │ (NEW)           │  │              │ │
│  │ - same API      │  │ - send_requests │  │ - same API   │ │
│  │                 │  │ - recv_requests │  │ - async ops  │ │
│  │                 │──│ - pending count │──│ - request[]  │ │
│  └─────────────────┘  └─────────────────┘  └──────────────┘ │
│           │                    │                   │         │
│           └────────────────────┼───────────────────┘         │
│                                ▼                             │
│                    ┌─────────────────────┐                   │
│                    │  Non-Blocking MPI   │                   │
│                    │  MPI_Isend/Irecv    │                   │
│                    │  MPI_Waitany/Testany│                   │
│                    └─────────────────────┘                   │
└─────────────────────────────────────────────────────────────┘
```

### Target Flow (Async)

1. Root: Post `MPI_Irecv` for all workers (pre-posted receives)
2. Root: Send task via `MPI_Isend`, don't wait
3. Root: Send more tasks (fill worker pipeline)
4. Root: `MPI_Waitany` — process first completion
5. Root: Re-post receive, send next task
6. Overlap: local work between sends/waits

## Component Changes

### 1. Worker struct (`src/worker.h`)

**Current:**
```cpp
void send_task(Queue& queue) {
    MPI_Send(&queue, 1, queue_type_mpi, this->rank, QUEUE_TAG, MPI_COMM_WORLD);
    on_duty = true;
}

void callback() {
    MPI_Recv(&return_value, 1, MPI_INT, this->rank, TERMINATE_TAG, ...);
    on_duty = false;
    // ...
}
```

**After:**
```cpp
MPI_Request send_request;
MPI_Request recv_request;
int result_buffer;  // Dedicated buffer per worker

void send_task_async(Queue& queue) {
    MPI_Isend(&queue, 1, queue_type_mpi, this->rank, QUEUE_TAG, ..., &send_request);
    on_duty = true;
}

void post_receive() {
    MPI_Irecv(&result_buffer, 1, MPI_INT, this->rank, TERMINATE_TAG, ..., &recv_request);
}

// callback() becomes simpler - just process result_buffer
```

### 2. QueueScheduler (`src/queue_scheduler.h`)

**Current:**
```cpp
Worker* waitQueue(int type) {
    if (type == 0) {  // blocking
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &_status);
        workers[_rank].callback();
        // ...
    }
}
```

**After:**
```cpp
// Manage array of recv requests
MPI_Request recv_requests[MAX_WORKERS];

Worker* waitQueueAsync() {
    int index;
    MPI_Waitany(num_workers, recv_requests, &index, MPI_STATUS_IGNORE);
    Worker* worker = &workers[index + 1];  // rank offset
    // Process completion
    // Re-post receive
    return worker;
}

Worker* testQueueAsync() {
    int index, flag;
    MPI_Testany(num_workers, recv_requests, &index, &flag, MPI_STATUS_IGNORE);
    if (flag) {
        // Same as waitQueueAsync
    }
    return flag ? &workers[index + 1] : nullptr;
}
```

### 3. IrregularRoutines (`src/irregular_routines.cpp`)

**Current main loop:**
```cpp
do {
    queue_scheduler.assignQueueAuto();
    queue_scheduler.runQueueAuto();
    do {
        worker = queue_scheduler.waitQueue(1);
        // CM particle handling...
    } while (worker == nullptr);
    queue_scheduler.callback(worker);
} while (queue_scheduler.isComplete());
```

**After:**
```cpp
// Initial: post all receives
queue_scheduler.postAllReceives();

do {
    // Send tasks to all free workers (no waiting)
    queue_scheduler.assignAndSendAll();
    
    // Process completions as they arrive
    while ((worker = queue_scheduler.waitQueueAsync()) != nullptr) {
        queue_scheduler.processCompletion(worker);
        // CM particle handling...
        if (more_tasks) {
            queue_scheduler.sendNextTask(worker);
        }
    }
} while (queue_scheduler.isComplete());
```

## Build Order (Phases)

### Phase 11: Async Infrastructure
1. Add request arrays to Worker struct
2. Add `send_task_async()` and `post_receive()`
3. Add profiler timers for async operations

### Phase 12: QueueScheduler Async
1. Add `waitQueueAsync()` using `MPI_Waitany`
2. Add `testQueueAsync()` using `MPI_Testany`
3. Maintain backward compatibility (keep blocking methods)

### Phase 13: IrregularRoutines Integration
1. Modify main loop to use async methods
2. Handle CM particle dependencies correctly
3. Validate correctness (energy conservation)

### Phase 14: Measurement & Tuning
1. Profile async vs blocking
2. Measure overlap effectiveness
3. Tune pipeline depth

## Data Flow Changes

### Message Ordering

**Before:** Strict send-wait-send-wait ordering, messages naturally ordered.

**After:** Multiple in-flight messages, order determined by `MPI_Waitany` completion.

**Impact:** 
- Particle processing order may change within a timestep
- CM particle dependency logic must still work
- Physics should be deterministic (same results for same input)

### Buffer Management

**Before:** Reuse same buffer after blocking send completes.

**After:** Need separate buffers for each in-flight message, or wait before reuse.

**Strategy:**
- Each Worker has its own `Queue` buffer for sending
- Each Worker has its own `result_buffer` for receiving
- No buffer sharing between workers

## Profiler Integration

Add new timers:
- `TimerID::MPIIsend` — Time in MPI_Isend
- `TimerID::MPIIrecv` — Time in MPI_Irecv  
- `TimerID::MPIWaitany` — Time in MPI_Waitany
- `TimerID::MPITestany` — Time in MPI_Testany (may be very small)
- `TimerID::AsyncOverlap` — Time doing local work during async window

## Sources

- [LLNL HPC Tutorials](https://hpc-tutorials.llnl.gov/mpi/non_blocking/)
- [Towards Efficient HPC: Exploring Overlap Strategies](https://www.mdpi.com/2227-7390/13/11/1848)
- Current codebase: `src/queue_scheduler.h`, `src/worker.h`, `src/irregular_routines.cpp`

---
*Research completed: 2026-01-18*
