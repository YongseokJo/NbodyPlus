# Phase 2 Research: MPI Integration

## Current MPI Implementation Analysis

### Existing Structure (mpi_routines.cpp)

The current implementation uses **4 MPI windows** for shared memory:

| Window | Variable | Size | Purpose |
|--------|----------|------|---------|
| `win` | `particles` | `sizeof(Particle) * MAX_NUM_PARTICLE` | Main particle AoS array |
| `win2` | `g_state` | `sizeof(GlobalVariable)` | Global simulation state |
| `win3` | `neighbors` | `sizeof(int) * MAX_NUM_PARTICLE * MAX_NUM_NEIGHBOR` | Neighbor indices |
| `win4` | `new_neighbors` | `sizeof(int) * MAX_NUM_PARTICLE * MAX_NUM_NEIGHBOR` | New neighbor indices |

### Allocation Pattern (Lines 51-61)

```cpp
// Root (shared_rank == 0) allocates actual memory
if (shared_rank == 0) {
    MPI_Win_allocate_shared(sizeof(Particle) * MAX_NUM_PARTICLE, sizeof(Particle),
                            MPI_INFO_NULL, shared_comm, &particles, &win);
    // ... more allocations
} else {
    // Workers allocate 0 bytes but get pointer via shared_query
    MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, shared_comm,
                            &particles, &win);
}

// All ranks query to get the actual pointer
MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
```

### Synchronization Pattern (Lines 142-143)

```cpp
MPI_Win_sync(win);  // Memory synchronization
MPI_Barrier(shared_comm);  // Process synchronization
```

---

## SoA Migration Strategy

### Window Count

From Phase 1 research, ParticleData has **66 arrays**:
- 50 double arrays (43 critical + 5 important + 2 low priority)
- 6 ull_t arrays
- 8 int arrays
- 3 bool arrays

Plus existing non-particle windows:
- `win2` (g_state) — keep as-is
- `win3` (neighbors) — keep as-is
- `win4` (new_neighbors) — keep as-is

**Total windows needed: 66 (SoA) + 3 (existing) = 69 windows**

### Window Organization

Group windows by type for cleaner management:

| Type | Count | Example Arrays |
|------|-------|----------------|
| double | 50 | pos_x, vel_y, mass, acc_total[0][0], ... |
| ull_t | 6 | current_block_irr, time_block_reg, ... |
| int | 8 | num_neighbors, pid, particle_type, ... |
| bool | 3 | is_active, is_up_to_date, is_cm_particle |

### MPI_Win Array Declaration

```cpp
// In global.h or particle_data_mpi.h
// Window arrays organized by type
MPI_Win win_double[50];   // All double arrays
MPI_Win win_ull[6];       // All ull_t arrays
MPI_Win win_int[8];       // All int arrays
MPI_Win win_bool[3];      // All bool arrays

// Or flat array with indices
MPI_Win win_particle_data[66];
```

---

## Implementation Patterns

### Pattern 1: Separate MPI Allocation Function

Create `particle_data_mpi.cpp` with allocation/deallocation:

```cpp
void ParticleDataMPI::allocate_shared(size_t capacity, MPI_Comm comm) {
    int shared_rank;
    MPI_Comm_rank(comm, &shared_rank);

    size_t padded = padded_capacity(capacity);

    // Double arrays (example for pos_x)
    if (shared_rank == 0) {
        MPI_Win_allocate_shared(sizeof(double) * padded, sizeof(double),
                                MPI_INFO_NULL, comm, &pos_x_, &win_pos_x_);
    } else {
        MPI_Win_allocate_shared(0, sizeof(double),
                                MPI_INFO_NULL, comm, &pos_x_, &win_pos_x_);
    }
    MPI_Win_shared_query(win_pos_x_, 0, &size, &disp, &pos_x_);

    // ... repeat for all 66 arrays
}

void ParticleDataMPI::deallocate_shared() {
    MPI_Win_free(&win_pos_x_);
    // ... repeat for all 66 windows
}
```

### Pattern 2: Loop-Based Allocation

Use arrays of pointers and windows for cleaner code:

```cpp
// Array descriptors
struct ArrayDesc {
    void** ptr;        // Pointer to store allocated address
    MPI_Win* win;      // Window handle
    size_t elem_size;  // Element size (sizeof(double), etc.)
};

ArrayDesc double_arrays[] = {
    {(void**)&pos_x_, &win_double[0], sizeof(double)},
    {(void**)&pos_y_, &win_double[1], sizeof(double)},
    // ...
};

// Single loop handles all
for (auto& desc : double_arrays) {
    allocate_shared_array(desc.ptr, desc.win, desc.elem_size, capacity, comm);
}
```

---

## Synchronization Strategy

### Current Sync Point (ParticleSynchronization)

```cpp
void ParticleSynchronization() {
    InitialAssignmentOfTasks(task, num_workers, TASK_TAG);
    MPI_Win_sync(win);  // Sync particle data
    MPI_Barrier(shared_comm);
    // ... wait for workers
}
```

### SoA Sync Pattern

Option A: Sync all windows explicitly
```cpp
void SyncAllParticleWindows() {
    for (int i = 0; i < NUM_DOUBLE_ARRAYS; i++) {
        MPI_Win_sync(win_double[i]);
    }
    for (int i = 0; i < NUM_ULL_ARRAYS; i++) {
        MPI_Win_sync(win_ull[i]);
    }
    // ... etc
    MPI_Barrier(shared_comm);
}
```

Option B: Lock/Unlock for explicit synchronization
```cpp
void LockAllParticleWindows(int lock_type) {
    for (auto& win : all_windows) {
        MPI_Win_lock_all(lock_type, win);
    }
}

void UnlockAllParticleWindows() {
    for (auto& win : all_windows) {
        MPI_Win_unlock_all(win);
    }
}
```

**Recommendation**: Use `MPI_Win_sync()` + `MPI_Barrier()` pattern (matches existing code). Add helper function to sync all windows.

---

## Risk Assessment: Window Explosion (Pitfall #6)

### Concern

69 MPI windows vs current 4 windows:
- Memory overhead per window (internal MPI bookkeeping)
- Synchronization overhead (calling MPI_Win_sync 69 times)
- Creation time (69 MPI_Win_allocate_shared calls)

### Mitigation Strategies

1. **Measure creation time**: Add timing around window allocation
2. **Test with target MPI**: OpenMPI on HPC cluster (check module)
3. **Batch synchronization**: Helper function to sync all windows efficiently
4. **Alternative if problematic**: Single large allocation, manual offset calculation (loses some SoA benefits)

### Expected Impact

Based on MPI-3 shared memory design:
- Window creation is one-time startup cost
- Per-window memory overhead is small (typically < 1KB per window)
- Synchronization scales linearly (69x overhead vs 1x, but sync is fast)

**Recommendation**: Proceed with 69-window approach. Measure creation time. Fall back to combined allocation only if creation time exceeds 1 second on target system.

---

## Integration with ParticleData Class

### Option A: Subclass

```cpp
class ParticleDataMPI : public ParticleData {
public:
    void allocate_shared(size_t capacity, MPI_Comm comm);
    void deallocate_shared();
    void sync_all();

private:
    MPI_Win win_double_[50];
    MPI_Win win_ull_[6];
    MPI_Win win_int_[8];
    MPI_Win win_bool_[3];
    MPI_Comm comm_;
};
```

### Option B: Composition

```cpp
class ParticleDataMPI {
public:
    ParticleData& data() { return data_; }
    void allocate_shared(...);
    void sync_all();

private:
    ParticleData data_;  // Uses MPI-allocated pointers
    MPI_Win windows_[66];
};
```

### Option C: Template Strategy

```cpp
template<typename Allocator>
class ParticleData {
    // Allocator provides allocate/deallocate
};

using ParticleDataStandard = ParticleData<StdAllocator>;
using ParticleDataMPI = ParticleData<MPISharedAllocator>;
```

**Recommendation**: Option A (Subclass) — simplest, keeps ParticleData usable standalone for testing, MPI version inherits accessors.

---

## Files to Create/Modify

### New Files

1. `src/particle_data_mpi.h` — MPI-aware ParticleData declaration
   - Window handles
   - allocate_shared() / deallocate_shared()
   - sync_all()

2. `src/particle_data_mpi.cpp` — Implementation
   - MPI_Win_allocate_shared calls for all 66 arrays
   - MPI_Win_free calls
   - Synchronization helpers

### Modified Files

1. `src/mpi_routines.cpp` — Update initializeMPI()
   - Remove old `MPI_Win_allocate_shared` for `particles` array
   - Call `ParticleDataMPI::allocate_shared()` instead
   - Keep `win2`, `win3`, `win4` for now (g_state, neighbors)

2. `src/global.h` — Update declarations
   - Remove `extern MPI_Win win;`
   - Remove `extern Particle* particles;`
   - Add `extern ParticleDataMPI particle_data;`

---

## Migration Path

### Step 1: Create ParticleDataMPI (non-breaking)
- New files, not yet used
- Test allocation in isolation

### Step 2: Wire into initializeMPI()
- Replace `particles` allocation with `particle_data.allocate_shared()`
- Both old and new exist temporarily

### Step 3: Update synchronization points
- Replace `MPI_Win_sync(win)` with `particle_data.sync_all()`

### Step 4: Remove old Particle array
- Remove `win`, `particles` global
- Update all access points (Phase 4 work)

---

## Dependencies

### Requires from Phase 1
- `ParticleData` class with all 66 arrays ✓ (complete)
- Accessor functions ✓ (complete)
- `padded_capacity()` function ✓ (in particle_data.cpp)

### Required by Later Phases
- Phase 3 (GPU): Will use CPU SoA pointers as source for transfers
- Phase 4 (CPU Routines): Will use accessors
- Phase 6 (I/O): Will use accessors

---

## Questions Resolved

1. **How many windows?** 66 for SoA + 3 existing = 69 total
2. **Class design?** Subclass ParticleData as ParticleDataMPI
3. **Sync pattern?** MPI_Win_sync() + Barrier, helper function for all
4. **Window explosion risk?** Acceptable, measure creation time, fallback available
5. **File organization?** New particle_data_mpi.h/.cpp, modify mpi_routines.cpp

---
*Generated: 2026-01-17*
