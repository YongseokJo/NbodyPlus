# Phase 4 Research: CPU Routines SoA Conversion

## 1. Overview

Phase 4 updates CPU force calculation and integration routines to use the SoA `ParticleData` container created in Phase 1, replacing direct AoS `Particle` struct access.

**Requirements:**
- CPU-01: Update regular force calculation for SoA
- CPU-02: Update irregular force calculation for SoA
- CPU-03: Update prediction routines for SoA
- CPU-04: Update correction routines for SoA
- CPU-05: Update timestep routines for SoA

## 2. Current AoS Access Patterns

### 2.1 Particle Struct Fields Accessed

The following fields are accessed in CPU routines:

| Field | Type | Size | Usage |
|-------|------|------|-------|
| `position[3]` | double[3] | 24B | Position in 3D space |
| `velocity[3]` | double[3] | 24B | Velocity in 3D space |
| `mass` | double | 8B | Particle mass |
| `acc_total[3][4]` | double[3][4] | 96B | Total acceleration + derivatives |
| `acc_regular[3][4]` | double[3][4] | 96B | Regular (far) acceleration |
| `acc_irregular[3][4]` | double[3][4] | 96B | Irregular (near) acceleration |
| `new_position[3]` | double[3] | 24B | Predicted position |
| `new_velocity[3]` | double[3] | 24B | Predicted velocity |
| `neighbor_radius_sq` | double | 8B | Squared neighbor radius |
| `num_neighbors` | int | 4B | Neighbor count |
| `new_num_neighbors` | int | 4B | Updated neighbor count |
| `neighbors_offset` | int | 4B | Offset in neighbor array |
| `current_time_irr/reg` | double | 16B | Current times |
| `time_step_irr/reg` | double | 16B | Time steps |
| `time_block_irr/reg` | ull_t | 16B | Time blocks |
| `time_level_irr/reg` | int | 8B | Time levels |
| `current_block_irr/reg` | ull_t | 16B | Current blocks |
| `new_current_block_irr` | ull_t | 8B | Updated block |
| `pid` | int | 4B | Particle ID |
| `particle_index` | int | 4B | Array index |
| `is_active` | bool | 1B | Active flag |
| `is_cm_particle` | bool | 1B | CM particle flag |
| `cm_particle_index` | int | 4B | Parent CM index |
| `members[10]` | int[10] | 40B | Group member indices |
| `num_members` | int | 4B | Member count |
| `new_members[10]` | int[10] | 40B | New group members |
| `new_num_members` | int | 4B | New member count |

### 2.2 File Analysis

#### compute_acceleration.cpp

Three main functions:
1. **`compute_acceleration_irr()`** — Irregular force calculation
   - Accesses: `this->` for target particle, `ptcl->` for neighbors
   - Fields: position, velocity, mass, acc_irregular, acc_regular, acc_total, new_position, new_velocity, num_neighbors, neighbors_offset, current_time_irr, time_step_irr, current_time_reg, particle_index, is_active, cm_particle_index, new_members, new_num_members

2. **`compute_acceleration_reg()`** — Regular force calculation (CPU-only path)
   - Fields: Same as above plus neighbor_radius_sq, new_neighbors, is_cm_particle, num_members, members

3. **`update_regular_particle_cuda()`** — Post-GPU regular update
   - Processes GPU results, handles neighbor transitions
   - Same field set as regular calculation

#### regular_routines.cpp

Orchestration file — dispatches to:
- `calculateRegAccelerationOnGPU()` when CUDA enabled
- Queue scheduler for CPU path
- Accesses particles array via `particles[index]` pattern

#### irregular_routines.cpp

Complex file with many functions:
- `IrregularRoutines()` — Main irregular loop
- Skip list management
- Few-body handling (FEWBODY code paths)
- Binary detection and termination

Access patterns:
- `particles[ParticleList[i]]` — indexed access
- `ptcl = &particles[id]; ptcl->field` — pointer dereference
- Fields: is_active, current_time_irr, time_step_irr, num_neighbors, current_block_irr, new_current_block_irr, next_block_irr, time_block_reg, current_block_reg, particle_index, is_cm_particle, members, num_members, new_num_members, binary interrupt state, etc.

#### timestep_routines.cpp

Standalone functions (not methods):
- `getNewTimeStepReg(v[3], df[3][4])` — Uses velocity and acceleration derivatives
- `getNewTimeStepIrr(f[3][4], df[3][4])` — Uses total and irregular accelerations
- `getNewTimeStep(f[3][4], df[3][4])` — Generic timestep calculation
- `getBlockTimeStep(dt, TimeLevel, TimeBlock, TimeStep)` — Block quantization

**Note:** These functions take raw arrays, not Particle references. The calling code (in update_particle.cpp) passes `this->velocity`, `this->acc_total`, etc.

#### update_particle.cpp

`Particle::` methods:
- `correct_particle_fourth_order()` — 4th order correction
- `update_radius()` — Neighbor radius adjustment
- `calculate_time_step_irr()` / `_v2()` — Irregular timestep
- `calculate_time_step_reg()` — Regular timestep

Access pattern: All use `this->` to access fields.

## 3. ParticleData SoA Container (Phase 1)

The container from Phase 1 provides:

```cpp
class ParticleData {
    // Position (3 arrays)
    double* pos_x_, *pos_y_, *pos_z_;
    double get_pos_x(size_t i) const;
    void set_pos_x(size_t i, double val);

    // Velocity (3 arrays)
    double* vel_x_, *vel_y_, *vel_z_;

    // Mass (1 array)
    double* mass_;

    // Accelerations (36 arrays: 3 types × 3 dims × 4 orders)
    double* acc_total_[3][4];
    double* acc_reg_[3][4];
    double* acc_irr_[3][4];

    // New pos/vel (6 arrays)
    double* new_pos_x_, *new_pos_y_, *new_pos_z_;
    double* new_vel_x_, *new_vel_y_, *new_vel_z_;

    // Time stepping (10 arrays)
    double* current_time_irr_, *current_time_reg_;
    double* time_step_irr_, *time_step_reg_;
    ull_t* current_block_irr_, *new_current_block_irr_, *current_block_reg_;
    ull_t* time_block_irr_, *time_block_reg_, *next_block_irr_;

    // Neighbor info (3 arrays)
    int* num_neighbors_, *new_num_neighbors_, *neighbors_offset_;

    // IDs and flags (8 arrays)
    int* pid_, *particle_index_, *particle_type_;
    int* time_level_irr_, *time_level_reg_;
    bool* is_active_, *is_up_to_date_, *is_cm_particle_;

    // Other (3 arrays)
    double* neighbor_radius_sq_, *radius_, *delta_mass_;
};
```

## 4. Conversion Strategy

### 4.1 Key Challenge: Method-Based Access

Current code uses Particle methods like `this->position[dim]`. The SoA conversion requires:

1. **Pass particle index explicitly** — Methods need to know which particle
2. **Access via ParticleData** — Replace `this->field` with `data.get_field(i)`

### 4.2 Recommended Approach: Parallel Structures

**Keep Particle struct for now** (Phase 4), but add:
1. Copy data FROM Particle struct TO ParticleData at start of force calculation
2. Operate on ParticleData during hot loops
3. Copy results BACK to Particle struct

This maintains backward compatibility while getting SoA benefits in hot loops.

**Alternative:** Full conversion would require:
- Converting all `Particle` methods to free functions taking `ParticleData&` and index
- Massive refactoring across FewBody, SDAR, etc.
- Risk of breaking energy conservation

### 4.3 Hot Loop Identification

**Most critical loops (O(N²) or O(N×Nneighbor)):**
1. `compute_acceleration_irr()` — inner neighbor loop
2. `compute_acceleration_reg()` — all-particle loop
3. `update_regular_particle_cuda()` — neighbor transition loops

These are where SoA conversion provides the most benefit.

## 5. File Dependencies

```
regular_routines.cpp
    └── compute_acceleration.cpp (via queue scheduler)
    └── cuda/calculate_regular_acceleration.cpp (CUDA path)

irregular_routines.cpp
    └── compute_acceleration.cpp (via queue scheduler)
    └── timestep_routines.cpp (via calculate_time_step_irr)
    └── update_particle.cpp (update_particle, calculate_time_step_irr/reg)
    └── FewBody/* (SDAR integration)

update_particle.cpp
    └── timestep_routines.cpp (getNewTimeStep*, getBlockTimeStep)
```

## 6. Recommended Order of Changes

Given dependencies and risk:

1. **Plan 1: timestep_routines.cpp** (Low risk)
   - Already takes raw arrays, minimal change
   - Update function signatures to take SoA pointers

2. **Plan 2: update_particle.cpp** (Medium risk)
   - Convert Particle methods to free functions
   - Add particle index parameter
   - Use ParticleData accessors

3. **Plan 3: compute_acceleration.cpp** (High risk, high reward)
   - Convert inner loops to SoA access
   - Maintain AoS interface for caller
   - Critical for energy conservation

4. **Plan 4: regular_routines.cpp** (Low risk)
   - Mainly orchestration
   - Update to pass ParticleData to subfunctions

5. **Plan 5: irregular_routines.cpp** (Medium risk)
   - Complex control flow
   - Update particle access patterns

## 7. Integration with GPU Code (Phase 3)

Phase 3 created `ParticleDataGPU` with SoA arrays. The connection point is in `cuda_acceleration.cu`:

- `_ReceiveFromHost()` currently uses AoS vectors (`i_particle_t`, `j_particle_t`)
- Added AoS→SoA bridge for kernel call
- Phase 4 should update callers to provide SoA data directly

Key files:
- `src/cuda/calculate_regular_acceleration.cpp` — Calls GPU functions
- `src/cuda/cuda_acceleration.cu` — GPU dispatch

## 8. Risk Areas for Energy Conservation

**Critical:** Energy conservation is the primary acceptance criterion.

Risk factors:
1. **Floating-point order changes** — SoA may change loop order, affecting FP accumulation
2. **Acceleration storage** — acc_total = acc_regular + acc_irregular must remain exact
3. **Time stepping** — Any change to timestep calculation affects integration accuracy
4. **Prediction/correction** — 4th order Hermite integration is sensitive

Mitigation:
- Keep identical mathematical operations
- Run baseline comparison after each plan
- Focus on data layout, not algorithm changes

## 9. Fields NOT in ParticleData

The following Particle fields are **not** in the ParticleData container:

- `background_acc[3]` — Background acceleration
- `binary_state` — Binary interrupt state + pair ID
- `spin_param[3]` — Spin parameters
- `group_info` — Group pointer
- `cm_particle_index` — CM parent index
- `members[10]` — Group member array
- `num_members` — Member count
- `new_members[10]` — New member array
- `new_num_members` — New member count
- `time_check` — Interrupt check time
- `radius` — Physical radius (in ParticleData)
- `delta_mass` — Mass delta (in ParticleData)

**Implication:** FewBody and SDAR code will continue to use Particle struct for group-related operations. Phase 5 (SDAR Compatibility) will address this with a proxy pattern.

## 10. Summary

| Requirement | Files | Complexity | Risk |
|-------------|-------|------------|------|
| CPU-01 (Regular) | compute_acceleration.cpp | High | High |
| CPU-02 (Irregular) | compute_acceleration.cpp | High | High |
| CPU-03 (Prediction) | update_particle.cpp, particle.h | Medium | Medium |
| CPU-04 (Correction) | update_particle.cpp | Medium | Medium |
| CPU-05 (Timestep) | update_particle.cpp, timestep_routines.cpp | Low | Low |

**Recommended approach:** Incremental conversion with bridge layers, similar to Phase 3's AoS→SoA bridge in cuda_acceleration.cu.

---
*Generated: 2026-01-17*
