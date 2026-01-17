# Phase 1 Research: Core SoA Container

## Field Inventory from Particle Struct

Analyzed `src/particle.h`. Constants: `DIM=3`, `HERMITE_ORDER=4`.

### Critical Fields (Hot Path) — 43 arrays

| Field | Type | SoA Arrays | Notes |
|-------|------|------------|-------|
| `position[3]` | double[3] | `pos_x`, `pos_y`, `pos_z` | Force calc inner loop |
| `velocity[3]` | double[3] | `vel_x`, `vel_y`, `vel_z` | Prediction/correction |
| `mass` | double | `mass` | Force calculation |
| `acc_total[3][4]` | double[3][4] | 12 arrays: `acc_total_x[0..3]`, `_y`, `_z` | Hermite integration |
| `acc_regular[3][4]` | double[3][4] | 12 arrays: `acc_reg_x[0..3]`, `_y`, `_z` | Regular force |
| `acc_irregular[3][4]` | double[3][4] | 12 arrays: `acc_irr_x[0..3]`, `_y`, `_z` | Irregular force |
| `new_position[3]` | double[3] | `new_pos_x`, `new_pos_y`, `new_pos_z` | Prediction output |
| `new_velocity[3]` | double[3] | `new_vel_x`, `new_vel_y`, `new_vel_z` | Prediction output |
| `neighbor_radius_sq` | double | `neighbor_radius_sq` | Neighbor finding |

**Total: 43 double arrays**

### Important Fields (Warm Path) — 13 arrays

| Field | Type | SoA Array | Notes |
|-------|------|-----------|-------|
| `current_time_irr` | double | `current_time_irr` | Timestep selection |
| `current_time_reg` | double | `current_time_reg` | Timestep selection |
| `time_step_irr` | double | `time_step_irr` | Timestep calculation |
| `time_step_reg` | double | `time_step_reg` | Timestep calculation |
| `current_block_irr` | ull_t | `current_block_irr` | Block time |
| `new_current_block_irr` | ull_t | `new_current_block_irr` | New block (found in struct!) |
| `current_block_reg` | ull_t | `current_block_reg` | Block time |
| `time_block_irr` | ull_t | `time_block_irr` | Block time |
| `time_block_reg` | ull_t | `time_block_reg` | Block time |
| `next_block_irr` | ull_t | `next_block_irr` | Scheduling |
| `num_neighbors` | int | `num_neighbors` | Neighbor count |
| `new_num_neighbors` | int | `new_num_neighbors` | New neighbor count |
| `neighbors_offset` | int | `neighbors_offset` | Neighbor array index |

**Total: 13 arrays (5 double, 5 ull_t, 3 int)**

### Low Priority Fields — 10 arrays

| Field | Type | SoA Array | Notes |
|-------|------|-----------|-------|
| `pid` | int | `pid` | Particle ID |
| `particle_index` | int | `particle_index` | Array index |
| `particle_type` | int | `particle_type` | Star type |
| `is_active` | bool | `is_active` | Active flag |
| `is_up_to_date` | bool | `is_up_to_date` | Sync flag |
| `is_cm_particle` | bool | `is_cm_particle` | CM flag |
| `time_level_irr` | int | `time_level_irr` | Level tracking |
| `time_level_reg` | int | `time_level_reg` | Level tracking |
| `radius` | double | `radius` | Physical radius |
| `delta_mass` | double | `delta_mass` | Mass distribution |

**Total: 10 arrays (5 int, 3 bool, 2 double)**

### Keep As-Is (Don't Convert)

| Field | Type | Rationale |
|-------|------|-----------|
| `group_info` | Group* | Pointer to SDAR group |
| `cm_particle_index` | int | CM reference (sparse) |
| `num_members` | int | CM only |
| `members[10]` | int[10] | Small fixed array, CM only |
| `new_num_members` | int | CM only |
| `new_members[10]` | int[10] | Small fixed array, CM only |
| `binary_state` | long long | State tracking |
| `time_check` | double | Interrupt check |
| `spin_param[3]` | double[3] | Rarely accessed |
| `background_acc[3]` | double[3] | External force |
| `stellar_evolution` | StarSEVN* | SEVN object (ifdef SEVN) |
| `binary_evolution` | Binstar* | SEVN binary (ifdef SEVN) |
| `formation_time` | double | SEVN only |
| `world_time` | double | SEVN only |

---

## Complete Array Count

| Category | Type | Count |
|----------|------|-------|
| Critical | double | 43 |
| Important | double | 5 |
| Important | ull_t | 5 |
| Important | int | 3 |
| Low Priority | int | 5 |
| Low Priority | bool | 3 |
| Low Priority | double | 2 |
| **Total** | — | **66** |

---

## Accessor Patterns Needed

### Pattern 1: Scalar fields (mass, radius, etc.)
```cpp
double* mass() { return mass_; }
const double* mass() const { return mass_; }
double get_mass(size_t i) const { return mass_[i]; }
void set_mass(size_t i, double val) { mass_[i] = val; }
```

### Pattern 2: 3D vector fields (position, velocity)
```cpp
// Component arrays
double* pos_x() { return pos_x_; }
double* pos_y() { return pos_y_; }
double* pos_z() { return pos_z_; }

// Component access
double get_pos_x(size_t i) const { return pos_x_[i]; }
void set_pos_x(size_t i, double val) { pos_x_[i] = val; }

// Convenience: 3D access
void get_position(size_t i, double& x, double& y, double& z) const {
    x = pos_x_[i]; y = pos_y_[i]; z = pos_z_[i];
}
void set_position(size_t i, double x, double y, double z) {
    pos_x_[i] = x; pos_y_[i] = y; pos_z_[i] = z;
}
```

### Pattern 3: Acceleration arrays [dim][order]
```cpp
// Naming: acc_{type}_{dim}_{order}  e.g. acc_total_x_0, acc_total_x_1
// Or: arrays of pointers for indexed access

// Option A: Named arrays (verbose but clear)
double* acc_total_x_0() { return acc_total_x_0_; }
double* acc_total_x_1() { return acc_total_x_1_; }
// ... 36 more

// Option B: Array of pointers (compact)
double* acc_total_[3][4];  // [dim][order]
double* acc_total(int dim, int order) { return acc_total_[dim][order]; }
double get_acc_total(size_t i, int dim, int order) const {
    return acc_total_[dim][order][i];
}
void set_acc_total(size_t i, int dim, int order, double val) {
    acc_total_[dim][order][i] = val;
}
```

**Recommendation**: Use Option B (array of pointers) for accelerations to avoid 36 separate named arrays and match the original `acc_total[dim][order]` access pattern.

---

## Memory Allocation Strategy

### Padding for Cache Conflict Avoidance

```cpp
size_t padded_capacity(size_t n) {
    // Avoid exact powers of 2 (cache associativity issues)
    if (n > 64 && (n & (n - 1)) == 0) {
        return n + 64;  // Add padding
    }
    // Align to 64 bytes (cache line)
    return ((n + 7) / 8) * 8;
}
```

### Allocation Order

Allocate arrays in groups for locality:
1. Position arrays (pos_x, pos_y, pos_z)
2. Velocity arrays (vel_x, vel_y, vel_z)
3. Mass
4. Acceleration arrays (36 total)
5. New position/velocity (6 arrays)
6. Timestep data
7. Neighbor data
8. Flags/indices

### Single Allocation Block (Alternative)

Could allocate one large block and partition:
```cpp
// Total bytes needed
size_t total = capacity * (
    43 * sizeof(double) +    // critical doubles
    5 * sizeof(double) +     // important doubles
    5 * sizeof(ull_t) +      // important ull_t
    3 * sizeof(int) +        // important ints
    // ... etc
);
char* block = new char[total];
// Partition into arrays...
```

**Recommendation**: Separate allocations per array for simplicity and MPI compatibility (Phase 2 needs separate MPI_Win per array).

---

## Phase 1 Scope Decision

**For Phase 1 (Core Container):**
- Include all 66 arrays (critical + important + low priority)
- Exclude pointer/object fields (group_info, stellar_evolution, etc.)
- Use standard `new[]`/`delete[]` allocation (MPI shared memory comes in Phase 2)

**Rationale**: Build complete container now, swap allocation method in Phase 2. Avoids partial container that needs expansion later.

---

## Files to Create

1. `src/particle_data.h` — Class declaration
   - Private array pointers
   - Public accessor functions
   - allocate()/deallocate() methods
   - size()/capacity() methods

2. `src/particle_data.cpp` — Implementation
   - Memory allocation with padding
   - Deallocation
   - (Optional) copy/move semantics

---

## Key Implementation Notes

1. **Array naming convention**: `{field}_{dim}_{order}_` for internal, `{field}_{dim}_{order}()` for accessor
2. **No runtime AoS↔SoA conversion**: Store permanently in SoA
3. **Const correctness**: Provide both const and non-const accessors
4. **Type aliases**: Use `ull_t` consistently (defined in def.h as `unsigned long long`)
5. **Include guards**: Standard `#ifndef PARTICLE_DATA_H`

---

## Questions Resolved

1. **How many arrays?** 66 total (43 critical + 13 important + 10 low priority)
2. **Acceleration layout?** Array of pointers `double* acc_total_[3][4]` for indexed access
3. **Allocation strategy?** Separate arrays, padded capacity, standard new/delete for Phase 1
4. **What to exclude?** Pointers (Group*, StarSEVN*), fixed small arrays (members[10]), SEVN-only fields

---
*Generated: 2026-01-17*
