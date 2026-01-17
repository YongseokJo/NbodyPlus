# SoA Conversion Scope

## Field Classification

Based on access patterns in ABYSS, fields are classified by performance impact.

### Critical (Must Convert) — Hot Path Fields

These are accessed in O(N²) force loops and CUDA kernels:

| Field | Current Type | SoA Type | Rationale |
|-------|--------------|----------|-----------|
| `position[3]` | double[3] | 3× double* | Force calculation inner loop |
| `velocity[3]` | double[3] | 3× double* | Prediction, correction |
| `mass` | double | double* | Force calculation |
| `acc_total[3][4]` | double[3][4] | 12× double* | Hermite integration |
| `acc_regular[3][4]` | double[3][4] | 12× double* | Regular force update |
| `acc_irregular[3][4]` | double[3][4] | 12× double* | Irregular force update |
| `new_position[3]` | double[3] | 3× double* | Prediction output |
| `new_velocity[3]` | double[3] | 3× double* | Prediction output |
| `neighbor_radius_sq` | double | double* | Neighbor finding |

**Total critical arrays: ~42 double arrays**

### Important (Should Convert) — Warm Path Fields

Accessed per-timestep for subset of particles:

| Field | Current Type | SoA Type | Rationale |
|-------|--------------|----------|-----------|
| `current_time_irr` | double | double* | Timestep selection |
| `current_time_reg` | double | double* | Timestep selection |
| `time_step_irr` | double | double* | Timestep calculation |
| `time_step_reg` | double | double* | Timestep calculation |
| `current_block_irr` | ull_t | ull_t* | Block time management |
| `current_block_reg` | ull_t | ull_t* | Block time management |
| `time_block_irr` | ull_t | ull_t* | Block time |
| `time_block_reg` | ull_t | ull_t* | Block time |
| `next_block_irr` | ull_t | ull_t* | Scheduling |
| `num_neighbors` | int | int* | Neighbor list |
| `new_num_neighbors` | int | int* | Neighbor update |
| `neighbors_offset` | int | int* | Neighbor array index |

**Total important arrays: ~12 arrays**

### Low Priority (Can Convert) — Cold Path Fields

Accessed occasionally or per-particle during specific events:

| Field | Current Type | Notes |
|-------|--------------|-------|
| `pid` | int | Lookup only |
| `particle_index` | int | Identity |
| `particle_type` | int | SEVN classification |
| `is_active` | bool | Sparse access |
| `is_up_to_date` | bool | Sync flag |
| `is_cm_particle` | bool | Few-body only |
| `time_level_irr` | int | Level tracking |
| `time_level_reg` | int | Level tracking |

### Keep As-Is — Pointer/Object Fields

Don't convert to SoA (remain as separate per-particle data):

| Field | Type | Rationale |
|-------|------|-----------|
| `group_info` | Group* | Pointer to SDAR group |
| `stellar_evolution` | StarSEVN* | SEVN object |
| `binary_evolution` | Binstar* | SEVN binary |
| `members[10]` | int[10] | Small fixed array |
| `new_members[10]` | int[10] | Small fixed array |
| `spin_param[3]` | double[3] | Rarely accessed |
| `background_acc[3]` | double[3] | External force |
| `binary_state` | long long | State tracking |

## GPU Types Conversion

### j_particle_t (source particles in force calc)

Current:
```cpp
struct j_particle_t {
    cuda_real_t pos_x, pos_y, pos_z;
    cuda_real_t vel_x, vel_y, vel_z;
    cuda_real_t mass;
    int index;
};
```

SoA:
```cpp
struct JParticleArrays {
    cuda_real_t* pos_x;  // [N]
    cuda_real_t* pos_y;  // [N]
    cuda_real_t* pos_z;  // [N]
    cuda_real_t* vel_x;  // [N]
    cuda_real_t* vel_y;  // [N]
    cuda_real_t* vel_z;  // [N]
    cuda_real_t* mass;   // [N]
    int* index;          // [N]
};
```

### i_particle_t (target particles in force calc)

Similar conversion pattern.

## Accessor Function Requirements

For each SoA array, provide:

1. **Direct array pointer** — for CUDA memcpy, MPI operations
2. **Single element access** — `get_pos_x(i)`, `set_pos_x(i, val)`
3. **Bulk access** — for vectorized CPU loops

Example:
```cpp
// Direct access for GPU
double* pos_x() { return pos_x_; }
const double* pos_x() const { return pos_x_; }

// Single element
double get_pos_x(size_t i) const { return pos_x_[i]; }
void set_pos_x(size_t i, double val) { pos_x_[i] = val; }

// 3D position helpers
void get_position(size_t i, double& x, double& y, double& z) const;
void set_position(size_t i, double x, double y, double z);
```

## Conversion Scope Summary

| Category | Array Count | Priority |
|----------|-------------|----------|
| Critical (hot path) | ~42 | Phase 1 |
| Important (warm path) | ~12 | Phase 2 |
| Low priority (cold) | ~8 | Phase 3 |
| Keep as-is (pointers) | ~10 | Skip |

**Total new arrays: ~62 for full conversion**
