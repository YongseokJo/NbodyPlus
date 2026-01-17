# Plan 01: ParticleData Class Declaration

```yaml
wave: 1
depends_on: []
files_modified:
  - src/particle_data.h
autonomous: true
```

## Objective

Create the `ParticleData` SoA container class header with all array declarations and accessor functions.

## Tasks

<task id="1">
Create `src/particle_data.h` with:
- Include guards
- Necessary includes (cstddef for size_t, def.h for ull_t)
- Forward declarations if needed
</task>

<task id="2">
Declare private member arrays for critical fields (43 double arrays):
- Position: `pos_x_`, `pos_y_`, `pos_z_`
- Velocity: `vel_x_`, `vel_y_`, `vel_z_`
- Mass: `mass_`
- Accelerations: `acc_total_[3][4]`, `acc_reg_[3][4]`, `acc_irr_[3][4]` (array of pointers)
- New position/velocity: `new_pos_x_`, `new_pos_y_`, `new_pos_z_`, `new_vel_x_`, `new_vel_y_`, `new_vel_z_`
- `neighbor_radius_sq_`
</task>

<task id="3">
Declare private member arrays for important fields (13 arrays):
- Timestep doubles: `current_time_irr_`, `current_time_reg_`, `time_step_irr_`, `time_step_reg_`
- Block times (ull_t): `current_block_irr_`, `new_current_block_irr_`, `current_block_reg_`, `time_block_irr_`, `time_block_reg_`, `next_block_irr_`
- Neighbor ints: `num_neighbors_`, `new_num_neighbors_`, `neighbors_offset_`
</task>

<task id="4">
Declare private member arrays for low priority fields (10 arrays):
- Ints: `pid_`, `particle_index_`, `particle_type_`, `time_level_irr_`, `time_level_reg_`
- Bools: `is_active_`, `is_up_to_date_`, `is_cm_particle_`
- Doubles: `radius_`, `delta_mass_`
</task>

<task id="5">
Declare capacity and count members:
```cpp
size_t capacity_;
size_t count_;
```
</task>

<task id="6">
Declare public lifecycle methods:
```cpp
void allocate(size_t capacity);
void deallocate();
size_t size() const;
size_t capacity() const;
void set_size(size_t n);
```
</task>

<task id="7">
Declare accessor functions for scalar fields using pattern:
```cpp
// For each scalar field X:
double* X();
const double* X() const;
double get_X(size_t i) const;
void set_X(size_t i, double val);
```
Apply to: mass, neighbor_radius_sq, current_time_irr, current_time_reg, time_step_irr, time_step_reg, radius, delta_mass
</task>

<task id="8">
Declare accessor functions for 3D vector fields (position, velocity, new_position, new_velocity):
```cpp
// Component accessors
double* pos_x();
double* pos_y();
double* pos_z();
// ... const versions, get/set versions

// Convenience 3D accessors
void get_position(size_t i, double& x, double& y, double& z) const;
void set_position(size_t i, double x, double y, double z);
```
</task>

<task id="9">
Declare accessor functions for acceleration arrays:
```cpp
// Indexed access to acceleration array pointers
double* acc_total(int dim, int order);
const double* acc_total(int dim, int order) const;
double get_acc_total(size_t i, int dim, int order) const;
void set_acc_total(size_t i, int dim, int order, double val);

// Same pattern for acc_reg, acc_irr
```
</task>

<task id="10">
Declare accessor functions for ull_t fields:
```cpp
ull_t* current_block_irr();
ull_t get_current_block_irr(size_t i) const;
void set_current_block_irr(size_t i, ull_t val);
// ... same pattern for other ull_t fields
```
</task>

<task id="11">
Declare accessor functions for int fields:
```cpp
int* pid();
int get_pid(size_t i) const;
void set_pid(size_t i, int val);
// ... same pattern for other int fields
```
</task>

<task id="12">
Declare accessor functions for bool fields:
```cpp
bool* is_active();
bool get_is_active(size_t i) const;
void set_is_active(size_t i, bool val);
// ... same pattern for other bool fields
```
</task>

<task id="13">
Add constructor and destructor declarations:
```cpp
ParticleData();
~ParticleData();
// Delete copy constructor/assignment (manage memory manually)
ParticleData(const ParticleData&) = delete;
ParticleData& operator=(const ParticleData&) = delete;
```
</task>

## Verification

- [ ] File compiles with `g++ -std=c++11 -c src/particle_data.h -o /dev/null`
- [ ] All 66 array pointers declared
- [ ] All accessor functions declared for each array type
- [ ] Include guards present
- [ ] Copy semantics deleted (memory safety)

## must_haves

- [ ] ParticleData class declared with all 66 SoA array pointers
- [ ] Accessor functions declared for all field types (scalar, 3D vector, acceleration, ull_t, int, bool)
- [ ] Memory lifecycle methods declared (allocate, deallocate, size, capacity)
