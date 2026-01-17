# Plan 02: ParticleData Implementation

```yaml
wave: 2
depends_on: [01]
files_modified:
  - src/particle_data.cpp
autonomous: true
```

## Objective

Implement the `ParticleData` class methods: memory allocation/deallocation with cache-friendly padding, and all accessor function bodies.

## Tasks

<task id="1">
Create `src/particle_data.cpp` with:
- Include `particle_data.h`
- Include `<cstring>` for memset if needed
</task>

<task id="2">
Implement helper function for padded capacity:
```cpp
namespace {
size_t padded_capacity(size_t n) {
    // Avoid exact powers of 2 (cache associativity issues)
    if (n > 64 && (n & (n - 1)) == 0) {
        return n + 64;
    }
    // Align to 8 elements (64 bytes for doubles)
    return ((n + 7) / 8) * 8;
}
}
```
</task>

<task id="3">
Implement constructor - initialize all pointers to nullptr:
```cpp
ParticleData::ParticleData()
    : capacity_(0), count_(0),
      pos_x_(nullptr), pos_y_(nullptr), pos_z_(nullptr),
      // ... all other pointers
{
    // Initialize acceleration pointer arrays
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            acc_total_[d][o] = nullptr;
            acc_reg_[d][o] = nullptr;
            acc_irr_[d][o] = nullptr;
        }
    }
}
```
</task>

<task id="4">
Implement destructor - call deallocate:
```cpp
ParticleData::~ParticleData() {
    deallocate();
}
```
</task>

<task id="5">
Implement `allocate(size_t capacity)`:
- Calculate padded capacity
- Allocate all 66 arrays with `new`
- Initialize all arrays to zero
- Store capacity and set count to 0

Structure:
```cpp
void ParticleData::allocate(size_t requested_capacity) {
    if (capacity_ > 0) {
        deallocate();  // Clean up existing allocation
    }

    capacity_ = padded_capacity(requested_capacity);
    count_ = 0;

    // Critical double arrays (43)
    pos_x_ = new double[capacity_]();
    pos_y_ = new double[capacity_]();
    pos_z_ = new double[capacity_]();
    // ... velocity, mass, accelerations, etc.

    // Acceleration arrays (36 total)
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            acc_total_[d][o] = new double[capacity_]();
            acc_reg_[d][o] = new double[capacity_]();
            acc_irr_[d][o] = new double[capacity_]();
        }
    }

    // Important ull_t arrays (5)
    current_block_irr_ = new ull_t[capacity_]();
    // ...

    // Important int arrays (3)
    num_neighbors_ = new int[capacity_]();
    // ...

    // Low priority arrays
    // ...
}
```
</task>

<task id="6">
Implement `deallocate()`:
```cpp
void ParticleData::deallocate() {
    if (capacity_ == 0) return;

    // Delete all arrays
    delete[] pos_x_; pos_x_ = nullptr;
    delete[] pos_y_; pos_y_ = nullptr;
    // ... all other arrays

    // Acceleration arrays
    for (int d = 0; d < 3; d++) {
        for (int o = 0; o < 4; o++) {
            delete[] acc_total_[d][o]; acc_total_[d][o] = nullptr;
            delete[] acc_reg_[d][o]; acc_reg_[d][o] = nullptr;
            delete[] acc_irr_[d][o]; acc_irr_[d][o] = nullptr;
        }
    }

    capacity_ = 0;
    count_ = 0;
}
```
</task>

<task id="7">
Implement size/capacity methods:
```cpp
size_t ParticleData::size() const { return count_; }
size_t ParticleData::capacity() const { return capacity_; }
void ParticleData::set_size(size_t n) { count_ = n; }
```
</task>

<task id="8">
Implement inline accessor functions (can be in header or .cpp):

For scalar double fields:
```cpp
double* ParticleData::mass() { return mass_; }
const double* ParticleData::mass() const { return mass_; }
double ParticleData::get_mass(size_t i) const { return mass_[i]; }
void ParticleData::set_mass(size_t i, double val) { mass_[i] = val; }
```
Repeat for: neighbor_radius_sq, current_time_irr, current_time_reg, time_step_irr, time_step_reg, radius, delta_mass
</task>

<task id="9">
Implement 3D vector accessors for position:
```cpp
double* ParticleData::pos_x() { return pos_x_; }
double* ParticleData::pos_y() { return pos_y_; }
double* ParticleData::pos_z() { return pos_z_; }
// const versions...
double ParticleData::get_pos_x(size_t i) const { return pos_x_[i]; }
// set versions...

void ParticleData::get_position(size_t i, double& x, double& y, double& z) const {
    x = pos_x_[i];
    y = pos_y_[i];
    z = pos_z_[i];
}

void ParticleData::set_position(size_t i, double x, double y, double z) {
    pos_x_[i] = x;
    pos_y_[i] = y;
    pos_z_[i] = z;
}
```
Repeat pattern for: velocity, new_position, new_velocity
</task>

<task id="10">
Implement acceleration array accessors:
```cpp
double* ParticleData::acc_total(int dim, int order) {
    return acc_total_[dim][order];
}
const double* ParticleData::acc_total(int dim, int order) const {
    return acc_total_[dim][order];
}
double ParticleData::get_acc_total(size_t i, int dim, int order) const {
    return acc_total_[dim][order][i];
}
void ParticleData::set_acc_total(size_t i, int dim, int order, double val) {
    acc_total_[dim][order][i] = val;
}
```
Repeat for: acc_reg, acc_irr
</task>

<task id="11">
Implement ull_t field accessors:
```cpp
ull_t* ParticleData::current_block_irr() { return current_block_irr_; }
const ull_t* ParticleData::current_block_irr() const { return current_block_irr_; }
ull_t ParticleData::get_current_block_irr(size_t i) const { return current_block_irr_[i]; }
void ParticleData::set_current_block_irr(size_t i, ull_t val) { current_block_irr_[i] = val; }
```
Repeat for: new_current_block_irr, current_block_reg, time_block_irr, time_block_reg, next_block_irr
</task>

<task id="12">
Implement int field accessors:
```cpp
int* ParticleData::num_neighbors() { return num_neighbors_; }
int ParticleData::get_num_neighbors(size_t i) const { return num_neighbors_[i]; }
void ParticleData::set_num_neighbors(size_t i, int val) { num_neighbors_[i] = val; }
```
Repeat for: new_num_neighbors, neighbors_offset, pid, particle_index, particle_type, time_level_irr, time_level_reg
</task>

<task id="13">
Implement bool field accessors:
```cpp
bool* ParticleData::is_active() { return is_active_; }
bool ParticleData::get_is_active(size_t i) const { return is_active_[i]; }
void ParticleData::set_is_active(size_t i, bool val) { is_active_[i] = val; }
```
Repeat for: is_up_to_date, is_cm_particle
</task>

## Verification

- [ ] File compiles with `g++ -std=c++11 -c src/particle_data.cpp -I src/`
- [ ] `allocate()` allocates all 66 arrays
- [ ] `deallocate()` frees all memory and nulls pointers
- [ ] Padded capacity avoids powers of 2
- [ ] All accessor functions implemented

## must_haves

- [ ] Memory allocation with cache-friendly padding implemented
- [ ] All 66 arrays allocated in `allocate()` and freed in `deallocate()`
- [ ] All accessor function bodies implemented
- [ ] Constructor initializes all pointers to nullptr
- [ ] Destructor calls deallocate()
