# SoA Architecture for ABYSS

## Container Class Design

### Primary Container: ParticleData

```cpp
// src/particle_data.h
class ParticleData {
public:
    // Lifecycle
    void allocate(size_t capacity);
    void deallocate();
    size_t size() const { return count_; }
    size_t capacity() const { return capacity_; }

    // ========================================
    // Position accessors
    // ========================================
    double* pos_x() { return pos_x_; }
    double* pos_y() { return pos_y_; }
    double* pos_z() { return pos_z_; }

    void get_position(size_t i, double& x, double& y, double& z) const {
        x = pos_x_[i]; y = pos_y_[i]; z = pos_z_[i];
    }
    void set_position(size_t i, double x, double y, double z) {
        pos_x_[i] = x; pos_y_[i] = y; pos_z_[i] = z;
    }

    // ... similar for velocity, mass, accelerations, etc.

private:
    size_t capacity_;
    size_t count_;

    // Core kinematics
    double* pos_x_;
    double* pos_y_;
    double* pos_z_;
    double* vel_x_;
    double* vel_y_;
    double* vel_z_;
    double* mass_;

    // Acceleration arrays (3 components × 4 derivatives × 3 types)
    double* acc_total_x_[4];   // [derivative][particle]
    double* acc_total_y_[4];
    double* acc_total_z_[4];
    double* acc_regular_x_[4];
    double* acc_regular_y_[4];
    double* acc_regular_z_[4];
    double* acc_irregular_x_[4];
    double* acc_irregular_y_[4];
    double* acc_irregular_z_[4];

    // Timestep data
    double* current_time_irr_;
    double* current_time_reg_;
    // ... etc.
};
```

## Memory Allocation Strategy

### CPU/MPI Shared Memory

```cpp
void ParticleData::allocate_shared(MPI_Comm shared_comm, size_t capacity) {
    capacity_ = capacity;

    // Allocate each array as separate MPI window
    MPI_Win_allocate_shared(capacity * sizeof(double), sizeof(double),
        MPI_INFO_NULL, shared_comm, &pos_x_, &win_pos_x_);
    MPI_Win_allocate_shared(capacity * sizeof(double), sizeof(double),
        MPI_INFO_NULL, shared_comm, &pos_y_, &win_pos_y_);
    // ... repeat for all arrays
}
```

### GPU Memory

```cpp
// Separate device arrays
struct ParticleDataGPU {
    cuda_real_t* d_pos_x;
    cuda_real_t* d_pos_y;
    cuda_real_t* d_pos_z;
    cuda_real_t* d_vel_x;
    cuda_real_t* d_vel_y;
    cuda_real_t* d_vel_z;
    cuda_real_t* d_mass;
    // ... etc.

    void allocate(size_t n);
    void copy_from_host(const ParticleData& host, size_t n);
    void copy_to_host(ParticleData& host, size_t n);
};
```

## Data Flow Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    ParticleData (SoA)                    │
│  MPI Shared Memory: multiple MPI_Win objects            │
├─────────────────────────────────────────────────────────┤
│  pos_x[N], pos_y[N], pos_z[N]                          │
│  vel_x[N], vel_y[N], vel_z[N]                          │
│  mass[N]                                                │
│  acc_total_x[4][N], acc_total_y[4][N], ...             │
│  ... other arrays ...                                   │
└─────────────────────────────────────────────────────────┘
           │                              │
           │ Direct pointer access        │ cudaMemcpy per array
           ▼                              ▼
┌─────────────────────┐    ┌─────────────────────────────┐
│   CPU Routines      │    │      GPU (ParticleDataGPU)  │
│   - Irregular       │    │  d_pos_x[N], d_pos_y[N]...  │
│   - SDAR (via proxy)│    │                             │
│   - I/O             │    │  CUDA kernels access        │
└─────────────────────┘    │  arrays directly            │
                           └─────────────────────────────┘
```

## SDAR Compatibility Layer

SDAR expects a `Particle` object with fields. Create a lightweight proxy:

```cpp
// Proxy for SDAR compatibility
class ParticleProxy {
public:
    ParticleProxy(ParticleData& data, size_t index)
        : data_(data), i_(index) {}

    // SDAR-compatible accessors
    double* getPos() { return pos_buf_; }
    double* getVel() { return vel_buf_; }
    double& mass { return data_.mass()[i_]; }

    // Sync buffer to SoA (call after SDAR modifies)
    void sync_to_soa() {
        data_.set_position(i_, pos_buf_[0], pos_buf_[1], pos_buf_[2]);
        data_.set_velocity(i_, vel_buf_[0], vel_buf_[1], vel_buf_[2]);
    }

    // Sync SoA to buffer (call before SDAR reads)
    void sync_from_soa() {
        data_.get_position(i_, pos_buf_[0], pos_buf_[1], pos_buf_[2]);
        data_.get_velocity(i_, vel_buf_[0], vel_buf_[1], vel_buf_[2]);
    }

private:
    ParticleData& data_;
    size_t i_;
    double pos_buf_[3];  // Temporary buffer for SDAR
    double vel_buf_[3];
};
```

## File Organization

```
src/
├── particle_data.h          # SoA container declaration
├── particle_data.cpp        # SoA container implementation
├── particle_data_gpu.h      # GPU SoA container
├── particle_data_gpu.cu     # GPU allocation/transfer
├── particle_proxy.h         # SDAR compatibility proxy
├── particle.h               # Keep for SDAR (minimal)
└── cuda/
    ├── cuda_defs.h          # GPU types (SoA pointers)
    └── cuda_kernels.cu      # Updated for SoA access
```

## Kernel Access Pattern

### Before (AoS)
```cpp
__global__ void compute_forces(j_particle_t* particles, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float px = particles[i].pos_x;  // Non-coalesced
    float py = particles[i].pos_y;
    float pz = particles[i].pos_z;
}
```

### After (SoA)
```cpp
__global__ void compute_forces(
    float* pos_x, float* pos_y, float* pos_z,
    float* mass, int N
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float px = pos_x[i];  // Coalesced!
    float py = pos_y[i];
    float pz = pos_z[i];
}
```

## Build Order

Suggested implementation sequence based on dependencies:

1. **ParticleData container** — core SoA class
2. **MPI allocation** — shared memory windows
3. **GPU container** — device arrays
4. **Accessor integration** — update access patterns
5. **CUDA kernels** — convert to SoA access
6. **CPU routines** — update force calculations
7. **SDAR proxy** — compatibility layer
8. **I/O** — HDF5 output from SoA
