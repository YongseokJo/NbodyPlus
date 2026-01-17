# Plan 04-02: Update Timestep Routines for SoA Compatibility

## Frontmatter

```yaml
wave: 1
depends_on: []
files_modified:
  - src/timestep_routines.cpp
  - src/Particle/update_particle.cpp
autonomous: true
```

## Objective

Update timestep calculation routines to accept SoA-style array pointers instead of relying on Particle member access. The functions in timestep_routines.cpp already take raw arrays, so changes are minimal.

## Context

Current functions in timestep_routines.cpp:
- `getNewTimeStepReg(double v[3], double df[3][4])` — Takes raw arrays
- `getNewTimeStepIrr(double f[3][4], double df[3][4])` — Takes raw arrays
- `getNewTimeStep(double f[3][4], double df[3][4])` — Takes raw arrays
- `getBlockTimeStep(double dt, int& TimeLevel, ull_t& TimeBlock, double& TimeStep)` — Pure calculation

These are called from Particle methods in update_particle.cpp which access `this->velocity`, `this->acc_total`, etc.

## Tasks

<task id="1">
Add overloaded versions of getNewTimeStepReg/Irr that accept ParticleData& and particle index
</task>

<task id="2">
Add free function versions of calculate_time_step_irr and calculate_time_step_reg that take ParticleData& and index
</task>

<task id="3">
Keep existing Particle:: method implementations as wrappers calling the new functions
</task>

<task id="4">
Add helper to extract velocity and acceleration as contiguous arrays from ParticleData
</task>

## Verification

- [ ] All existing tests pass (timestep calculations unchanged)
- [ ] New SoA-compatible functions produce identical results to existing methods
- [ ] No changes to timestep quantization logic

## must_haves

- Timestep calculations produce bit-identical results
- Block time step quantization unchanged
- Backward compatibility maintained via wrapper methods
