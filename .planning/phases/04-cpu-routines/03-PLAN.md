# Plan 04-03: Update Prediction and Correction Routines

## Frontmatter

```yaml
wave: 2
depends_on: [01-PLAN.md]
files_modified:
  - src/Particle/update_particle.cpp
  - src/particle.h
autonomous: true
```

## Objective

Add SoA-compatible versions of prediction and correction routines. The key methods are:
- `predict_particle_second_order()` — Predicts position/velocity
- `correct_particle_fourth_order()` — 4th order correction
- `update_particle()` — Applies new_position/new_velocity
- `update_radius()` — Adjusts neighbor radius

## Context

These methods currently use `this->` to access Particle fields. We need:
1. Free function versions taking ParticleData& and index
2. Keep Particle methods as wrappers for backward compatibility

The prediction function is performance-critical as it's called for every particle pair in force calculations.

## Tasks

<task id="1">
Add predict_second_order() free function that operates on ParticleData
Signature: void predict_second_order(const ParticleData& data, size_t i, double dt, double pos_out[3], double vel_out[3])
</task>

<task id="2">
Add correct_fourth_order() free function
Signature: void correct_fourth_order(ParticleData& data, size_t i, double dt, const double pos[3], const double vel[3], const double a[3][4])
</task>

<task id="3">
Add update_particle_state() free function
Signature: void update_particle_state(ParticleData& data, size_t i)
Copies new_pos/vel to pos/vel
</task>

<task id="4">
Add update_neighbor_radius() free function
Signature: void update_neighbor_radius(ParticleData& data, size_t i, int fixed_num_neighbors)
</task>

<task id="5">
Keep existing Particle methods as wrappers
They can extract data to local arrays, call free functions, then write back
</task>

## Verification

- [ ] Free functions produce identical results to Particle methods
- [ ] Prediction output matches original for same input
- [ ] 4th order correction preserves energy conservation

## must_haves

- predict_second_order produces bit-identical position/velocity output
- correct_fourth_order applies exact same numerical corrections
- No change to integration order or coefficients
