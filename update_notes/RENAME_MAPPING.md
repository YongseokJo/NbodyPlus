# Variable/Function Rename Mapping (PascalCase/camelCase → snake_case)

## Refactoring Status

| File | Status | Notes |
|------|--------|-------|
| `src/def.h` | DONE | Constants, types with legacy aliases |
| `src/particle.h` | DONE | All members/methods with legacy aliases |
| `src/GlobalVariable.h` | DONE | Struct members with legacy aliases |
| `src/global.h` | DONE | All extern declarations with legacy aliases |
| `src/DefaultGlobal.cpp` | DONE | Variable definitions |
| `src/cuda/cuda_defs.h` | DONE | GPU structures with legacy aliases |
| `src/cuda/cuda_kernels.h` | DONE | Kernel declarations |
| `src/cuda/cuda_kernels.cu` | PARTIAL | Headers updated, kernels use macros |
| `src/Worker.h` | DONE | Worker struct with legacy aliases |
| `src/Queue.h` | DONE | TaskName enum with legacy aliases |
| Other .cpp files | PENDING | Should work via legacy macros |

---

## Constants (def.h)

| Old Name | New Name |
|----------|----------|
| MaxNumParticle | max_num_particle |
| MaxNumNeighbor | max_num_neighbor |
| MaxNeighborRadius | max_neighbor_radius |
| BatchSize | batch_size |
| GridDimY | grid_dim_y |
| NNB_per_block | nnb_per_block |
| Dim | DIM |
| HERMITE_ORDER | HERMITE_ORDER |
| MIN_LEVEL_BUFFER | MIN_LEVEL_BUFFER |
| CUDA_REAL | cuda_real_t |
| nbodymax | nbody_max |

## Particle Struct Members (particle.h)

| Old Name | New Name |
|----------|----------|
| PID | pid |
| ParticleIndex | particle_index |
| ParticleType | particle_type |
| Position | position |
| Velocity | velocity |
| Mass | mass |
| a_tot | acc_total |
| a_reg | acc_regular |
| a_irr | acc_irregular |
| NumberOfNeighbor | num_neighbors |
| NewNumberOfNeighbor | new_num_neighbors |
| NeighborsOffset | neighbors_offset |
| CurrentTimeIrr | current_time_irr |
| CurrentTimeReg | current_time_reg |
| CurrentBlockIrr | current_block_irr |
| NewCurrentBlockIrr | new_current_block_irr |
| CurrentBlockReg | current_block_reg |
| NextBlockIrr | next_block_irr |
| TimeStepIrr | time_step_irr |
| TimeStepReg | time_step_reg |
| TimeBlockIrr | time_block_irr |
| TimeBlockReg | time_block_reg |
| TimeLevelIrr | time_level_irr |
| TimeLevelReg | time_level_reg |
| NewPosition | new_position |
| NewVelocity | new_velocity |
| RadiusOfNeighbor | neighbor_radius_sq |
| BackgroundAcceleration | background_acc |
| isActive | is_active |
| isUpdateToDate | is_up_to_date |
| radius | radius |
| dm | delta_mass |
| time_check | time_check |
| binary_state | binary_state |
| a_spin | spin_param |
| GroupInfo | group_info |
| isCMptcl | is_cm_particle |
| CMPtclIndex | cm_particle_index |
| NumberOfMember | num_members |
| Members | members |
| NewNumberOfMember | new_num_members |
| NewMembers | new_members |
| StellarEvolution | stellar_evolution |
| BinaryEvolution | binary_evolution |
| FormationTime | formation_time |
| WorldTime | world_time |

## Particle Methods (particle.h)

| Old Name | New Name |
|----------|----------|
| initialize | initialize |
| clear | clear |
| normalizeParticle | normalize_particle |
| updateParticle | update_particle |
| predictParticleSecondOrder | predict_particle_second_order |
| correctParticleFourthOrder | correct_particle_fourth_order |
| updateRadius | update_radius |
| initializeTimeStep | initialize_time_step |
| computeAccelerationIrr | compute_acceleration_irr |
| computeAccelerationReg | compute_acceleration_reg |
| calculateTimeStepIrr | calculate_time_step_irr |
| calculateTimeStepIrr2 | calculate_time_step_irr_v2 |
| calculateTimeStepReg | calculate_time_step_reg |
| updateRegularParticleCuda | update_regular_particle_cuda |
| checkNewGroup | check_new_group |
| checkNewGroup2 | check_new_group_v2 |
| checkNewGroup3 | check_new_group_v3 |
| checkNewGroup4 | check_new_group_v4 |
| setBinaryPairID | set_binary_pair_id |
| setBinaryInterruptState | set_binary_interrupt_state |
| getBinaryInterruptState | get_binary_interrupt_state |
| getBinaryPairID | get_binary_pair_id |
| getPos | get_position |
| getVel | get_velocity |
| printColumnTitle | print_column_title |
| printColumn | print_column |
| copyNewMembers | copy_new_members |
| printParticleInfo | print_particle_info |

## Global Variables (global.h)

| Old Name | New Name |
|----------|----------|
| MyRank | my_rank |
| NumberOfProcessor | num_processors |
| NumberOfWorker | num_workers |
| shared_comm | shared_comm |
| particles | particles |
| global_variable | g_state |
| Neighbors | neighbors |
| Neighbors_original | neighbors_original |
| NewNeighbors | new_neighbors |
| NewNeighbors_original | new_neighbors_original |
| QueueType | queue_type_mpi |
| IparticleType | iparticle_type_mpi |
| JparticleType | jparticle_type_mpi |
| LastParticleIndex | last_particle_index |
| NumberOfParticle | num_particles |
| NewCMPID | new_cm_pid |
| FixNumNeighbor | fixed_num_neighbors |
| InitialNeighborRadius | initial_neighbor_radius |
| RSearch | r_search |
| TSearch | t_search |
| UseCompression | use_compression |
| CompressionLevel | compression_level |
| RestartEnabled | restart_enabled |
| CheckpointFile | checkpoint_file |
| CMPtclWorker | cm_particle_worker_map |
| PrevCMPtclWorker | prev_cm_particle_worker_map |
| global_time | global_time |
| global_time_irr | global_time_irr |
| NextRegTimeBlock | next_reg_time_block |
| time_block | time_block |
| time_step | time_step |
| block_max | block_max |
| endTime | end_time |
| EnzoTimeStep | enzo_time_step |
| outputTime | output_time |
| outNum | output_num |
| outputTimeStep | output_time_step |
| E_binary | energy_binary |
| E_binary_SD | energy_binary_sd |
| E_merger | energy_merger |
| E_PN | energy_pn |
| binout | bin_output_file |
| mergerout | merger_output_file |
| SEVNout | sevn_output_file |
| workerout | worker_output_file |

## CUDA Structures (cuda_defs.h)

| Old Name | New Name |
|----------|----------|
| Jparticle | j_particle_t |
| Iparticle | i_particle_t |
| posx, posy, posz | pos_x, pos_y, pos_z |
| velx, vely, velz | vel_x, vel_y, vel_z |
| mass | mass |
| index | index |
| r2 | radius_sq |
| dtr | dt_regular |

## CUDA Functions (cuda_kernels.cu)

| Old Name | New Name |
|----------|----------|
| compute_forces | compute_forces_kernel |
| reduce_forces_kernel | reduce_forces_kernel |
| gather_neighbor | gather_neighbor_kernel |
| gather_numneighbor | gather_num_neighbor_kernel |

## Core Routine Functions

| Old Name | New Name |
|----------|----------|
| IrregularRoutines | irregular_routines |
| RegularRoutines | regular_routines |
| createSkipList | create_skip_list |
| updateSkipList | update_skip_list |
| formBinaries | form_binaries |
| FBTermination | fb_termination |
| calculateRegAccelerationOnGPU | calculate_reg_acceleration_gpu |
| sendAllParticlesToGPU | send_all_particles_to_gpu |
| sendAllParticlesToGPU_Worker | send_all_particles_to_gpu_worker |
| getRegularList | get_regular_list |
| updateRegularMap | update_regular_map |
| GetAcceleration | get_acceleration |
| InitializeDevice | initialize_device |
| OpenDevice | open_device |
| CloseDevice | close_device |
| SendToDevice | send_to_device |
| ProfileDevice | profile_device |
| CalculateAccelerationOnDevice | calculate_acceleration_on_device |

## Queue/Task Names

| Old Name | New Name |
|----------|----------|
| TaskName | task_name_t |
| IrrForce | IRR_FORCE |
| IrrUpdate | IRR_UPDATE |
| RegForce | REG_FORCE |
| RegUpdate | REG_UPDATE |
| RegCuda | REG_CUDA |
| ARIntegration | AR_INTEGRATION |
| SearchGroup | SEARCH_GROUP |
| MakeGroup | MAKE_GROUP |
| DeleteGroup | DELETE_GROUP |
| MergeManyBody | MERGE_MANY_BODY |
| PrepareGPUCalc | PREPARE_GPU_CALC |
| Ends | TASK_END |

## Class Names

| Old Name | New Name |
|----------|----------|
| QueueScheduler | queue_scheduler_t |
| Worker | worker_t |
| SkipList | skip_list_t |
| Node | skip_list_node_t |
| GlobalVariable | global_state_t |
| Performance | performance_t |

## Enum Names

| Old Name | New Name |
|----------|----------|
| BinaryInterruptState | binary_interrupt_state_t |
| TimerID | timer_id_t |

---

## Files to Modify (in order)

1. `src/def.h` - Constants and type definitions
2. `src/particle.h` - Particle struct
3. `src/GlobalVariable.h` - Global state struct
4. `src/global.h` - Global variable declarations
5. `src/cuda/cuda_defs.h` - CUDA structures
6. `src/cuda/cuda_kernels.cu` - CUDA kernels
7. `src/cuda/cuda_kernels.h` - CUDA kernel declarations
8. `src/cuda/cuda_my_acceleration.cu` - CUDA acceleration
9. `src/cuda/cuda_routines.cu` - CUDA routines
10. `src/cuda/CalculateRegularAcceleration.cpp`
11. `src/Worker.h` - Worker class
12. `src/QueueScheduler.h` - Queue scheduler
13. `src/SkipList.h` - Skip list
14. `src/performance.h` - Performance tracking
15. `src/profiler.h` - Profiler
16. `src/MPIRoutines.cpp` - MPI routines
17. `src/IrregularRoutines.cpp` - Irregular force routines
18. `src/RegularRoutines.cpp` - Regular force routines
19. `src/RootRoutines.cpp` - Root process routines
20. `src/WorkerRoutines.cpp` - Worker process routines
21. `src/main.cpp` - Main entry point
22. All other .cpp files
