#ifdef FEWBODY
#include "../global.h"
#include "../particle_data.h"
#include <random>
#include <map>

// External SoA container
extern ParticleDataMPI particle_data;

#ifdef SEVN
void Mix(StarSEVN* star1, StarSEVN* star2);
void SetRadius(Particle* ptcl);
#endif


void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time);
void remnantSpinMass(Particle* p1, Particle* p2);
void recoilKick(Particle* p1, Particle* p2);


// made 2024.08.12 by Eunwoo Chung

// Reference: SDAR/sample/AR/ar.cxx & PeTar/src/hard.hpp
void Group::ARIntegration(double next_time) {
    PROFILE_START(TimerID::FewBodyIntegration);

    // Note: particles[] already has current data - no sync needed at entry
    // We sync TO SoA at the end after modifications

    for (int dim=0; dim<DIM; dim++) {
        sym_int.particles.cm.position[dim] = groupCM->position[dim];
        sym_int.particles.cm.velocity[dim] = groupCM->velocity[dim];
        for (int j=0; j<HERMITE_ORDER; j++)
            sym_int.particles.cm.acc_total[dim][j] = groupCM->acc_total[dim][j];
    }
// /*
    if (groupCM->current_time_reg >= groupCM->current_time_irr) { // neighbors were updated in regular routine
        sym_int.particles.cm.num_neighbors = groupCM->num_neighbors;
    }
// */


#ifdef SEVN
#ifdef SEVN_BINARY // this code is not tested yet!!!!!
    // Check if RLOF is triggered
    if (groupCM->spin_param[2] < 0.0) {
        assert(groupCM->num_members == 2);
        assert(groupCM->spin_param[1] == 0.0);
        assert(groupCM->spin_param[0] > 0.0);
        
        auto& bin_root = sym_int.info.getBinaryTreeRoot();

        double ecc_old = bin_root.ecc;
        double a_old = bin_root.semi;

        // New semi-major axis with calcularization conserving the binary angular momentum
        if (a_old > groupCM->spin_param[0])
            bin_root.semi = groupCM->spin_param[0];
        else
            bin_root.semi = a_old * (1 - ecc_old * ecc_old);
        
        bin_root.ecc = 0.0; // Eccentricity should be zero if RLOF is triggered
        bin_root.calcParticles(double(1.0));
        for (int dim = 0; dim < DIM; dim++) {
            bin_root.getLeftMember()->position[dim] += bin_root.position[dim];
            bin_root.getLeftMember()->velocity[dim] += bin_root.velocity[dim];
            bin_root.getRightMember()->position[dim] += bin_root.position[dim];
            bin_root.getRightMember()->velocity[dim] += bin_root.velocity[dim];
        }
        sym_int.initialIntegration(CurrentTime*enzo_time_step); // commented out by EW 2025.9.9 // not tested yet!!!!!
        groupCM->spin_param[2] = 1.0; // This means that RLOF is applied to SDAR
    }
#endif
    bool evolved = false;
    bool kicked = false;
    for (int i=0; i < sym_int.particles.getSize(); i++) {
        Particle* members = &sym_int.particles[i];
        if (members->mass != particles[members->particle_index].mass) {
            members->mass = particles[members->particle_index].mass;
            evolved = true;
            if (particles[members->particle_index].get_binary_interrupt_state() == BinaryInterruptState::kicked) {
                kicked = true;
                // break; // Let's update all particles even if one of them is kicked by stellar evolution by EW 2025.7.5
            }
        }
    }
    if (!kicked && evolved) { // Eunwoo: orbital parameters should be re-calculated due to mass changes during stellar evolution!
        sym_int.particles.shiftToOriginFrame();
        sym_int.info.generateBinaryTree(sym_int.particles,manager.interaction.gravitational_constant);
        sym_int.initialIntegration(CurrentTime*enzo_time_step); // commented out by EW 2025.9.9 // not tested yet!!!!!
    }
    if (kicked) {
        for (int i = 0; i < sym_int.particles.getSize(); i++) {
            Particle* members = &sym_int.particles[i];
            particles[members->particle_index].current_time_irr = CurrentTime;
            // Sync member to SoA
            particle_data.sync_from_particle(particles[members->particle_index],
                                             static_cast<size_t>(members->particle_index));
        }
        isTerminate = true;
        PROFILE_STOP(TimerID::FewBodyIntegration);
        return;
    }
#endif

    assert(next_time > CurrentTime);
    // auto bin_interrupt = sym_int.integrateToTime(next_time*enzo_time_step); // original AR integrator

    energy_binary -= sym_int.getEtot();
    energy_binary_sd -= sym_int.getEtotSlowDown();

// /* // Let's use Kepler solver for unperturbed binary. AR integrator might be very slow if there is a hard binary.
    AR::InterruptBinary<Particle> bin_interrupt;
    if (groupCM->num_neighbors == 0 && groupCM->num_members == 2) { // for unperturbed binary

        bin_interrupt.status = AR::InterruptStatus::none;

        auto& bin_root = sym_int.info.getBinaryTreeRoot();
        bin_root.calcOrbit(double(1.0));

        bin_root.evolve((next_time - CurrentTime)*enzo_time_step);
        bin_root.calcParticles(double(1.0));
        bin_interrupt.time_now = next_time*enzo_time_step;

        if (manager.interrupt_detection_option > 0) {
            Interaction interaction;
            interaction.modifyAndInterruptKepler(bin_interrupt, bin_root, (next_time - CurrentTime)*enzo_time_step);
        }
        // /* // commented out by EW 2025.9.9
        if (bin_interrupt.status == AR::InterruptStatus::none)
            sym_int.initialIntegration(next_time*enzo_time_step);
        // */
    }
    else {
        bin_interrupt = sym_int.integrateToTime(next_time*enzo_time_step);
    }
// */

    energy_binary += sym_int.getEtot();
    energy_binary_sd += sym_int.getEtotSlowDown();

// /* PN corrections
    if (bin_interrupt.status == AR::InterruptStatus::none) { // Every bound orbit

        double delta_Ebin = -sym_int.getEtot();
        double delta_Ebin_SD = -sym_int.getEtotSlowDown();
        
        auto& bin_root = sym_int.info.getBinaryTreeRoot();

        GR_energy_loss_iter(bin_interrupt, bin_root, CurrentTime, next_time);

        // /* // commented out by EW 2025.9.9
        if (bin_interrupt.status == AR::InterruptStatus::none)
            sym_int.initialIntegration(next_time*enzo_time_step); // Eunwoo: this should be fixed later // Eunwoo: I don't think so!
        // */
        
        groupCM->spin_param[1] = bin_root.ecc;

        delta_Ebin += sym_int.getEtot();
        delta_Ebin_SD += sym_int.getEtotSlowDown();

        energy_binary += delta_Ebin;
        energy_binary_sd += delta_Ebin_SD;
        energy_pn -= delta_Ebin;
    }    
// */

    if (bin_interrupt.status != AR::InterruptStatus::none) {

        isMerger = true;
        groupCM->set_binary_interrupt_state(BinaryInterruptState::merger);

        energy_merger += sym_int.getEtot();

        if (sym_int.particles.getSize() == 2) {

            /* // Let's re-calculate member & CM pos/vel in FBTermination!
            double pos[DIM], vel[DIM];

            groupCM->predict_particle_second_order(bin_interrupt.time_now/enzo_time_step - CurrentTime, pos, vel);
            // This might be changed later because changing Pos & Vel during Irregular Acceleration calculation is not good
            // But if SDAR integration is done after Irregular Acceleration calculation, this is fine
            // (Query) by EW 2025.1.6
            for (int dim=0; dim<DIM; dim++) {
                groupCM->position[dim] = pos[dim];
                groupCM->velocity[dim] = vel[dim];
            }
            */
            CurrentTime = bin_interrupt.time_now/enzo_time_step;
            groupCM->current_time_irr = CurrentTime;

            assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
            for (int i = 0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];

                for (int dim=0; dim<DIM; dim++) {
                    particles[members->particle_index].position[dim] = groupCM->position[dim] + members->position[dim];
                    particles[members->particle_index].velocity[dim] = groupCM->velocity[dim] + members->velocity[dim];
                }
                particles[members->particle_index].mass = members->mass;
                particles[members->particle_index].binary_state = members->binary_state;
                particles[members->particle_index].current_time_irr = CurrentTime;

                // Sync member to SoA
                particle_data.sync_from_particle(particles[members->particle_index],
                                                 static_cast<size_t>(members->particle_index));
            }

            isTerminate = true;
            PROFILE_STOP(TimerID::FewBodyIntegration);
            return;
        }
        else {

            CurrentTime = bin_interrupt.time_now/enzo_time_step;

            assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
            for (int i = 0; i < sym_int.particles.getSize(); i++) {
                Particle* members = &sym_int.particles[i];

                for (int dim=0; dim<DIM; dim++) {
                    particles[members->particle_index].position[dim] = groupCM->position[dim] + members->position[dim];
                    particles[members->particle_index].velocity[dim] = groupCM->velocity[dim] + members->velocity[dim];
                }
                particles[members->particle_index].mass = members->mass;
                particles[members->particle_index].binary_state = members->binary_state;

                // Sync member to SoA
                particle_data.sync_from_particle(particles[members->particle_index],
                                                 static_cast<size_t>(members->particle_index));
            }

            // NewFBInitialization3(this);
            PROFILE_STOP(TimerID::FewBodyIntegration);
            return;
        }
    }

    // for write_out_group function by EW 2025.1.6
    assert(!sym_int.particles.isOriginFrame()); // for debugging by EW 2025.1.6
    for (int i = 0; i < sym_int.particles.getSize(); i++) {
        Particle* members = &sym_int.particles[i];

        for (int dim=0; dim<DIM; dim++) {
            particles[members->particle_index].position[dim] = groupCM->position[dim] + members->position[dim];
            particles[members->particle_index].velocity[dim] = groupCM->velocity[dim] + members->velocity[dim];
        }
        particles[members->particle_index].mass = members->mass;
        particles[members->particle_index].current_time_irr = next_time;

        // Sync modified member particle to SoA
        particle_data.sync_from_particle(particles[members->particle_index],
                                         static_cast<size_t>(members->particle_index));
    }

    CurrentTime = next_time;
    PROFILE_STOP(TimerID::FewBodyIntegration);
    return;
}

// made 2024.09.19 by Eunwoo Chung

// reference: tides3.f from Nbody6++GPU (Rizzuto et al. 2020, Arca Sedda et al. 2023)
// Gravitational wave energy loss of hard binaries according to the orbit averaged approximation of Peters & Mathews 1963.
// Calculate average change of energy and angular momentum per orbit.
// Modify semi-major axis, eccentricity, omega(argument of periapsis) per time step in SDAR.
// Orbit shrinking by PN2.5
// Precession by PN1.0 & PN2.0

void GR_energy_loss(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time) {

    const double c = 299752.458 / (velocity_unit / yr * pc / 1e5); // speed of light in code unit
    const double m1 = _bin.m1;
    const double m2 = _bin.m2;
    const double mtot = m1 + m2;
    const double cost = pow(c, -5) * m1 * m2 * mtot;

    double dt = (next_time - current_time) * enzo_time_step;

    double e = _bin.ecc;
    double semi = _bin.semi;
    
    // Define the derivative functions for `e`, `semi`, and `rot_self`
    auto de_dt = [&](double e, double semi) {
        double e2 = e * e;
        double FE2 = e * pow((1 - e2), -2.5) * (1 + 121.0 / 304.0 * e2);
        return 304.0 / 15.0 * cost * pow(semi, -4.0) * FE2;
    };
    auto dsemi_dt = [&](double e, double semi) {
        double e2 = e * e;
        double e4 = e2 * e2;
        double FE1 = pow((1 - e2), -3.5) * (1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e4);
        return 64.0 / 5.0 * cost * pow(semi, -3.0) * FE1;
    };
    auto domega_dt = [&](double e, double semi) {
        double e2 = e * e;
        return (6.0 * M_PI / (c * c * _bin.period) * mtot / (semi * (1 - e2)) +
                3.0 * (18.0 + e2) * M_PI / (2 * pow(c, 4) * _bin.period) * pow((mtot / (semi * (1 - e2))), 2.0));
    };

    // Runge-Kutta 4th Order Method
    // k1
    double k1_de = de_dt(e, semi) * dt;
    double k1_dsemi = dsemi_dt(e, semi) * dt;
    double k1_domega = domega_dt(e, semi) * dt;

    // k2
    double k2_de = de_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;
    double k2_dsemi = dsemi_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;
    double k2_domega = domega_dt(e - 0.5 * k1_de, semi - 0.5 * k1_dsemi) * dt;

    // k3
    double k3_de = de_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;
    double k3_dsemi = dsemi_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;
    double k3_domega = domega_dt(e - 0.5 * k2_de, semi - 0.5 * k2_dsemi) * dt;

    // k4
    double k4_de = de_dt(e - k3_de, semi - k3_dsemi) * dt;
    double k4_dsemi = dsemi_dt(e - k3_de, semi - k3_dsemi) * dt;
    double k4_domega = domega_dt(e - k3_de, semi - k3_dsemi) * dt;

    // Update variables using weighted sum of Runge-Kutta increments
    _bin.ecc -= (k1_de + 2 * k2_de + 2 * k3_de + k4_de) / 6.0;
    _bin.semi -= (k1_dsemi + 2 * k2_dsemi + 2 * k3_dsemi + k4_dsemi) / 6.0;
    _bin.rot_self += (k1_domega + 2 * k2_domega + 2 * k3_domega + k4_domega) / 6.0;  

    // Check for invalid state
    if (!(_bin.ecc > 0) || !(_bin.semi > 0)) {
        fprintf(worker_output_file, "GW driven Merger happened! (a < da)\n");
        fprintf(worker_output_file, "PID: %d and %d\n", _bin.getLeftMember()->pid, _bin.getRightMember()->pid);
        fprintf(worker_output_file, "ecc: %e, semi: %e pc, dsemi: %e pc, timestep: %e Myr\n", _bin.ecc, _bin.semi*position_unit, (k1_dsemi + 2 * k2_dsemi + 2 * k3_dsemi + k4_dsemi) / 6.0*position_unit, dt*1e4);
        fflush(worker_output_file);

        // _bin_interrupt.time_now = current_time * enzo_time_step + dt * num;
        _bin_interrupt.time_now = next_time * enzo_time_step;

        auto* p1 = _bin.getLeftMember();
        auto* p2 = _bin.getRightMember();

        p1->set_binary_interrupt_state(BinaryInterruptState::collision);
        p2->set_binary_interrupt_state(BinaryInterruptState::collision);
        p1->set_binary_pair_id(p2->particle_index);
        p2->set_binary_pair_id(p1->particle_index);
        _bin_interrupt.status = AR::InterruptStatus::merge;
        _bin_interrupt.adr = &_bin;
        return;
    }

    _bin.calcParticles(double(1.0));
    for (int dim = 0; dim < DIM; dim++) {
        _bin.getLeftMember()->position[dim] += _bin.position[dim];
        _bin.getLeftMember()->velocity[dim] += _bin.velocity[dim];
        _bin.getRightMember()->position[dim] += _bin.position[dim];
        _bin.getRightMember()->velocity[dim] += _bin.velocity[dim];
    }
}


void GR_energy_loss_iter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double current_time, double next_time) {
    if (_bin.getMemberN() == 2) {
        _bin.calcOrbit(double(1.0));
        if (_bin.ecc < 1)
            GR_energy_loss(_bin_interrupt, _bin, current_time, next_time);
    }
    else {
        for (int k=0; k<2; k++) {
            if (_bin.isMemberTree(k)) {
                auto memberTree = _bin.getMemberAsTree(k);
                GR_energy_loss_iter(_bin_interrupt, *memberTree, current_time, next_time);
            }
        }
    }
}


void Merge(Particle* p1, Particle* p2) { // Stellar merger

    if (p1->mass < p2->mass)
        std::swap(p1, p2); // p1 should have the larger mass than p2 (p1->mass > p2->mass)

    p1->set_binary_interrupt_state(BinaryInterruptState::none);
    p2->set_binary_interrupt_state(BinaryInterruptState::none);

    double radius;

    if (p1->particle_type > REMNANT && p2->particle_type > REMNANT) {

        // radius = (p1->radius > p2->radius) ? 3*p1->radius : 3*p2->radius; 
        // r_ISCO == 3 * Schwartzschild radius
        radius = (p1->mass >= p2->mass) ? 6*p1->mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2) : 6*p2->mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2);
        // fprintf(merger_output_file, "Separation: %e pc\n", dist(p1->position, p2->position)*position_unit);
        // fprintf(merger_output_file, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
        fprintf(merger_output_file, "r_ISCO: %e pc\n", radius*position_unit);

        fprintf(merger_output_file, "GW driven merger happens!!! (PID: %d, PID: %d)\n", p1->pid, p2->pid);
        fprintf(merger_output_file, "Time: %e Myr\n", p1->current_time_irr*enzo_time_step*1e4);
        // fprintf(merger_output_file, "In center-of-mass frame...\n");
        fprintf(merger_output_file, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->pid, p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
        fprintf(merger_output_file, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->pid, p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "PID: %d. Mass (Msol) - %e, \n", p1->pid, p1->mass*mass_unit);
        fprintf(merger_output_file, "PID: %d. ParticleType - %d\n", p1->pid, p1->particle_type);
        fprintf(merger_output_file, "PID: %d. DIMensionless spin - %e, %e, %e\n", p1->pid, p1->spin_param[0], p1->spin_param[1], p1->spin_param[2]);
        fprintf(merger_output_file, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->pid, p2->position[0]*position_unit, p2->position[1]*position_unit, p2->position[2]*position_unit);
        fprintf(merger_output_file, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->pid, p2->velocity[0]*velocity_unit/yr*pc/1e5, p2->velocity[1]*velocity_unit/yr*pc/1e5, p2->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "PID: %d. Mass (Msol) - %e, \n", p2->pid, p2->mass*mass_unit);
        fprintf(merger_output_file, "PID: %d. ParticleType - %d\n", p2->pid, p2->particle_type);
        fprintf(merger_output_file, "PID: %d. DIMensionless spin - %e, %e, %e\n", p2->pid, p2->spin_param[0], p2->spin_param[1], p2->spin_param[2]);


        double mcm = p1->mass + p2->mass;
        for (int k=0; k<3; k++) {
            p1->position[k] = (p1->mass*p1->position[k] + p2->mass*p2->position[k])/mcm;
            p1->velocity[k] = (p1->mass*p1->velocity[k] + p2->mass*p2->velocity[k])/mcm;
            p2->position[k] = 0.0;
            p2->velocity[k] = 0.0;
        }
        recoilKick(p1, p2);
        remnantSpinMass(p1, p2);
#ifdef SEVN
        SetRadius(p1); // Set radius of the remnant
#else
        p1->radius = 2*p1->mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzschild radius
#endif

        if (p1->particle_type < p2->particle_type) {
            fprintf(merger_output_file, "Warning: p1->particle_type < p2->particle_type in merger! (p1: %d, p2: %d)\n", p1->particle_type, p2->particle_type);
            fprintf(merger_output_file, "We can't trust SEVN data for PID %d!!!\n", p1->pid);
            p1->particle_type = p2->particle_type;
        }

        p2->mass = -1.0;
        // Sync merged particles to SoA
        particle_data.sync_from_particle(*p1, static_cast<size_t>(p1->particle_index));
        particle_data.sync_from_particle(*p2, static_cast<size_t>(p2->particle_index));
        fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
        fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
        fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "Mass (Msol) - %e, \n", p1->mass*mass_unit);
        fprintf(merger_output_file, "ParticleType - %d\n", p1->particle_type);
        fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
    }
    else if ((p1->particle_type > REMNANT && p2->particle_type < REMNANT) ||
            (p1->particle_type < REMNANT && p2->particle_type > REMNANT)) {

        if (p2->particle_type > REMNANT)
            std::swap(p1, p2); // p1 should be compact object

        radius = 1.3*pow((p1->mass + p2->mass)/p2->mass, 1./3)*p2->radius; // TDE radius

        // fprintf(merger_output_file, "Separation: %e pc\n", dist(p1->position, p2->position)*position_unit);
        // fprintf(merger_output_file, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
        fprintf(merger_output_file, "r_TDE: %e pc\n", radius*position_unit);

        fprintf(merger_output_file, "TDE happens!!! (PID: %d, PID: %d)\n", p1->pid, p2->pid);
        fprintf(merger_output_file, "Time: %e Myr\n", p1->current_time_irr*enzo_time_step*1e4);
        // fprintf(merger_output_file, "In center-of-mass frame...\n");
        fprintf(merger_output_file, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->pid, p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
        fprintf(merger_output_file, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->pid, p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "PID: %d. Mass (Msol) - %e, \n", p1->pid, p1->mass*mass_unit);
        fprintf(merger_output_file, "PID: %d. ParticleType - %d\n", p1->pid, p1->particle_type);
        fprintf(merger_output_file, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->pid, p2->position[0]*position_unit, p2->position[1]*position_unit, p2->position[2]*position_unit);
        fprintf(merger_output_file, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->pid, p2->velocity[0]*velocity_unit/yr*pc/1e5, p2->velocity[1]*velocity_unit/yr*pc/1e5, p2->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "PID: %d. Mass (Msol) - %e, \n", p2->pid, p2->mass*mass_unit);
        fprintf(merger_output_file, "PID: %d. ParticleType - %d\n", p2->pid, p2->particle_type);

        double mcm = p1->mass + p2->mass * 0.5; // If TDE happens, the half of the mass of star is accreted to a BH.
        for (int k=0; k<3; k++) {
            p1->position[k] = (p1->mass*p1->position[k] + p2->mass*p2->position[k])/mcm;
            p1->velocity[k] = (p1->mass*p1->velocity[k] + p2->mass*p2->velocity[k])/mcm;
        }

        p1->radius = 2*p1->mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzschild radius in code unit

        p1->delta_mass += 0.5 * p2->mass;
        p1->mass = mcm;
        p2->mass = -1.0;
        // Sync TDE remnants to SoA
        particle_data.sync_from_particle(*p1, static_cast<size_t>(p1->particle_index));
        particle_data.sync_from_particle(*p2, static_cast<size_t>(p2->particle_index));
        fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
        fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
        fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "Mass (Msol) - %e, \n", p1->mass*mass_unit);
        fprintf(merger_output_file, "ParticleType - %d\n", p1->particle_type);
        fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
    }
    else if (p1->particle_type < REMNANT && p2->particle_type < REMNANT) { // Stellar merger

        radius = p1->radius + p2->radius; // Sum of two stellar radius
        // fprintf(merger_output_file, "Separation: %e pc\n", dist(p1->position, p2->position)*position_unit);
        // fprintf(merger_output_file, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
        fprintf(merger_output_file, "r1 + r2: %e pc\n", radius*position_unit);

        fprintf(merger_output_file, "Stellar merger happens!!! (PID: %d, PID: %d)\n", p1->pid, p2->pid);
        fprintf(merger_output_file, "Time: %e Myr\n", p1->current_time_irr*enzo_time_step*1e4);
        // fprintf(merger_output_file, "In center-of-mass frame...\n");
        fprintf(merger_output_file, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->pid, p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
        fprintf(merger_output_file, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->pid, p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "PID: %d. Mass (Msol) - %e, \n", p1->pid, p1->mass*mass_unit);
        fprintf(merger_output_file, "PID: %d. ParticleType - %d\n", p1->pid, p1->particle_type);
        fprintf(merger_output_file, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->pid, p2->position[0]*position_unit, p2->position[1]*position_unit, p2->position[2]*position_unit);
        fprintf(merger_output_file, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->pid, p2->velocity[0]*velocity_unit/yr*pc/1e5, p2->velocity[1]*velocity_unit/yr*pc/1e5, p2->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "PID: %d. Mass (Msol) - %e, \n", p2->pid, p2->mass*mass_unit);
        fprintf(merger_output_file, "PID: %d. ParticleType - %d\n", p2->pid, p2->particle_type);

        double mcm = p1->mass + p2->mass;
        for (int k=0; k<3; k++) {
            p1->position[k] = (p1->mass*p1->position[k] + p2->mass*p2->position[k])/mcm;
            p1->velocity[k] = (p1->mass*p1->velocity[k] + p2->mass*p2->velocity[k])/mcm;
            p2->position[k] = p1->position[k];
            p2->velocity[k] = p1->velocity[k];
        }

#ifdef SEVN
        if (p1->particle_type == NO_FEEDBACK_STAR && p2->particle_type == NO_FEEDBACK_STAR) {

            // p1->delta_mass = mcm - p1->mass;
            // p2->delta_mass = -p2->mass;
            p1->mass = mcm;
            p2->mass = -1.0;

            p1->radius = 2.25461e-8/position_unit*pow(p1->mass*mass_unit, 1./3); // stellar radius in code unit

            if (mcm*mass_unit > 2.2 && mcm*mass_unit < 600) {

                std::stringstream mass;
                mass << std::setprecision(17) << p1->mass*mass_unit; // convert to Msun
                
                std::vector<std::string> init_params{mass.str(), "0.0002", "0.0", "delayed", "zams", "end", "events"};
                size_t id = p1->pid;

                p1->stellar_evolution = new StarSEVN(sevnio, init_params, id, false);
                p1->particle_type = (int)p1->stellar_evolution->getp(Phase::ID);
                p1->formation_time = p1->current_time_irr*enzo_time_step*1e4;
                p1->world_time = p1->current_time_irr*enzo_time_step*1e4;

                SEVNList.insert({p1->world_time + p1->stellar_evolution->getp(Timestep::ID), p1->particle_index});

                SetRadius(p1);
                fprintf(stdout, "New Star class made!\n");
                fprintf(stdout, "PID: %d. Mass: %e Msol, Radius: %e pc\n", p1->pid, p1->mass*mass_unit, p1->radius*position_unit);
            }
            fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
            fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
            fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(merger_output_file, "Mass (Msol) - %e, \n", p1->mass*mass_unit);
            fprintf(merger_output_file, "ParticleType - %d\n", p1->particle_type);
            fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
        }
        else if (p1->particle_type != NO_FEEDBACK_STAR && p2->particle_type != NO_FEEDBACK_STAR) {

            fprintf(stdout, "Before Mix... p1 (PID: %d). ParticleType: %d, p2 (PID: %d). ParticleType: %d\n", p1->pid, p1->particle_type, p2->pid, p2->particle_type);

            Mix(p1->stellar_evolution, p2->stellar_evolution);
            fprintf(stdout, "Mix done!\n");

            if (p1->stellar_evolution->amiremnant())
                p1->particle_type = REMNANT + (int)p1->stellar_evolution->getp(RemnantType::ID);
            else
                p1->particle_type = (int)p1->stellar_evolution->getp(Phase::ID);

            if (p2->stellar_evolution->amiremnant())
                p2->particle_type = REMNANT + (int)p2->stellar_evolution->getp(RemnantType::ID);
            else
                p2->particle_type = (int)p2->stellar_evolution->getp(Phase::ID);

            if (p1->stellar_evolution->amiempty() && !p2->stellar_evolution->amiempty()) {
                p1->mass = -1.0;

                p2->mass = p2->stellar_evolution->getp(Mass::ID)/mass_unit;

                fprintf(stdout, "After Mix... p1 (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", p1->pid, p1->mass*mass_unit, p1->stellar_evolution->get_zams());
                fprintf(stdout, "After Mix... p2 (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", p2->pid, p2->mass*mass_unit, p2->stellar_evolution->get_zams());
                fprintf(stdout, "p1: amiempty(): %d\n", p1->stellar_evolution->amiempty());
                fprintf(stdout, "p2: amiempty(): %d\n", p2->stellar_evolution->amiempty());
                // p2->delta_mass += p1->delta_mass // not yet by EW 2025.1.20
                // p1->delta_mass = 0.0; // not yet by EW 2025.1.20
                SetRadius(p2);

                fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
                fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p2->position[0]*position_unit, p2->position[1]*position_unit, p2->position[2]*position_unit);
                fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->velocity[0]*velocity_unit/yr*pc/1e5, p2->velocity[1]*velocity_unit/yr*pc/1e5, p2->velocity[2]*velocity_unit/yr*pc/1e5);
                fprintf(merger_output_file, "Mass (Msol) - %e, \n", p2->mass*mass_unit);
                fprintf(merger_output_file, "ParticleType - %d\n", p2->particle_type);
                fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else if (!p1->stellar_evolution->amiempty() && p2->stellar_evolution->amiempty()) {
                p2->mass = -1.0;

                p1->mass = p1->stellar_evolution->getp(Mass::ID)/mass_unit;

                fprintf(stdout, "After Mix... p1 (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", p1->pid, p1->mass*mass_unit, p1->stellar_evolution->get_zams());
                fprintf(stdout, "After Mix... p2 (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", p2->pid, p2->mass*mass_unit, p2->stellar_evolution->get_zams());
                fprintf(stdout, "p1: amiempty(): %d\n", p1->stellar_evolution->amiempty());
                fprintf(stdout, "p2: amiempty(): %d\n", p2->stellar_evolution->amiempty());
                // p1->delta_mass += p2->delta_mass // not yet by EW 2025.1.20
                // p2->delta_mass = 0.0; // not yet by EW 2025.1.20
                SetRadius(p1);

                fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
                fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
                fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
                fprintf(merger_output_file, "Mass (Msol) - %e, \n", p1->mass*mass_unit);
                fprintf(merger_output_file, "ParticleType - %d\n", p1->particle_type);
                fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else if (p1->stellar_evolution->amiempty() && p2->stellar_evolution->amiempty()) { // Type Ia supernova
                p1->delta_mass += p1->mass;
                p1->mass = -1.0;
                p2->delta_mass += p2->mass;
                p2->mass = -1.0;
                fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
                fprintf(merger_output_file, "Type Ia Supernova event! Both of the stars becomes empty!\n");
                fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else
                throw std::runtime_error("None of stars are empty: Something wrong in stellar merger!");
        }
        else if (p1->particle_type == NO_FEEDBACK_STAR && p2->particle_type != NO_FEEDBACK_STAR) {

            if (!p2->stellar_evolution->amiremnant()) {
                p2->stellar_evolution->update_from_binary(Mass::ID, p1->mass*mass_unit);
                p2->stellar_evolution->update_from_binary(dMcumul_binary::ID, p1->mass*mass_unit);
                if (p2->stellar_evolution->aminakedhelium())
                    p2->stellar_evolution->jump_to_normal_tracks();
                else
                    p2->stellar_evolution->find_new_track_after_merger();

                p2->mass = p2->stellar_evolution->getp(Mass::ID)/mass_unit;
                SetRadius(p2);
                p2->particle_type = (int)p2->stellar_evolution->getp(Phase::ID);
                p1->mass = -1.0;
            }            

            fprintf(stdout, "Mix with no done!\n");
            fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
            fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p2->position[0]*position_unit, p2->position[1]*position_unit, p2->position[2]*position_unit);
            fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->velocity[0]*velocity_unit/yr*pc/1e5, p2->velocity[1]*velocity_unit/yr*pc/1e5, p2->velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(merger_output_file, "Mass (Msol) - %e, \n", p2->mass*mass_unit);
            fprintf(merger_output_file, "ParticleType - %d\n", p2->particle_type);
            fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
        }
        else if (p1->particle_type != NO_FEEDBACK_STAR && p2->particle_type == NO_FEEDBACK_STAR) {

            if (!p1->stellar_evolution->amiremnant()) {
                p1->stellar_evolution->update_from_binary(Mass::ID, p2->mass*mass_unit);
                p1->stellar_evolution->update_from_binary(dMcumul_binary::ID, p2->mass*mass_unit);
                if (p1->stellar_evolution->aminakedhelium())
                    p1->stellar_evolution->jump_to_normal_tracks();
                else
                    p1->stellar_evolution->find_new_track_after_merger();

                p1->mass = p1->stellar_evolution->getp(Mass::ID)/mass_unit;
                SetRadius(p1);
                p1->particle_type = (int)p1->stellar_evolution->getp(Phase::ID);
                p2->mass = -1.0;
            }

            fprintf(stdout, "Mix with no done!\n");
            fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
            fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
            fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(merger_output_file, "Mass (Msol) - %e, \n", p1->mass*mass_unit);
            fprintf(merger_output_file, "ParticleType - %d\n", p1->particle_type);
            fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
        }
    }
    fflush(merger_output_file);
    fflush(stdout);
#else
        p1->delta_mass = mcm - p1->mass;
        p2->delta_mass = -p2->mass;
        p1->mass = mcm;
        p2->mass = -1.0;
        p1->radius = 2.25461e-8/position_unit*pow(p1->mass*mass_unit, 1./3); // stellar radius in code unit

        fprintf(merger_output_file, "---------------Merger remnant properties---------------\n");
        fprintf(merger_output_file, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->position[0]*position_unit, p1->position[1]*position_unit, p1->position[2]*position_unit);
        fprintf(merger_output_file, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
        fprintf(merger_output_file, "Mass (Msol) - %e, \n", p1->mass*mass_unit);
        fprintf(merger_output_file, "---------------------END-OF-MERGER---------------------\n\n");
        // Sync stellar merger remnants to SoA (non-SEVN path)
        particle_data.sync_from_particle(*p1, static_cast<size_t>(p1->particle_index));
        particle_data.sync_from_particle(*p2, static_cast<size_t>(p2->particle_index));
    }
    fflush(merger_output_file);
#endif
}

// Reference for remnant spin: Hofmann et al. (2016) (https://iopscience.iop.org/article/10.3847/2041-8205/825/2/L19/pdf)
// Using eq (2)-(6), (13)-(16)
// n_M = 3, n_J = 4
// Reference for remnant mass: Barausse et al. (2012) (https://iopscience.iop.org/article/10.1088/0004-637X/758/1/63/pdf)
// Using eq (1)-(5), (12), (15)-(18)
void remnantSpinMass(Particle* p1, Particle* p2) {

    double k[4][5] = {{-5.9, 3.39221, 4.48865, -5.77101, -13.0459},
                    {35.1287, -72.9336, -86.0036, 93.7371, 200.975},
                    {-146.822, 387.184, 447.009, -467.383, -884.339},
                    {223.911, -648.502, -697.177, 753.738, 1166.89}};
    double ksi = 0.474046;

    double a1 = sqrt(p1->spin_param[0]*p1->spin_param[0] + p1->spin_param[1]*p1->spin_param[1] + p1->spin_param[2]*p1->spin_param[2]);
    double a2 = sqrt(p2->spin_param[0]*p2->spin_param[0] + p2->spin_param[1]*p2->spin_param[1] + p2->spin_param[2]*p2->spin_param[2]);
    double Mtot = p1->mass + p2->mass;
    double q = p2->mass/p1->mass;
    assert(q <= 1);
    double nu = q/(1+q)/(1+q);
    const double c = 299752.458 / (velocity_unit / yr * pc / 1e5);

    double pos_rel[3];
    double vel_rel[3];
    double L_ang[3]; // specific angular momentum (r_rel x v_rel)

    for (int dim=0; dim<DIM; dim++) {
        pos_rel[dim] = p1->position[dim] - p2->position[dim];
        vel_rel[dim] = p1->velocity[dim] - p2->velocity[dim];
    }
    L_ang[0] = pos_rel[1] * vel_rel[2] - pos_rel[2] * vel_rel[1];
    L_ang[1] = pos_rel[2] * vel_rel[0] - pos_rel[0] * vel_rel[2];
    L_ang[2] = pos_rel[0] * vel_rel[1] - pos_rel[1] * vel_rel[0];

    fprintf(merger_output_file, "L_orbit: (%e, %e, %e)\n", (p1->mass*p2->mass/Mtot) * L_ang[0],
                                                    (p1->mass*p2->mass/Mtot) * L_ang[1],
                                                    (p1->mass*p2->mass/Mtot) * L_ang[2]);
    fprintf(merger_output_file, "S_1: (%e, %e, %e)\n", (p1->mass*p1->mass/c) * p1->spin_param[0], 
                                                (p1->mass*p1->mass/c) * p1->spin_param[1], 
                                                (p1->mass*p1->mass/c) * p1->spin_param[2]);
    fprintf(merger_output_file, "S_2: (%e, %e, %e)\n", (p2->mass*p2->mass/c) * p2->spin_param[0], 
                                                (p2->mass*p2->mass/c) * p2->spin_param[1], 
                                                (p2->mass*p2->mass/c) * p2->spin_param[2]);
    fprintf(merger_output_file, "In code unit!\n");

    auto cosine = [&](double a[3], double b[3]) {
        double a_mag = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
        double b_mag = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
        if (a_mag == 0. || b_mag == 0)
            return 0.;
        else
            return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])/a_mag/b_mag;
    };

    double cosa = cosine(p1->spin_param, p2->spin_param); // cos(alpha): This angle should be changed to the initial value. I will change this later.
    double cosb = cosine(L_ang, p1->spin_param);      // cos(beta)
    double cosg = cosine(L_ang, p2->spin_param);      // cos(gamma)

    double atot = (a1*cosb + a2*cosg*q*q)/(1 + q)/(1 + q); // atilde in Barausse et al. (2012)
    double aeff = atot + ksi*nu*(a1*cosb + a2*cosg);

    auto Z1 = [&](double a) {
        return 1 + pow(1 - a*a, 1./3)*(pow(1 + a, 1./3) + pow(1 - a, 1./3));
    };
    auto Z2 = [&](double a, double Z1) {
        return sqrt(3*a*a + Z1*Z1);
    };
    auto r_ISCO = [&](double a, double Z1, double Z2) {
        if (a > 0)
            return 3 + Z2 - sqrt(3 - Z1)*sqrt(3 + Z1 + 2*Z2);
        else if (a < 0)
            return 3 + Z2 + sqrt(3 - Z1)*sqrt(3 + Z1 + 2*Z2);
        else
            return 3 + Z2;
    };
    auto E_ISCO = [&](double r) {
        return sqrt(1 - 2./3/r);
    };

    double Z1_aeff = Z1(aeff);
    double Z2_aeff = Z2(aeff, Z1_aeff);
    double r_ISCO_aeff = r_ISCO(aeff, Z1_aeff, Z2_aeff);
    double E_ISCO_aeff = E_ISCO(r_ISCO_aeff);
    double L_ISCO_aeff = 2/3/sqrt(3)*(1 + 2*sqrt(3*r_ISCO_aeff - 2));

    double l = L_ISCO_aeff - 2 * atot * (E_ISCO_aeff - 1);
    for (int i = 0; i <= 3; ++i) {      // n_M = 3
        for (int j = 0; j <= 4; ++j) {  // n_J = 4
            l += k[i][j] * pow(nu, 1 + i) * pow(aeff, j);
        }
    }
    l = abs(l);

    double afin = 1./pow(1 + q, 2)*sqrt(a1*a1 + a2*a2*pow(q, 4) + 2*a1*a2*q*q*cosa + 2*(a1*cosb + a2*q*q*cosg)*l*q + l*l*q*q);

    double J_tot[3];
    for (int i=0; i<3; i++)
        J_tot[i] = (p1->mass*p2->mass/Mtot) * L_ang[i] + (p1->mass*p1->mass/c) * p1->spin_param[i] + (p2->mass*p2->mass/c) * p2->spin_param[i];

    double J_tot_norm[3];
    double J_tot_mag = sqrt(J_tot[0]*J_tot[0] + J_tot[1]*J_tot[1] + J_tot[2]*J_tot[2]);
    for (int i=0; i<3; i++)
        J_tot_norm[i] = J_tot[i]/J_tot_mag;

    for (int i=0; i<3; i++)
        p1->spin_param[i] = afin*J_tot_norm[i];

    fprintf(merger_output_file, "DIMensionless spin of remnant BH: (%e, %e, %e)\n", p1->spin_param[0], p1->spin_param[1], p1->spin_param[2]);

    double Z1_atot = Z1(atot);
    double Z2_atot = Z2(atot, Z1_atot);
    double r_ISCO_atot = r_ISCO(atot, Z1_atot, Z2_atot);
    double E_ISCO_atot = E_ISCO(r_ISCO_atot);

    double Erad = (1 - E_ISCO_atot)*nu + 4*nu*nu*(4*0.04827 + 16*0.01707*atot*(atot+1) + E_ISCO_atot - 1);

    p1->mass = (1 - Erad) * Mtot;

    fprintf(merger_output_file, "Mass of remnant BH: %e Msol\n", p1->mass*mass_unit);
    fflush(merger_output_file);
}

// Reference: Arca Sedda et al. (2020) (https://iopscience.iop.org/article/10.3847/1538-4357/ab88b2/pdf)
void recoilKick(Particle* p1, Particle* p2) {
    
    double A      = 1.2e4; // km/s
    double B      = -0.93;
    double H      = 6.9e3; // km/s
    double ksi    = 145 * M_PI / 180; // rad
    double V11    = 3677.76; // km/s
    double VA     = 2.481e3; // km/s
    double VB     = 1.793e3; // km/s
    double VC     = 1.507e3; // km/s

    auto cosine = [&](double a[3], double b[3]) {
        double a_mag = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
        double b_mag = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
        if (a_mag == 0. || b_mag == 0.)
            return 0.;
        else
            return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])/a_mag/b_mag;
    };

    double a1 = sqrt(p1->spin_param[0]*p1->spin_param[0] + p1->spin_param[1]*p1->spin_param[1] + p1->spin_param[2]*p1->spin_param[2]);
    double a2 = sqrt(p2->spin_param[0]*p2->spin_param[0] + p2->spin_param[1]*p2->spin_param[1] + p2->spin_param[2]*p2->spin_param[2]);
    double q = p2->mass/p1->mass;
    // fprintf(merger_output_file, "q: %e\n", q);
    assert(q <= 1);
    double nu = q/(1+q)/(1+q);
    // fprintf(merger_output_file, "nu: %e\n", nu);

    double pos_rel[3];
    double vel_rel[3];
    double L_ang[3]; // specific angular momentum (r_rel x v_rel)

    for (int dim=0; dim<DIM; dim++) {
        pos_rel[dim] = p1->position[dim] - p2->position[dim];
        vel_rel[dim] = p1->velocity[dim] - p2->velocity[dim];
    }
    L_ang[0] = pos_rel[1] * vel_rel[2] - pos_rel[2] * vel_rel[1];
    L_ang[1] = pos_rel[2] * vel_rel[0] - pos_rel[0] * vel_rel[2];
    L_ang[2] = pos_rel[0] * vel_rel[1] - pos_rel[1] * vel_rel[0];


    double cosa = cosine(p1->spin_param, p2->spin_param); // cos(alpha): This angle should be changed to the initial value. I will change this later.
    double sina = sin(acos(cosa));
    double cosb = cosine(L_ang, p1->spin_param);      // cos(beta)
    double sinb = sin(acos(cosb));
    double cosg = cosine(L_ang, p2->spin_param);      // cos(gamma)
    double sing = sin(acos(cosg));
    // fprintf(merger_output_file, "cosa: %e, sina: %e, cosb: %e, sinb: %e, cosg: %e, sing: %e\n", cosa, sina, cosb, sinb, cosg, sing);

    double a2par  = a2 * cosg;
    // fprintf(merger_output_file, "a2par: %e\n", a2par);
    double a2per1 = a2 * sing;
    // fprintf(merger_output_file, "a2per1, %e\n", a2per1);
    double a2per2 = 0.0;
    // fprintf(merger_output_file, "a2per2, %e\n", a2per2);
    
    double a1par  = a1 * cosb;
    // fprintf(merger_output_file, "a1par: %e\n", a1par);
    double a1per1 = a1 * sinb*cosa;
    // fprintf(merger_output_file, "a1per1: %e\n", a1per1);
    double a1per2 = a1 * sinb*sina;
    // fprintf(merger_output_file, "a1per2: %e\n", a1per2);

    double KSIpar = 2 * (a2par + q*q*a1par) / (1 + q) / (1 + q);
    // fprintf(merger_output_file, "KSIpar: %e\n", KSIpar);
    std::random_device rd; // Obtain a random number from hardware
    std::mt19937 mt(rd()); // Seed the generator
    std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
    double phi = 2 * M_PI * distr(mt); // phi_Delta - phi_1
    // fprintf(merger_output_file, "phi: %e\n", phi);

    double vm = A*nu*nu*sqrt(1 - 4*nu) * (1 + B*nu);
    // fprintf(merger_output_file, "vm: %e\n", vm);
    double vper = H*nu*nu / (1 + q) * (a2par - q * a1par);
    // fprintf(merger_output_file, "vper: %e\n", vper);
    double vpar = 16*nu*nu / (1 + q) * (V11 + VA*KSIpar + VB*KSIpar*KSIpar + VC*KSIpar*KSIpar*KSIpar);
    vpar *= sqrt((a2per1 - q * a1per1) * (a2per1 - q * a1per1) + (a2per2 - q * a1per2) * (a2per2 - q * a1per2)) * cos(phi);
    // fprintf(merger_output_file, "vpar: %e\n", vpar);

    double L_ang_mag = sqrt(L_ang[0]*L_ang[0] + L_ang[1]*L_ang[1] + L_ang[2]*L_ang[2]);
    double e_par[3]   = {L_ang[0]/L_ang_mag, L_ang[1]/L_ang_mag, L_ang[2]/L_ang_mag};
    double e_per1[3]  = {p2->spin_param[0] - p2->spin_param[0]*cosg, p2->spin_param[1] - p2->spin_param[1]*cosg, p2->spin_param[2] - p2->spin_param[2]*cosg};
    double norm1      = sqrt(e_per1[0]*e_per1[0] + e_per1[1]*e_per1[1] + e_per1[2]*e_per1[2]);
    if (norm1 != 0) {
        for (int i=0; i<3; i++)
            e_per1[i]   /= norm1; 
    }
            
    double e_per2[3]  =   {e_par[1] * e_per1[2] - e_par[2] * e_per1[1], 
                        e_par[2] * e_per1[0] - e_par[0] * e_per1[2], 
                        e_par[0] * e_per1[1] - e_par[1] * e_per1[0]}; // cross product: e2 = e3 x e1
    // fprintf(merger_output_file, "e1: %e, %e, %e\n", e_per1[0], e_per1[1], e_per1[2]);
    // fprintf(merger_output_file, "e2: %e, %e, %e\n", e_per2[0], e_per2[1], e_per2[2]);
    // fprintf(merger_output_file, "e3: %e, %e, %e\n", e_par[0], e_par[1], e_par[2]);

    double vkick[3];
    for (int i=0; i<3; i++)
        vkick[i] = (vm + vper*cos(ksi)) * e_per1[i] + vper*sin(ksi) * e_per2[i] + vpar * e_par[i];

    fprintf(merger_output_file, "GW recoil kick: (%e, %e, %e) km/s\n", vkick[0], vkick[1], vkick[2]);
    fprintf(merger_output_file, "\t magnitude: %e km/s\n", sqrt(vkick[0]*vkick[0] + vkick[1]*vkick[1] + vkick[2]*vkick[2]));

    for (int i=0; i<3; i++)
        p1->velocity[i] += vkick[i]/(velocity_unit/yr*pc/1e5); // km/s to code unit

    // fprintf(merger_output_file, "Remnant velocity: (%e, %e, %e) km/s\n", p1->velocity[0]*velocity_unit/yr*pc/1e5, p1->velocity[1]*velocity_unit/yr*pc/1e5, p1->velocity[2]*velocity_unit/yr*pc/1e5);
    fflush(merger_output_file);
}

#endif