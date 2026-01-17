#ifdef FEWBODY
#include "../global.h"
#include <unordered_set>

// calculate dr dv of a pair
Float calcDrDv(const double *pos1, const double *pos2, const double *vel1, const double *vel2) {
    Float dx[3],dv[3];
    dx[0] = pos1[0] - pos2[0];
    dx[1] = pos1[1] - pos2[1];
    dx[2] = pos1[2] - pos2[2];

    dv[0] = vel1[0] - vel2[0];
    dv[1] = vel1[1] - vel2[1];
    dv[2] = vel1[2] - vel2[2];

    Float drdv= dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];

    return drdv;
}

//! check the pair with distance below r_crit for ptcl in adr_dt_sorted_
/*! First check nearest neighbor distance r_min
    If r_min<r_crit, check the direction, if income, accept as group
*/
void Particle::check_new_group() {

    const Float kappa_org_crit = 1e-2; // kappa_org criterion for new group kappa_org>kappa_org_crit
    const double r_crit = r_search; // distance criterion

    double pos1[DIM], vel1[DIM];

    Particle* ptcl2;

    std::unordered_set<int> CMPtclsSet;

    // check only active particles 
    // single case
    for (int i=0; i < this->num_neighbors; i++) {
        ptcl2 = &particles[neighbors[this->neighbors_offset + i]];
        if (!ptcl2->is_active) {
            if (ptcl2->cm_particle_index != -1) {
                CMPtclsSet.insert(ptcl2->cm_particle_index);
            }
            continue;
        }

        // if (ptcl2->time_step_irr > this->time_step_irr) // test_1e5_4 & 5: this must make the same result!
        if (ptcl2->time_step_irr > t_search) // fiducial: 1e-5 but for RSEARCH = 0.00025 pc, 1e-6 Myr seems good
            continue;

        double time_ptcl2 = (ptcl2->num_neighbors == 0) ? ptcl2->current_time_reg : ptcl2->current_time_irr;

        double dt = std::abs(this->current_time_irr - time_ptcl2);

        double pos2[DIM], vel2[DIM];
        
		this->predict_particle_second_order(dt, pos1, vel1);
		ptcl2->predict_particle_second_order(dt, pos2, vel2);

        const Float dr = sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
        
        if (dr < r_crit) {

            Float drdv = calcDrDv(pos1, pos2, vel1, vel2);
            // only inwards
            if(drdv<0.0) {
// /* // test_1e4_2
                Float fcm[3] = {this->mass*this->acc_total[0][0] + ptcl2->mass*ptcl2->acc_total[0][0], 
                                this->mass*this->acc_total[1][0] + ptcl2->mass*ptcl2->acc_total[1][0], 
                                this->mass*this->acc_total[2][0] + ptcl2->mass*ptcl2->acc_total[2][0]};

                AR::SlowDown sd;
                Interaction interaction;
                Float mcm = this->mass + ptcl2->mass;

                // sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
                sd.initialSlowDownReference(1e-6, NUMERIC_FLOAT_MAX);

                sd.pert_in = interaction.calcPertFromMR(dr, this->mass, ptcl2->mass);
                sd.pert_out = interaction.calcPertFromForce(fcm, mcm, mcm);

                sd.calcSlowDownFactor();
                Float kappa_org = sd.getSlowDownFactorOrigin();

                // avoid strong perturbed case, estimate perturbation
                // if kappa_org < criterion, avoid to form new group, should be consistent as checkbreak
                if(kappa_org<kappa_org_crit) continue;
// */ // test_1e4_2
                this->new_members[this->new_num_members++] = ptcl2->particle_index;
            }
        }
    }
    for (int i: CMPtclsSet) {
        ptcl2 = &particles[i];

        if (this->pid == ptcl2->pid) {
            continue;
        }

        if (!ptcl2->is_active) {
            fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", this->pid, ptcl2->pid);
            assert(ptcl2->is_active);
        }

        // if (ptcl2->time_step_irr > this->time_step_irr) // test_1e5_4 & 5: this must make the same result!
        if (ptcl2->time_step_irr > t_search) // fiducial: 1e-5 but for RSEARCH = 0.00025 pc, 1e-6 Myr seems good
            continue;

        double time_ptcl2 = (ptcl2->num_neighbors == 0) ? ptcl2->current_time_reg : ptcl2->current_time_irr;

        double dt = std::abs(this->current_time_irr - time_ptcl2);

        double pos2[DIM], vel2[DIM];
        
        this->predict_particle_second_order(dt, pos1, vel1);
        ptcl2->predict_particle_second_order(dt, pos2, vel2);

        const Float dr = sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
        
        if (dr < r_crit) {

            Float drdv = calcDrDv(pos1, pos2, vel1, vel2);
            // only inwards
            if(drdv<0.0) {
// /* // test_1e4_2
                Float fcm[3] = {this->mass*this->acc_total[0][0] + ptcl2->mass*ptcl2->acc_total[0][0], 
                                this->mass*this->acc_total[1][0] + ptcl2->mass*ptcl2->acc_total[1][0], 
                                this->mass*this->acc_total[2][0] + ptcl2->mass*ptcl2->acc_total[2][0]};

                AR::SlowDown sd;
                Interaction interaction;
                Float mcm = this->mass + ptcl2->mass;

                // sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
                sd.initialSlowDownReference(1e-6, NUMERIC_FLOAT_MAX);

                sd.pert_in = interaction.calcPertFromMR(dr, this->mass, ptcl2->mass);
                sd.pert_out = interaction.calcPertFromForce(fcm, mcm, mcm);

                sd.calcSlowDownFactor();
                Float kappa_org = sd.getSlowDownFactorOrigin();

                // avoid strong perturbed case, estimate perturbation
                // if kappa_org < criterion, avoid to form new group, should be consistent as checkbreak
                if(kappa_org<kappa_org_crit) continue;
// */ // test_1e4_2
                this->new_members[this->new_num_members++] = i;
            }
        }
    }
}

// Use when many-body (>3) group broke
// Or primordial binary search
void Particle::check_new_group_v2() {
    
    const Float kappa_org_crit = 1e-2; // kappa_org criterion for new group kappa_org>kappa_org_crit
    const double r_crit = r_search; // distance criterion

    double pos1[DIM], vel1[DIM];

    Particle* ptcl2;

    std::unordered_set<int> CMPtclsSet;

    // check only active particles 
    // single case
    for (int i=0; i < this->num_neighbors; i++) {
        ptcl2 = &particles[neighbors[this->neighbors_offset + i]];
        if (!ptcl2->is_active) {
            if (ptcl2->cm_particle_index != -1) {
                CMPtclsSet.insert(ptcl2->cm_particle_index);
            }
            continue;
        }

        double time_ptcl2 = (ptcl2->num_neighbors == 0) ? ptcl2->current_time_reg : ptcl2->current_time_irr;

        double dt = std::abs(this->current_time_irr - time_ptcl2);

        double pos2[DIM], vel2[DIM];
        
		this->predict_particle_second_order(dt, pos1, vel1);
		ptcl2->predict_particle_second_order(dt, pos2, vel2);

        const Float dr = sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]));

        double v2 = (vel1[0]-vel2[0])*(vel1[0]-vel2[0]) + (vel1[1]-vel2[1])*(vel1[1]-vel2[1]) + (vel1[2]-vel2[2])*(vel1[2]-vel2[2]);
        double energy = v2/2 - (this->mass + ptcl2->mass)/dr; // determine they are bound or not
        
        if (dr < r_crit && energy < 0) {
            this->new_members[this->new_num_members++] = ptcl2->particle_index;
        }
    }
    for (int i: CMPtclsSet) {
        ptcl2 = &particles[i];

        if (this->pid == ptcl2->pid) {
            continue;
        }

        if (!ptcl2->is_active) {
            fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", this->pid, ptcl2->pid);
            assert(ptcl2->is_active);
        }

        double time_ptcl2 = (ptcl2->num_neighbors == 0) ? ptcl2->current_time_reg : ptcl2->current_time_irr;

        double dt = std::abs(this->current_time_irr - time_ptcl2);

        double pos2[DIM], vel2[DIM];
        
        this->predict_particle_second_order(dt, pos1, vel1);
        ptcl2->predict_particle_second_order(dt, pos2, vel2);

        const Float dr = sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]));

        double v2 = (vel1[0]-vel2[0])*(vel1[0]-vel2[0]) + (vel1[1]-vel2[1])*(vel1[1]-vel2[1]) + (vel1[2]-vel2[2])*(vel1[2]-vel2[2]);
        double energy = v2/2 - (this->mass + ptcl2->mass)/dr; // determine they are bound or not
        
        if (dr < r_crit && energy < 0) {
            this->new_members[this->new_num_members++] = i;
        }
    }
}

void Particle::check_new_group_v3() {

    const Float kappa_org_crit = 1e-2; // kappa_org criterion for new group kappa_org>kappa_org_crit
    Particle* ptcl2;

    int NumberOfGroupCandidate = this->new_num_members;

    this->new_num_members = 0;

    double pos1[DIM], vel1[DIM];
    double pos2[DIM], vel2[DIM];

    for (int i=0; i < NumberOfGroupCandidate; i++) {
        ptcl2 = &particles[this->new_members[i]];

        double time_ptcl2 = (ptcl2->num_neighbors == 0) ? ptcl2->current_time_reg : ptcl2->current_time_irr;

        double dt = std::abs(this->current_time_irr - time_ptcl2);

        this->predict_particle_second_order(dt, pos1, vel1);
        ptcl2->predict_particle_second_order(dt, pos2, vel2);

        const Float dr = sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]));

        if (dr > r_search) continue;

        Float fcm[3] = {this->mass*this->acc_total[0][0] + ptcl2->mass*ptcl2->acc_total[0][0], 
        this->mass*this->acc_total[1][0] + ptcl2->mass*ptcl2->acc_total[1][0], 
        this->mass*this->acc_total[2][0] + ptcl2->mass*ptcl2->acc_total[2][0]};
    
        AR::SlowDown sd;
        Interaction interaction;
        Float mcm = this->mass + ptcl2->mass;
    
        // sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
        sd.initialSlowDownReference(1e-6, NUMERIC_FLOAT_MAX);
    
        sd.pert_in = interaction.calcPertFromMR(dr, this->mass, ptcl2->mass);
        sd.pert_out = interaction.calcPertFromForce(fcm, mcm, mcm);
    
        sd.calcSlowDownFactor();
        Float kappa_org = sd.getSlowDownFactorOrigin();
    
        // avoid strong perturbed case, estimate perturbation
        // if kappa_org < criterion, avoid to form new group, should be consistent as checkbreak
        if(kappa_org<kappa_org_crit) continue;

        this->new_members[this->new_num_members++] = ptcl2->particle_index;
    }
}

void Particle::check_new_group_v4() {

    const double r_crit = r_search; // distance criterion

    double pos1[DIM], vel1[DIM];

    Particle* ptcl2;

    std::unordered_set<int> CMPtclsSet;

    // check only active particles 
    // single case
    for (int i=0; i < this->num_neighbors; i++) {
        ptcl2 = &particles[neighbors[this->neighbors_offset + i]];
        if (!ptcl2->is_active) {
            if (ptcl2->cm_particle_index != -1) {
                CMPtclsSet.insert(ptcl2->cm_particle_index);
            }
            continue;
        }

        // if (ptcl2->time_step_irr > this->time_step_irr) // test_1e5_4 & 5: this must make the same result!
        if (ptcl2->time_step_irr > t_search) // fiducial: 1e-5 but for RSEARCH = 0.00025 pc, 1e-6 Myr seems good
            continue;

        double time_ptcl2 = (ptcl2->num_neighbors == 0) ? ptcl2->current_time_reg : ptcl2->current_time_irr;

        double dt = std::abs(this->current_time_irr - time_ptcl2);

        double pos2[DIM], vel2[DIM];
        
		this->predict_particle_second_order(dt, pos1, vel1);
		ptcl2->predict_particle_second_order(dt, pos2, vel2);

        const Float dr = sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
        
        if (dr < r_crit) {

            Float drdv = calcDrDv(pos1, pos2, vel1, vel2);
            // only inwards
            if(drdv<0.0) {
                this->new_members[this->new_num_members++] = ptcl2->particle_index;
            }
        }
    }
    for (int i: CMPtclsSet) {
        ptcl2 = &particles[i];

        if (this->pid == ptcl2->pid) {
            continue;
        }

        if (!ptcl2->is_active) {
            fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", this->pid, ptcl2->pid);
            assert(ptcl2->is_active);
        }

        // if (ptcl2->time_step_irr > this->time_step_irr) // test_1e5_4 & 5: this must make the same result!
        if (ptcl2->time_step_irr > t_search) // fiducial: 1e-5 but for RSEARCH = 0.00025 pc, 1e-6 Myr seems good
            continue;

        double time_ptcl2 = (ptcl2->num_neighbors == 0) ? ptcl2->current_time_reg : ptcl2->current_time_irr;

        double dt = std::abs(this->current_time_irr - time_ptcl2);

        double pos2[DIM], vel2[DIM];
        
        this->predict_particle_second_order(dt, pos1, vel1);
        ptcl2->predict_particle_second_order(dt, pos2, vel2);

        const Float dr = sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]));
        
        if (dr < r_crit) {

            Float drdv = calcDrDv(pos1, pos2, vel1, vel2);
            // only inwards
            if(drdv<0.0) {
                this->new_members[this->new_num_members++] = i;
            }
        }
    }
}

// reference: checkBreak in hermite_integrator.h (SDAR)
bool Group::CheckBreak() {

    const Float kappa_org_crit = 1e-2;
    const int n_member = sym_int.particles.getSize();

    // generate binary tree
    sym_int.info.generateBinaryTree(sym_int.particles, manager.interaction.gravitational_constant);

    auto& bin_root = sym_int.info.getBinaryTreeRoot();
    bool outgoing_flag = false; // Indicate whether it is a outgoing case or income case

    // check whether periapsis distance >  2e-3 pc
    // Periapsis is too far and binary is not close enough to use regularized technique.
    // if (bin_root.semi*(1-bin_root.ecc) > 2e-3/position_unit){ // test7
    if (bin_root.semi*(1-bin_root.ecc) > r_search && bin_root.r > 2 * r_search) { // test8 // fiducial
    // if (bin_root.semi*(1-bin_root.ecc) > 1.2e-3/position_unit){ // test12
        fprintf(worker_output_file, "Break group: too far periapsis! (CM PID: %d)\n\t", groupCM->pid);
        fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
        fprintf(worker_output_file, "N_member: %d\n\t", n_member);
        fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit);
        fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
        fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
        fprintf(worker_output_file, "ecca: %e\n\t", bin_root.ecca);
        fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
        fprintf(worker_output_file, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
        fprintf(worker_output_file, "r_crit: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
        /*
        if (n_member > 2) {
            for (int i=0; i<n_member; i++) {
                Particle* ptcl1 = &particles[groupCM->members[i]];
                if (ptcl1->pid == bin_root.getLeftMember()->pid || ptcl1->pid == bin_root.getRightMember()->pid)
                    ptcl1->set_binary_interrupt_state(BinaryInterruptState::none);
                else
                    ptcl1->set_binary_interrupt_state(BinaryInterruptState::manybody);
            }
        }
        */
        /*
        if (n_member == 3) {
            fprintf(worker_output_file, "Left PID: %d, Right PID: %d\n", bin_root.getLeftMember()->pid, bin_root.getRightMember()->pid);
            int outgoingPID = bin_root.getLeftMember()->pid != -1 ? bin_root.getLeftMember()->pid : bin_root.getRightMember()->pid;
            for (int i=0; i<n_member; i++) {
                Particle* ptcl1 = &particles[groupCM->members[i]];
                if (ptcl1->pid == outgoingPID)
                    continue;
                else {
                    fprintf(worker_output_file, "ptcl1 PID: %d\n", ptcl1->pid);
                    ptcl1->new_num_members = 0;
                    ptcl1->set_binary_interrupt_state(BinaryInterruptState::threebody);
                    for (int j=0; j<n_member; j++) {
                        Particle* ptcl2 = &particles[groupCM->members[j]];
                        if (ptcl2->pid != ptcl1->pid && ptcl2->pid != outgoingPID) {
                            ptcl2->set_binary_interrupt_state(BinaryInterruptState::threebody);
                            ptcl2->new_num_members = 0;
                            ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                            fprintf(worker_output_file, "ptcl2 PID: %d, ptcl1 NewNumberOfMember: %d\n", ptcl2->pid, ptcl1->new_num_members);
                        }
                    }
                    break;
                }
            }
        }
        */
        for (int i=0; i<n_member; i++) {
            Particle* ptcl1 = &particles[groupCM->members[i]];
            ptcl1->new_num_members = 0;
        }
        if (sym_int.particles.getSize() == 3) {
            for (int k=0; k<2; k++) {
                if (bin_root.isMemberTree(k)) {
                    auto memberTree = bin_root.getMemberAsTree(k);
                    Particle* ptcl1 = &particles[memberTree->getLeftMember()->particle_index];
                    Particle* ptcl2 = &particles[memberTree->getRightMember()->particle_index];
                    ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                }
            }
        }
        if (sym_int.particles.getSize() == 4) {
            if (bin_root.isMemberTree(0) && bin_root.isMemberTree(1)) {
                auto memberTree1 = bin_root.getMemberAsTree(0);
                auto memberTree2 = bin_root.getMemberAsTree(1);
                Particle* ptcl1 = &particles[memberTree1->getLeftMember()->particle_index];
                Particle* ptcl2 = &particles[memberTree1->getRightMember()->particle_index];
                Particle* ptcl3 = &particles[memberTree2->getLeftMember()->particle_index];
                Particle* ptcl4 = &particles[memberTree2->getRightMember()->particle_index];
                ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                ptcl3->new_members[ptcl3->new_num_members++] = ptcl4->particle_index;
            } else {
                int outgoingPID = bin_root.getLeftMember()->pid != -1 ? bin_root.getLeftMember()->pid : bin_root.getRightMember()->pid;
                for (int i=0; i<4; i++) {
                    Particle* ptcl1 = &particles[groupCM->members[i]];
                    if (ptcl1->pid == outgoingPID)
                        continue;
                    else {
                        for (int j=0; j<4; j++) {
                            Particle* ptcl2 = &particles[groupCM->members[j]];
                            if (ptcl2->pid != ptcl1->pid && ptcl2->pid != outgoingPID) {
                                ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                            }
                        }
                        break;
                    }
                }
            }
        }
        fflush(worker_output_file);
        return true;
    }

    // check binary case 
    // ecc anomaly indicates outgoing (ecca>0) or income (ecca<0)
    if (bin_root.semi>0.0 && bin_root.ecca>0.0) {
        outgoing_flag = true;
        // check whether separation is larger than distance criterion. 
        // /*
        if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 2 * r_search) { // test8 // fiducial
        // if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 1.2e-3/position_unit) { // test12
            fprintf(worker_output_file, "Break group: binary escape! (CM PID: %d)\n\t", groupCM->pid);
            fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
            fprintf(worker_output_file, "N_member: %d\n\t", n_member);
            fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit);
            fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
            fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
            fprintf(worker_output_file, "ecca: %e\n\t", bin_root.ecca);
            fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
            fprintf(worker_output_file, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
            fprintf(worker_output_file, "r_crit: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
            /*
            if (n_member > 2) {
                for (int i=0; i<n_member; i++) {
                    Particle* ptcl1 = &particles[groupCM->members[i]];
                    if (ptcl1->pid == bin_root.getLeftMember()->pid || ptcl1->pid == bin_root.getRightMember()->pid)
                        ptcl1->set_binary_interrupt_state(BinaryInterruptState::none);
                    else
                        ptcl1->set_binary_interrupt_state(BinaryInterruptState::manybody);
                }
            }
            */
            /*
            if (n_member == 3) {
                fprintf(worker_output_file, "Left PID: %d, Right PID: %d\n", bin_root.getLeftMember()->pid, bin_root.getRightMember()->pid);
                int outgoingPID = bin_root.getLeftMember()->pid != -1 ? bin_root.getLeftMember()->pid : bin_root.getRightMember()->pid;
                for (int i=0; i<n_member; i++) {
                    Particle* ptcl1 = &particles[groupCM->members[i]];
                    if (ptcl1->pid == outgoingPID)
                        continue;
                    else {
                        fprintf(worker_output_file, "ptcl1 PID: %d\n", ptcl1->pid);
                        ptcl1->new_num_members = 0;
                        ptcl1->set_binary_interrupt_state(BinaryInterruptState::threebody);
                        for (int j=0; j<n_member; j++) {
                            Particle* ptcl2 = &particles[groupCM->members[j]];
                            if (ptcl2->pid != ptcl1->pid && ptcl2->pid != outgoingPID) {
                                ptcl2->set_binary_interrupt_state(BinaryInterruptState::threebody);
                                ptcl2->new_num_members = 0;
                                ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                                fprintf(worker_output_file, "ptcl2 PID: %d, ptcl1 NewNumberOfMember: %d\n", ptcl2->pid, ptcl1->new_num_members);
                            }
                        }
                        break;
                    }
                }
            }
            */
            for (int i=0; i<n_member; i++) {
                Particle* ptcl1 = &particles[groupCM->members[i]];
                ptcl1->new_num_members = 0;
            }
            if (sym_int.particles.getSize() == 3) {
                for (int k=0; k<2; k++) {
                    if (bin_root.isMemberTree(k)) {
                        auto memberTree = bin_root.getMemberAsTree(k);
                        Particle* ptcl1 = &particles[memberTree->getLeftMember()->particle_index];
                        Particle* ptcl2 = &particles[memberTree->getRightMember()->particle_index];
                        ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                    }
                }
            }
            if (sym_int.particles.getSize() == 4) {
                if (bin_root.isMemberTree(0) && bin_root.isMemberTree(1)) {
                    auto memberTree1 = bin_root.getMemberAsTree(0);
                    auto memberTree2 = bin_root.getMemberAsTree(1);
                    Particle* ptcl1 = &particles[memberTree1->getLeftMember()->particle_index];
                    Particle* ptcl2 = &particles[memberTree1->getRightMember()->particle_index];
                    Particle* ptcl3 = &particles[memberTree2->getLeftMember()->particle_index];
                    Particle* ptcl4 = &particles[memberTree2->getRightMember()->particle_index];
                    ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                    ptcl3->new_members[ptcl3->new_num_members++] = ptcl4->particle_index;
                } else {
                    int outgoingPID = bin_root.getLeftMember()->pid != -1 ? bin_root.getLeftMember()->pid : bin_root.getRightMember()->pid;
                    for (int i=0; i<4; i++) {
                        Particle* ptcl1 = &particles[groupCM->members[i]];
                        if (ptcl1->pid == outgoingPID)
                            continue;
                        else {
                            for (int j=0; j<4; j++) {
                                Particle* ptcl2 = &particles[groupCM->members[j]];
                                if (ptcl2->pid != ptcl1->pid && ptcl2->pid != outgoingPID) {
                                    ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                                }
                            }
                            break;
                        }
                    }
                }
            }
            fflush(worker_output_file);
            return true;
        }
        // */
        /*
        if (n_member == 2) {
            if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 2e-3/position_unit) {
                fprintf(worker_output_file, "Break group: binary escape!\n\t");
                fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
                fprintf(worker_output_file, "N_member: %d\n\t", n_member);
                fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
                fprintf(worker_output_file, "ecca: %e\n\t", bin_root.ecca);
                fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(worker_output_file, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(worker_output_file, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(worker_output_file);
                return true;
            }
        }
        else {
            if (bin_root.r > 2e-3/position_unit) {
                fprintf(worker_output_file, "Break group: binary escape!\n\t");
                fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
                fprintf(worker_output_file, "N_member: %d\n\t", n_member);
                fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
                fprintf(worker_output_file, "ecca: %e\n\t", bin_root.ecca);
                fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(worker_output_file, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(worker_output_file, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(worker_output_file);
                return true;
            }
        }
        */
    }

    // check hyperbolic case
    if (bin_root.semi<0.0) {
        // hyperbolic case, ecca is not correctly calculated
        Float dr2, drdv;
        sym_int.info.getDrDv(dr2, drdv, *bin_root.getLeftMember(), *bin_root.getRightMember());
        if (drdv>0.0) {
            outgoing_flag = true;
            // check distance criterion
            if (bin_root.r > 2 * r_search) { // test8 // seems good! // fiducial
            // if (bin_root.r > 1.2e-3/position_unit) { // test12
                fprintf(worker_output_file, "Break group: hyperbolic escape! (CM PID: %d)\n\t", groupCM->pid);
                fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
                fprintf(worker_output_file, "N_member: %d\n\t", n_member);
                fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit); 
                fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
                fprintf(worker_output_file, "ecca: %e\n\t", bin_root.ecca);
                fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(worker_output_file, "r_crit: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
                /*
                if (n_member > 2) {
                    for (int i=0; i<n_member; i++) {
                        Particle* ptcl1 = &particles[groupCM->members[i]];
                        if (ptcl1->pid == bin_root.getLeftMember()->pid || ptcl1->pid == bin_root.getRightMember()->pid)
                            ptcl1->set_binary_interrupt_state(BinaryInterruptState::none);
                        else
                            ptcl1->set_binary_interrupt_state(BinaryInterruptState::manybody);
                    }
                }
                */
                /*
                if (n_member == 3) {
                    fprintf(worker_output_file, "Left PID: %d, Right PID: %d\n", bin_root.getLeftMember()->pid, bin_root.getRightMember()->pid);
                    int outgoingPID = bin_root.getLeftMember()->pid != -1 ? bin_root.getLeftMember()->pid : bin_root.getRightMember()->pid;
                    for (int i=0; i<n_member; i++) {
                        Particle* ptcl1 = &particles[groupCM->members[i]];
                        if (ptcl1->pid == outgoingPID)
                            continue;
                        else {
                            fprintf(worker_output_file, "ptcl1 PID: %d\n", ptcl1->pid);
                            ptcl1->new_num_members = 0;
                            ptcl1->set_binary_interrupt_state(BinaryInterruptState::threebody);
                            for (int j=0; j<n_member; j++) {
                                Particle* ptcl2 = &particles[groupCM->members[j]];
                                if (ptcl2->pid != ptcl1->pid && ptcl2->pid != outgoingPID) {
                                    ptcl2->set_binary_interrupt_state(BinaryInterruptState::threebody);
                                    ptcl2->new_num_members = 0;
                                    ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                                    fprintf(worker_output_file, "ptcl2 PID: %d, ptcl1 NewNumberOfMember: %d\n", ptcl2->pid, ptcl1->new_num_members);
                                }
                            }
                            break;
                        }
                    }
                }
                */
                for (int i=0; i<n_member; i++) {
                    Particle* ptcl1 = &particles[groupCM->members[i]];
                    ptcl1->new_num_members = 0;
                }
                if (sym_int.particles.getSize() == 3) {
                    for (int k=0; k<2; k++) {
                        if (bin_root.isMemberTree(k)) {
                            auto memberTree = bin_root.getMemberAsTree(k);
                            Particle* ptcl1 = &particles[memberTree->getLeftMember()->particle_index];
                            Particle* ptcl2 = &particles[memberTree->getRightMember()->particle_index];
                            ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                        }
                    }
                }
                if (sym_int.particles.getSize() == 4) {
                    if (bin_root.isMemberTree(0) && bin_root.isMemberTree(1)) {
                        auto memberTree1 = bin_root.getMemberAsTree(0);
                        auto memberTree2 = bin_root.getMemberAsTree(1);
                        Particle* ptcl1 = &particles[memberTree1->getLeftMember()->particle_index];
                        Particle* ptcl2 = &particles[memberTree1->getRightMember()->particle_index];
                        Particle* ptcl3 = &particles[memberTree2->getLeftMember()->particle_index];
                        Particle* ptcl4 = &particles[memberTree2->getRightMember()->particle_index];
                        ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                        ptcl3->new_members[ptcl3->new_num_members++] = ptcl4->particle_index;
                    } else {
                        int outgoingPID = bin_root.getLeftMember()->pid != -1 ? bin_root.getLeftMember()->pid : bin_root.getRightMember()->pid;
                        for (int i=0; i<4; i++) {
                            Particle* ptcl1 = &particles[groupCM->members[i]];
                            if (ptcl1->pid == outgoingPID)
                                continue;
                            else {
                                for (int j=0; j<4; j++) {
                                    Particle* ptcl2 = &particles[groupCM->members[j]];
                                    if (ptcl2->pid != ptcl1->pid && ptcl2->pid != outgoingPID) {
                                        ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                fflush(worker_output_file);
                return true;
            }
            /*
            if (n_member == 2) {
                if (bin_root.r > 1e-3/position_unit) { // original
                    // if (bin_root.r > 2e-3/position_unit) { // test
                    fprintf(worker_output_file, "Break group: hyperbolic escape!\n\t");
                    fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
                    fprintf(worker_output_file, "N_member: %d\n\t", n_member);
                    fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit); 
                    fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
                    fprintf(worker_output_file, "ecca: %e\n\t", bin_root.ecca);
                    fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                    fprintf(worker_output_file, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                    fflush(worker_output_file);
                    return true;
                }
            }
            else {
                if (bin_root.r > 2e-3/position_unit) { // original
                    // if (bin_root.r > 2e-3/position_unit) { // test
                    fprintf(worker_output_file, "Break group: hyperbolic escape!\n\t");
                    fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
                    fprintf(worker_output_file, "N_member: %d\n\t", n_member);
                    fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit); 
                    fprintf(worker_output_file, "ecc: %e\n\t", bin_root.ecc);
                    fprintf(worker_output_file, "ecca: %e\n\t", bin_root.ecca);
                    fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                    fprintf(worker_output_file, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                    fflush(worker_output_file);
                    return true;
                }
            }
            */
        }

    }
    // return false;
/* // test_1e4_2
    // check strong perturbed binary case (only check further if it is outgoing case)
    // calculate slowdown in a consistent way like in checknewgroup to avoid switching
    // fcm may not properly represent the perturbation force (perturber mass is unknown)
    if (outgoing_flag && bin_root.semi > 0.0 && n_member == 2) {

        AR::SlowDown sd;
        auto& sd_group = sym_int.info.getBinaryTreeRoot().slowdown;
        sd.initialSlowDownReference(sd_group.getSlowDownFactorReference(),sd_group.getSlowDownFactorMax());
        sd.timescale = sd_group.timescale;
        sd.period = sd_group.period;

        // sd.pert_in = manager.interaction.calcPertFromBinary(bin_root); // original
        sd.pert_in = manager.interaction.calcPertFromMR(bin_root.r, bin_root.m1, bin_root.m2);
        Float acc_cm[3];
        for (int i = 0; i < DIM; ++i) {
            acc_cm[i] = sym_int.particles.cm.acc_total[i][0]; // a_irr or a_tot? // I think a_tot is absolutely right by EW 2025.7.19
        }
        Float fcm[3] = {acc_cm[0]*bin_root.mass, acc_cm[1]*bin_root.mass, acc_cm[2]*bin_root.mass};
        sd.pert_out= manager.interaction.calcPertFromForce(fcm, bin_root.mass, bin_root.mass);
        sd.calcSlowDownFactor();
        Float kappa_org = sd.getSlowDownFactorOrigin();

        if (kappa_org<kappa_org_crit) {
            // in binary case, only break when apo is larger than distance criterion
            Float apo = bin_root.semi * (1.0 + bin_root.ecc);
            // if (apo>sym_int.info.r_break_crit||bin_root.semi<0.0) {
            if (apo>sym_int.info.r_break_crit && bin_root.r > 2 * r_search) { // test8 // fiducial
            // if ((apo>sym_int.info.r_break_crit && bin_root.r > 1.2e-3/position_unit)||bin_root.semi<0.0) { // test12
                auto& sd_root = sym_int.info.getBinaryTreeRoot().slowdown;

                fprintf(worker_output_file, "Break group: strong perturbed! (CM PID: %d)\n\t", groupCM->pid);
                fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
                fprintf(worker_output_file, "N_member: %d\n\t", n_member);
                fprintf(worker_output_file, "pert_in: %e \n\t", sd_root.pert_in);
                fprintf(worker_output_file, "pert_out: %e \n\t", sd_root.pert_out);
                fprintf(worker_output_file, "kappa_org: %e \n\t", kappa_org);
                fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(worker_output_file, "ecc: %e \n\t", bin_root.ecc);
                fprintf(worker_output_file, "ecca: %e \n\t", bin_root.ecca);
                fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(worker_output_file, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(worker_output_file, "r_break: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
                fflush(worker_output_file);
                return true;
            }
        }
    }
*/ // test_1e4_2
    return false;

}

bool Group::CheckBreak2() {

    sym_int.info.generateBinaryTree(sym_int.particles, manager.interaction.gravitational_constant);
    auto& bin_root = sym_int.info.getBinaryTreeRoot();

    if (bin_root.r > r_search) {

        if (bin_root.semi > 0.0 && bin_root.ecca < 0.0) // incoming binary
            return false;

        fprintf(worker_output_file, "Break group: r > RSEARCH! (CM PID: %d)\n\t", groupCM->pid);
        fprintf(worker_output_file, "time: %e Myr\n\t", CurrentTime*enzo_time_step*1e4);
        fprintf(worker_output_file, "N_member: %d\n\t", sym_int.particles.getSize());
        fprintf(worker_output_file, "separation: %e pc\n\t", bin_root.r*position_unit);
        fprintf(worker_output_file, "semi: %e pc\n\t", bin_root.semi*position_unit);
        fprintf(worker_output_file, "ecc: %e \n\t", bin_root.ecc);
        fprintf(worker_output_file, "ecca: %e \n\t", bin_root.ecca);
        fprintf(worker_output_file, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
        fprintf(worker_output_file, "apo: %e pc\n\n", bin_root.semi*(1+bin_root.ecc)*position_unit);
        fflush(worker_output_file);

        for (int i=0; i<sym_int.particles.getSize(); i++) {
            Particle* ptcl1 = &particles[groupCM->members[i]];
            ptcl1->new_num_members = 0;
        }

        if (sym_int.particles.getSize() == 3) {
            for (int k=0; k<2; k++) {
                if (bin_root.isMemberTree(k)) {
                    auto memberTree = bin_root.getMemberAsTree(k);
                    Particle* ptcl1 = &particles[memberTree->getLeftMember()->particle_index];
                    Particle* ptcl2 = &particles[memberTree->getRightMember()->particle_index];
                    ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                }
            }
        }
        if (sym_int.particles.getSize() == 4) {
            if (bin_root.isMemberTree(0) && bin_root.isMemberTree(1)) {
                auto memberTree1 = bin_root.getMemberAsTree(0);
                auto memberTree2 = bin_root.getMemberAsTree(1);
                Particle* ptcl1 = &particles[memberTree1->getLeftMember()->particle_index];
                Particle* ptcl2 = &particles[memberTree1->getRightMember()->particle_index];
                Particle* ptcl3 = &particles[memberTree2->getLeftMember()->particle_index];
                Particle* ptcl4 = &particles[memberTree2->getRightMember()->particle_index];
                ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                ptcl3->new_members[ptcl3->new_num_members++] = ptcl4->particle_index;
            } else {
                int outgoingPID = bin_root.getLeftMember()->pid != -1 ? bin_root.getLeftMember()->pid : bin_root.getRightMember()->pid;
                for (int i=0; i<4; i++) {
                    Particle* ptcl1 = &particles[groupCM->members[i]];
                    if (ptcl1->pid == outgoingPID)
                        continue;
                    else {
                        for (int j=0; j<4; j++) {
                            Particle* ptcl2 = &particles[groupCM->members[j]];
                            if (ptcl2->pid != ptcl1->pid && ptcl2->pid != outgoingPID) {
                                ptcl1->new_members[ptcl1->new_num_members++] = ptcl2->particle_index;
                            }
                        }
                        break;
                    }
                }
            }
        }
        return true;
    }
    return false;
}
#endif