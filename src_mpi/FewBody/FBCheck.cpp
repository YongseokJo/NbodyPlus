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
void Particle::checkNewGroup() {

    const Float kappa_org_crit = 1e-2; // kappa_org criterion for new group kappa_org>kappa_org_crit
    const double r_crit = RSEARCH/position_unit; // distance criterion

    double pos1[Dim], vel1[Dim];

    Particle* ptcl2;

    std::unordered_set<int> CMPtclsSet;

    // check only active particles 
    // single case
    for (int i=0; i < this->NumberOfNeighbor; i++) {
        ptcl2 = &particles[this->Neighbors[i]];
        if (!ptcl2->isActive) {
            if (ptcl2->CMPtclIndex != -1) {
                CMPtclsSet.insert(ptcl2->CMPtclIndex);
            }
            continue;
        }

        // if (ptcl2->TimeStepIrr > this->TimeStepIrr) // test_1e5_4 & 5: this must make the same result!
        if (ptcl2->TimeStepIrr*EnzoTimeStep*1e4 > TSEARCH) // fiducial: 1e-5 but for RSEARCH = 0.00025 pc, 1e-6 Myr seems good
            continue;

        double dt = this->CurrentTimeIrr > ptcl2->CurrentTimeIrr ? \
                        this->CurrentTimeIrr - ptcl2->CurrentTimeIrr : ptcl2->CurrentTimeIrr - this->CurrentTimeIrr;

        double pos2[Dim], vel2[Dim];
        
		this->predictParticleSecondOrder(dt, pos1, vel1);
		ptcl2->predictParticleSecondOrder(dt, pos2, vel2);

        const Float dr = dist(pos1, pos2);
        
        if (dr < r_crit) {

            Float drdv = calcDrDv(pos1, pos2, vel1, vel2);
            // only inwards
            if(drdv<0.0) {
// /* // test_1e4_2
                Float fcm[3] = {this->Mass*this->a_irr[0][0] + ptcl2->Mass*ptcl2->a_irr[0][0], 
                                this->Mass*this->a_irr[1][0] + ptcl2->Mass*ptcl2->a_irr[1][0], 
                                this->Mass*this->a_irr[2][0] + ptcl2->Mass*ptcl2->a_irr[2][0]};

                AR::SlowDown sd;
                Interaction interaction;
                Float mcm = this->Mass + ptcl2->Mass;

                // sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
                sd.initialSlowDownReference(1e-6, NUMERIC_FLOAT_MAX);

                sd.pert_in = interaction.calcPertFromMR(dr, this->Mass, ptcl2->Mass);
                sd.pert_out = interaction.calcPertFromForce(fcm, mcm, mcm);

                sd.calcSlowDownFactor();
                Float kappa_org = sd.getSlowDownFactorOrigin();

                // avoid strong perturbed case, estimate perturbation
                // if kappa_org < criterion, avoid to form new group, should be consistent as checkbreak
                if(kappa_org<kappa_org_crit) continue;
// */ // test_1e4_2
                this->NewNeighbors[this->NewNumberOfNeighbor] = this->Neighbors[i];
                this->NewNumberOfNeighbor++;
            }
        }
    }
    if (!CMPtclsSet.empty()) {
        for (int i: CMPtclsSet) {
			ptcl2 = &particles[i];

            if (this->PID == ptcl2->PID) {
				continue;
			}

			if (!ptcl2->isActive) {
				fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", this->PID, ptcl2->PID);
				assert(ptcl2->isActive);
			}

            // if (ptcl2->TimeStepIrr > this->TimeStepIrr) // test_1e5_4 & 5: this must make the same result!
            if (ptcl2->TimeStepIrr*EnzoTimeStep*1e4 > TSEARCH) // fiducial: 1e-5 but for RSEARCH = 0.00025 pc, 1e-6 Myr seems good
                continue;

            double dt = this->CurrentTimeIrr > ptcl2->CurrentTimeIrr ? \
                            this->CurrentTimeIrr - ptcl2->CurrentTimeIrr : ptcl2->CurrentTimeIrr - this->CurrentTimeIrr;

            double pos2[Dim], vel2[Dim];
            
            this->predictParticleSecondOrder(dt, pos1, vel1);
            ptcl2->predictParticleSecondOrder(dt, pos2, vel2);

            const Float dr = dist(pos1, pos2);
            
            if (dr < r_crit) {

                Float drdv = calcDrDv(pos1, pos2, vel1, vel2);
                // only inwards
                if(drdv<0.0) {
// /* // test_1e4_2
                    Float fcm[3] = {this->Mass*this->a_irr[0][0] + ptcl2->Mass*ptcl2->a_irr[0][0], 
                                    this->Mass*this->a_irr[1][0] + ptcl2->Mass*ptcl2->a_irr[1][0], 
                                    this->Mass*this->a_irr[2][0] + ptcl2->Mass*ptcl2->a_irr[2][0]};

                    AR::SlowDown sd;
                    Interaction interaction;
                    Float mcm = this->Mass + ptcl2->Mass;

                    // sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
                    sd.initialSlowDownReference(1e-6, NUMERIC_FLOAT_MAX);

                    sd.pert_in = interaction.calcPertFromMR(dr, this->Mass, ptcl2->Mass);
                    sd.pert_out = interaction.calcPertFromForce(fcm, mcm, mcm);

                    sd.calcSlowDownFactor();
                    Float kappa_org = sd.getSlowDownFactorOrigin();

                    // avoid strong perturbed case, estimate perturbation
                    // if kappa_org < criterion, avoid to form new group, should be consistent as checkbreak
                    if(kappa_org<kappa_org_crit) continue;
// */ // test_1e4_2
                    this->NewNeighbors[this->NewNumberOfNeighbor] = i;
                    this->NewNumberOfNeighbor++;
                }
            }
        }
    }
}

// Use when many-body (>3) group broke
// Or primordial binary search
void Particle::checkNewGroup2() {
    
    const Float kappa_org_crit = 1e-2; // kappa_org criterion for new group kappa_org>kappa_org_crit
    const double r_crit = RSEARCH/position_unit; // distance criterion

    double pos1[Dim], vel1[Dim];

    Particle* ptcl2;

    std::unordered_set<int> CMPtclsSet;

    // check only active particles 
    // single case
    for (int i=0; i < this->NumberOfNeighbor; i++) {
        ptcl2 = &particles[this->Neighbors[i]];
        if (!ptcl2->isActive) {
            if (ptcl2->CMPtclIndex != -1) {
                CMPtclsSet.insert(ptcl2->CMPtclIndex);
            }
            continue;
        }

        double dt = this->CurrentTimeIrr > ptcl2->CurrentTimeIrr ? \
                        this->CurrentTimeIrr - ptcl2->CurrentTimeIrr : ptcl2->CurrentTimeIrr - this->CurrentTimeIrr;

        double pos2[Dim], vel2[Dim];
        
		this->predictParticleSecondOrder(dt, pos1, vel1);
		ptcl2->predictParticleSecondOrder(dt, pos2, vel2);

        const Float dr = dist(pos1, pos2);

        double v2 = std::pow(dist(vel1, vel2), 2);
        double energy = v2/2 - (this->Mass + ptcl2->Mass)/dr; // determine they are bound or not
        
        if (dr < r_crit && energy < 0) {
            this->NewNeighbors[this->NewNumberOfNeighbor] = this->Neighbors[i];
            this->NewNumberOfNeighbor++;
        }
    }
    if (!CMPtclsSet.empty()) {
        for (int i: CMPtclsSet) {
			ptcl2 = &particles[i];

            if (this->PID == ptcl2->PID) {
				continue;
			}

			if (!ptcl2->isActive) {
				fprintf(stderr, "Why inactive CM ptcl? this PID: %d, neighbor PID: %d\n", this->PID, ptcl2->PID);
				assert(ptcl2->isActive);
			}

            double dt = this->CurrentTimeIrr > ptcl2->CurrentTimeIrr ? \
                            this->CurrentTimeIrr - ptcl2->CurrentTimeIrr : ptcl2->CurrentTimeIrr - this->CurrentTimeIrr;

            double pos2[Dim], vel2[Dim];
            
            this->predictParticleSecondOrder(dt, pos1, vel1);
            ptcl2->predictParticleSecondOrder(dt, pos2, vel2);

            const Float dr = dist(pos1, pos2);

            double v2 = std::pow(dist(vel1, vel2), 2);
            double energy = v2/2 - (this->Mass + ptcl2->Mass)/dr; // determine they are bound or not
            
            if (dr < r_crit && energy < 0) {
                this->NewNeighbors[this->NewNumberOfNeighbor] = i;
                this->NewNumberOfNeighbor++;
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
    if (bin_root.semi*(1-bin_root.ecc) > RSEARCH/position_unit && bin_root.r > 2 * RSEARCH/position_unit) { // test8 // fiducial
    // if (bin_root.semi*(1-bin_root.ecc) > 1.2e-3/position_unit){ // test12
        fprintf(workerout, "Break group: too far periapsis! (CM PID: %d)\n\t", groupCM->PID);
        fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
        fprintf(workerout, "N_member: %d\n\t", n_member);
        fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit);
        fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
        fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
        fprintf(workerout, "ecca: %e\n\t", bin_root.ecca);
        fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
        fprintf(workerout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
        fprintf(workerout, "r_crit: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
        /*
        if (n_member > 2) {
            for (int i=0; i<n_member; i++) {
                Particle* ptcl1 = &particles[groupCM->Members[i]];
                if (ptcl1->PID == bin_root.getLeftMember()->PID || ptcl1->PID == bin_root.getRightMember()->PID)
                    ptcl1->setBinaryInterruptState(BinaryInterruptState::none);
                else
                    ptcl1->setBinaryInterruptState(BinaryInterruptState::manybody);
            }
        }
        */
        // /*
        if (n_member == 3) {
            fprintf(workerout, "Left PID: %d, Right PID: %d\n", bin_root.getLeftMember()->PID, bin_root.getRightMember()->PID);
            int outgoingPID = bin_root.getLeftMember()->PID != -1 ? bin_root.getLeftMember()->PID : bin_root.getRightMember()->PID;
            for (int i=0; i<n_member; i++) {
                Particle* ptcl1 = &particles[groupCM->Members[i]];
                if (ptcl1->PID == outgoingPID)
                    continue;
                else {
                    fprintf(workerout, "ptcl1 PID: %d\n", ptcl1->PID);
                    ptcl1->NewNumberOfNeighbor = 0;
                    ptcl1->setBinaryInterruptState(BinaryInterruptState::threebody);
                    for (int j=0; j<n_member; j++) {
                        Particle* ptcl2 = &particles[groupCM->Members[j]];
                        if (ptcl2->PID != ptcl1->PID && ptcl2->PID != outgoingPID) {
                            ptcl2->setBinaryInterruptState(BinaryInterruptState::threebody);
                            ptcl2->NewNumberOfNeighbor = 0;
                            ptcl1->NewNeighbors[ptcl1->NewNumberOfNeighbor++] = ptcl2->ParticleIndex;
                            fprintf(workerout, "ptcl2 PID: %d, ptcl1 NewNumberOfNeighbor: %d\n", ptcl2->PID, ptcl1->NewNumberOfNeighbor);
                        }
                    }
                    break;
                }
            }
        }
        // */
        fflush(workerout);
        return true;
    }

    // check binary case 
    // ecc anomaly indicates outgoing (ecca>0) or income (ecca<0)
    if (bin_root.semi>0.0 && bin_root.ecca>0.0) {
        outgoing_flag = true;
        // check whether separation is larger than distance criterion. 
        // /*
        if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 2 * RSEARCH/position_unit) { // test8 // fiducial
        // if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 1.2e-3/position_unit) { // test12
            fprintf(workerout, "Break group: binary escape! (CM PID: %d)\n\t", groupCM->PID);
            fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
            fprintf(workerout, "N_member: %d\n\t", n_member);
            fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit);
            fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
            fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
            fprintf(workerout, "ecca: %e\n\t", bin_root.ecca);
            fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
            fprintf(workerout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
            fprintf(workerout, "r_crit: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
            /*
            if (n_member > 2) {
                for (int i=0; i<n_member; i++) {
                    Particle* ptcl1 = &particles[groupCM->Members[i]];
                    if (ptcl1->PID == bin_root.getLeftMember()->PID || ptcl1->PID == bin_root.getRightMember()->PID)
                        ptcl1->setBinaryInterruptState(BinaryInterruptState::none);
                    else
                        ptcl1->setBinaryInterruptState(BinaryInterruptState::manybody);
                }
            }
            */
            // /*
            if (n_member == 3) {
                fprintf(workerout, "Left PID: %d, Right PID: %d\n", bin_root.getLeftMember()->PID, bin_root.getRightMember()->PID);
                int outgoingPID = bin_root.getLeftMember()->PID != -1 ? bin_root.getLeftMember()->PID : bin_root.getRightMember()->PID;
                for (int i=0; i<n_member; i++) {
                    Particle* ptcl1 = &particles[groupCM->Members[i]];
                    if (ptcl1->PID == outgoingPID)
                        continue;
                    else {
                        fprintf(workerout, "ptcl1 PID: %d\n", ptcl1->PID);
                        ptcl1->NewNumberOfNeighbor = 0;
                        ptcl1->setBinaryInterruptState(BinaryInterruptState::threebody);
                        for (int j=0; j<n_member; j++) {
                            Particle* ptcl2 = &particles[groupCM->Members[j]];
                            if (ptcl2->PID != ptcl1->PID && ptcl2->PID != outgoingPID) {
                                ptcl2->setBinaryInterruptState(BinaryInterruptState::threebody);
                                ptcl2->NewNumberOfNeighbor = 0;
                                ptcl1->NewNeighbors[ptcl1->NewNumberOfNeighbor++] = ptcl2->ParticleIndex;
                                fprintf(workerout, "ptcl2 PID: %d, ptcl1 NewNumberOfNeighbor: %d\n", ptcl2->PID, ptcl1->NewNumberOfNeighbor);
                            }
                        }
                        break;
                    }
                }
            }
            // */
            fflush(workerout);
            return true;
        }
        // */
        /*
        if (n_member == 2) {
            if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 2e-3/position_unit) {
                fprintf(workerout, "Break group: binary escape!\n\t");
                fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(workerout, "N_member: %d\n\t", n_member);
                fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
                fprintf(workerout, "ecca: %e\n\t", bin_root.ecca);
                fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(workerout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(workerout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(workerout);
                return true;
            }
        }
        else {
            if (bin_root.r > 2e-3/position_unit) {
                fprintf(workerout, "Break group: binary escape!\n\t");
                fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(workerout, "N_member: %d\n\t", n_member);
                fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
                fprintf(workerout, "ecca: %e\n\t", bin_root.ecca);
                fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(workerout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(workerout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(workerout);
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
            if (bin_root.r > 2 * RSEARCH/position_unit) { // test8 // seems good! // fiducial
            // if (bin_root.r > 1.2e-3/position_unit) { // test12
                fprintf(workerout, "Break group: hyperbolic escape! (CM PID: %d)\n\t", groupCM->PID);
                fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(workerout, "N_member: %d\n\t", n_member);
                fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit); 
                fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
                fprintf(workerout, "ecca: %e\n\t", bin_root.ecca);
                fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(workerout, "r_crit: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
                /*
                if (n_member > 2) {
                    for (int i=0; i<n_member; i++) {
                        Particle* ptcl1 = &particles[groupCM->Members[i]];
                        if (ptcl1->PID == bin_root.getLeftMember()->PID || ptcl1->PID == bin_root.getRightMember()->PID)
                            ptcl1->setBinaryInterruptState(BinaryInterruptState::none);
                        else
                            ptcl1->setBinaryInterruptState(BinaryInterruptState::manybody);
                    }
                }
                */
                // /*
                if (n_member == 3) {
                    fprintf(workerout, "Left PID: %d, Right PID: %d\n", bin_root.getLeftMember()->PID, bin_root.getRightMember()->PID);
                    int outgoingPID = bin_root.getLeftMember()->PID != -1 ? bin_root.getLeftMember()->PID : bin_root.getRightMember()->PID;
                    for (int i=0; i<n_member; i++) {
                        Particle* ptcl1 = &particles[groupCM->Members[i]];
                        if (ptcl1->PID == outgoingPID)
                            continue;
                        else {
                            fprintf(workerout, "ptcl1 PID: %d\n", ptcl1->PID);
                            ptcl1->NewNumberOfNeighbor = 0;
                            ptcl1->setBinaryInterruptState(BinaryInterruptState::threebody);
                            for (int j=0; j<n_member; j++) {
                                Particle* ptcl2 = &particles[groupCM->Members[j]];
                                if (ptcl2->PID != ptcl1->PID && ptcl2->PID != outgoingPID) {
                                    ptcl2->setBinaryInterruptState(BinaryInterruptState::threebody);
                                    ptcl2->NewNumberOfNeighbor = 0;
                                    ptcl1->NewNeighbors[ptcl1->NewNumberOfNeighbor++] = ptcl2->ParticleIndex;
                                    fprintf(workerout, "ptcl2 PID: %d, ptcl1 NewNumberOfNeighbor: %d\n", ptcl2->PID, ptcl1->NewNumberOfNeighbor);
                                }
                            }
                            break;
                        }
                    }
                }
                // */
                fflush(workerout);
                return true;
            }
            /*
            if (n_member == 2) {
                if (bin_root.r > 1e-3/position_unit) { // original
                    // if (bin_root.r > 2e-3/position_unit) { // test
                    fprintf(workerout, "Break group: hyperbolic escape!\n\t");
                    fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                    fprintf(workerout, "N_member: %d\n\t", n_member);
                    fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit); 
                    fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
                    fprintf(workerout, "ecca: %e\n\t", bin_root.ecca);
                    fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                    fprintf(workerout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                    fflush(workerout);
                    return true;
                }
            }
            else {
                if (bin_root.r > 2e-3/position_unit) { // original
                    // if (bin_root.r > 2e-3/position_unit) { // test
                    fprintf(workerout, "Break group: hyperbolic escape!\n\t");
                    fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                    fprintf(workerout, "N_member: %d\n\t", n_member);
                    fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit); 
                    fprintf(workerout, "ecc: %e\n\t", bin_root.ecc);
                    fprintf(workerout, "ecca: %e\n\t", bin_root.ecca);
                    fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                    fprintf(workerout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                    fflush(workerout);
                    return true;
                }
            }
            */
        }

    }
    // return false;
// /* // test_1e4_2
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
        for (int i = 0; i < Dim; ++i) {
            acc_cm[i] = sym_int.particles.cm.a_irr[i][0]; // a_irr or a_tot?
        }
        Float fcm[3] = {acc_cm[0]*bin_root.Mass, acc_cm[1]*bin_root.Mass, acc_cm[2]*bin_root.Mass};
        sd.pert_out= manager.interaction.calcPertFromForce(fcm, bin_root.Mass, bin_root.Mass);
        sd.calcSlowDownFactor();
        Float kappa_org = sd.getSlowDownFactorOrigin();

        if (kappa_org<kappa_org_crit) {
            // in binary case, only break when apo is larger than distance criterion
            Float apo = bin_root.semi * (1.0 + bin_root.ecc);
            // if (apo>sym_int.info.r_break_crit||bin_root.semi<0.0) {
            if (apo>sym_int.info.r_break_crit && bin_root.r > 2 * RSEARCH/position_unit) { // test8 // fiducial
            // if ((apo>sym_int.info.r_break_crit && bin_root.r > 1.2e-3/position_unit)||bin_root.semi<0.0) { // test12
                auto& sd_root = sym_int.info.getBinaryTreeRoot().slowdown;

                fprintf(workerout, "Break group: strong perturbed! (CM PID: %d)\n\t", groupCM->PID);
                fprintf(workerout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(workerout, "N_member: %d\n\t", n_member);
                fprintf(workerout, "pert_in: %e \n\t", sd_root.pert_in);
                fprintf(workerout, "pert_out: %e \n\t", sd_root.pert_out);
                fprintf(workerout, "kappa_org: %e \n\t", kappa_org);
                fprintf(workerout, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(workerout, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(workerout, "ecc: %e \n\t", bin_root.ecc);
                fprintf(workerout, "ecca: %e \n\t", bin_root.ecca);
                fprintf(workerout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(workerout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(workerout, "r_break: %e pc\n\n", sym_int.info.r_break_crit*position_unit);
                fflush(workerout);
                return true;
            }
        }
    }
// */ // test_1e4_2
    return false;

}
#endif