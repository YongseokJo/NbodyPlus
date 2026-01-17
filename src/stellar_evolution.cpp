#ifdef SEVN

#include "global.h"
#include <random>
#include <map>

#ifdef SEVN_BINARY
void UpdateEvolution(Particle* ptcl, bool from_BSE);
void UpdateBinaryEvolution(Particle* ptcl);
void convertBinaryToSingle(Particle* ptclCM);
std::vector<std::string> getCustomInitParams(Particle* ptcl);
#else
void UpdateEvolution(Particle* ptcl);
#endif

// Start SEVN stellar evolution
// Set stellar radii, BH spin, etc
void initializeStellarEvolution() {

    std::ostringstream oss;
	oss << sevnio->svpar;
	fprintf(sevn_output_file, "%s\n\n\n\n\n", oss.str().c_str());
    // fprintf(sevn_output_file, "PID\tPhase\tMass (Msun)\tRadius (pc)\tTime (Myr)\tWorldtime (Myr)\n")
	fflush(sevn_output_file);

    assert(SEVNList.empty());

    for (int i=0; i<num_particles; i++) {

        Particle* ptcl = &particles[i];

        ptcl->particle_type = NO_FEEDBACK_STAR;
        ptcl->formation_time = 0.0;
        ptcl->world_time = 0.0;

        if (ptcl->mass*mass_unit < 2.2) {
			ptcl->radius = 2.25461e-8/position_unit*pow(ptcl->mass*mass_unit, 1/3);
            // stellar radius in code unit
            // extrapolated from the solar radius (it is assumed that low-mass stars have the same stellar density to the sun)

			// fprintf(sevn_output_file, "PID: %d. Mass: %e Msol, Radius: %e pc\n", ptcl->pid, ptcl->mass*mass_unit, ptcl->radius*position_unit);
			continue;
		}

        // Mass, metallicity, spin, sn model, tini, tf, dtout, random seed(optional)
        std::stringstream mass;
        mass << std::setprecision(17) << ptcl->mass*mass_unit; // convert to Msun

		std::vector<std::string> init_params{mass.str(), "0.0002", "0.0", "delayed", "zams", "end", "events"};

        size_t id = ptcl->pid;
        ptcl->stellar_evolution = new StarSEVN(sevnio, init_params, id, false);
        SEVNList.insert({ptcl->world_time + ptcl->stellar_evolution->getp(Timestep::ID), ptcl->particle_index});

		ptcl->radius = ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // stellar radius in code unit
        if (ptcl->stellar_evolution->amiremnant())
            ptcl->particle_type = REMNANT + (int)ptcl->stellar_evolution->getp(RemnantType::ID);
        else
            ptcl->particle_type = (int)ptcl->stellar_evolution->getp(Phase::ID);
    }
}

void setBHspin(Particle* ptcl) {
    std::random_device rd; // Obtain a random number from hardware
    std::mt19937 mt(rd()); // Seed the generator
    std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
    double phi = 2 * M_PI * distr(mt);
    double theta = M_PI * distr(mt);

    ptcl->spin_param[0] = ptcl->stellar_evolution->getp(Xspin::ID) * sin(theta)*cos(phi);
    ptcl->spin_param[1] = ptcl->stellar_evolution->getp(Xspin::ID) * sin(theta)*sin(phi);
    ptcl->spin_param[2] = ptcl->stellar_evolution->getp(Xspin::ID) * cos(theta);
}

void StellarEvolution() {

    Particle* ptcl;
    while (!SEVNList.empty()) {
         
        auto it = SEVNList.begin();
        ptcl = &particles[it->second];
        ptcl->world_time += ptcl->stellar_evolution->getp(Timestep::ID);
        ptcl->stellar_evolution->evolve();

        while (ptcl->world_time + ptcl->stellar_evolution->getp(Timestep::ID) <= global_time * enzo_time_step * 1e4) {
            ptcl->world_time += ptcl->stellar_evolution->getp(Timestep::ID);
            ptcl->stellar_evolution->evolve();
        }

        it = SEVNList.erase(it);
#ifdef SEVN_BINARY
        UpdateEvolution(ptcl, false);
#else
        UpdateEvolution(ptcl);
#endif

        if (ptcl->cm_particle_index != -1) {
            Particle* ptclCM = &particles[ptcl->cm_particle_index];
            double CMPtclMass = 0.0;
            bool kicked = true;
            for (int i=0; i<ptclCM->num_members; i++) {
                Particle* member = &particles[ptclCM->members[i]];
                if (member->mass > 0.0)
                    CMPtclMass += member->mass;
            }
            ptclCM->mass = CMPtclMass;
        }

        if (SEVNList.empty() || SEVNList.begin()->first > global_time * enzo_time_step * 1e4)
            break;
    }
    fflush(sevn_output_file);
}

#ifdef SEVN_BINARY
void UpdateEvolution(Particle* ptcl, bool from_BSE) {
#else
void UpdateEvolution(Particle* ptcl) {
#endif

    if (ptcl->stellar_evolution->amiremnant())
        ptcl->particle_type = REMNANT + (int)ptcl->stellar_evolution->getp(RemnantType::ID);
    else
        ptcl->particle_type = (int)ptcl->stellar_evolution->getp(Phase::ID);

    if (!ptcl->stellar_evolution->amiremnant()) {
#ifdef SEVN_BINARY
        if (!from_BSE)
            SEVNList.insert({ptcl->world_time + ptcl->stellar_evolution->getp(Timestep::ID), ptcl->particle_index});
#else
        SEVNList.insert({ptcl->world_time + ptcl->stellar_evolution->getp(Timestep::ID), ptcl->particle_index});
#endif
        // fprintf(sevn_output_file, "PID: %d, Phase: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, int(ptcl->stellar_evolution->getp(Phase::ID)), ptcl->stellar_evolution->getp(Mass::ID), ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        ptcl->delta_mass += ptcl->mass - ptcl->stellar_evolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->mass = ptcl->stellar_evolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->mass*mass_unit > ptcl->stellar_evolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->mass*mass_unit/ptcl->stellar_evolution->get_max_zams(), 1./3);
        fprintf(sevn_output_file, "PID: %d, Phase: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, int(ptcl->stellar_evolution->getp(Phase::ID)), ptcl->mass*mass_unit, ptcl->stellar_evolution->get_zams(), ptcl->radius*position_unit, ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
    }
    else if (ptcl->stellar_evolution->amiWD()) {
        // SEVNList.erase(ptcl->particle_index);
        // fprintf(sevn_output_file, "WD. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->stellar_evolution->getp(Mass::ID), ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        ptcl->delta_mass += ptcl->mass - ptcl->stellar_evolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->mass = ptcl->stellar_evolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        // ptcl->world_time = NUMERIC_FLOAT_MAX;
        fprintf(sevn_output_file, "WD. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->mass*mass_unit, ptcl->stellar_evolution->get_zams(), ptcl->radius*position_unit, ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        if (ptcl->stellar_evolution->vkick[3] > 0.0) {
            fprintf(sevn_output_file, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->stellar_evolution->vkick[0], ptcl->stellar_evolution->vkick[1], ptcl->stellar_evolution->vkick[2]);
            for(int i=0; i<DIM; i++)
                ptcl->velocity[i] += ptcl->stellar_evolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->cm_particle_index != -1)
                ptcl->set_binary_interrupt_state(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->stellar_evolution->amiNS()) {
        // SEVNList.erase(ptcl->particle_index);
        // fprintf(sevn_output_file, "NS. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->stellar_evolution->getp(Mass::ID), ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        ptcl->delta_mass += ptcl->mass - ptcl->stellar_evolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->mass = ptcl->stellar_evolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        // ptcl->world_time = NUMERIC_FLOAT_MAX;
        fprintf(sevn_output_file, "NS. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->mass*mass_unit, ptcl->stellar_evolution->get_zams(), ptcl->radius*position_unit, ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        if (ptcl->stellar_evolution->vkick[3] > 0.0) {
            fprintf(sevn_output_file, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->stellar_evolution->vkick[0], ptcl->stellar_evolution->vkick[1], ptcl->stellar_evolution->vkick[2]);
            for(int i=0; i<DIM; i++)
                ptcl->velocity[i] += ptcl->stellar_evolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->cm_particle_index != -1)
                ptcl->set_binary_interrupt_state(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->stellar_evolution->amiBH()) {
        // SEVNList.erase(ptcl->particle_index);
        // fprintf(sevn_output_file, "BH. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->stellar_evolution->getp(Mass::ID), ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        setBHspin(ptcl);
        ptcl->delta_mass += ptcl->mass - ptcl->stellar_evolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->mass = ptcl->stellar_evolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        // ptcl->world_time = NUMERIC_FLOAT_MAX;
        fprintf(sevn_output_file, "BH. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->mass*mass_unit, ptcl->stellar_evolution->get_zams(), ptcl->radius*position_unit, ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        fprintf(sevn_output_file, "\tDIMless spin. mag: %e, (%e, %e, %e)\n", ptcl->stellar_evolution->getp(Xspin::ID), ptcl->spin_param[0], ptcl->spin_param[1], ptcl->spin_param[2]);
        if (ptcl->stellar_evolution->vkick[3] > 0.0) {
            fprintf(sevn_output_file, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->stellar_evolution->vkick[0], ptcl->stellar_evolution->vkick[1], ptcl->stellar_evolution->vkick[2]);
            for(int i=0; i<DIM; i++)
                ptcl->velocity[i] += ptcl->stellar_evolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->cm_particle_index != -1)
                ptcl->set_binary_interrupt_state(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->stellar_evolution->amiempty()) {
        // SEVNList.erase(ptcl->particle_index);
        // fprintf(sevn_output_file, "Empty. PID: %d, Mass: %e Msol, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->stellar_evolution->getp(Mass::ID), ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        ptcl->delta_mass += ptcl->mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->mass = -1.0;
        fprintf(sevn_output_file, "Empty. PID: %d, ZAMS Mass: %e Msol, Time: %e Myr, Worldtime: %e Myr\n", ptcl->pid, ptcl->stellar_evolution->get_zams(), ptcl->world_time, ptcl->stellar_evolution->getp(Worldtime::ID));
        if (ptcl->cm_particle_index != -1) {
            ptcl->set_binary_interrupt_state(BinaryInterruptState::kicked);
        }
        else {
            if (!ptcl->is_active) {
                fprintf(stderr, "Particle %d is not active\n", ptcl->pid);
                fflush(stderr);
                assert(ptcl->is_active);
            }
            // assert(ptcl->is_active); // for debugging by EW 2025.1.20
            ptcl->is_active = false;
            num_particles--;
        }
        delete ptcl->stellar_evolution;
        ptcl->stellar_evolution = nullptr;
    }
#ifdef SEVN_BINARY
    if (from_BSE && ptcl->stellar_evolution->amiremnant()) {
        auto it = SEVNList.begin();
        while (it != SEVNList.end()) {
            if (it->second == ptcl->particle_index) {
                it = SEVNList.erase(it);
                fprintf(sevn_output_file, "Remnant particle (PID: %d) is deleted from SEVNList\n", ptcl->pid);
                break;
            }
            else
                it++;
        }
    }
#endif
}

// Reference: int Mix::special_evolve(Binstar *binstar) in Processes.cpp of SEVN
void Mix(StarSEVN* star1, StarSEVN* star2) {

    // utilities::wait("Hey I have to mix",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);

    StarSEVN *donor = star1; //Donor is the star that will set to empty
    StarSEVN *accretor = star2; //Accretor is the star that will remain as results of the mix

    ///Choose the star that remains and the one that is set to empty
    //First handle the general case: the accretor is the more evolved star,
    // or the more compact (massive) remnant if both are remnant (handled in get_id_more_evolved_star).

    //If one of the star is empty return the other one
    int get_id_more_evolved_star;

    //In the stars have the same phase (not remnant handled before), return the more evolved is the one with the larger plife
    if (donor->getp(Phase::ID)==accretor->getp(Phase::ID))
        get_id_more_evolved_star = donor->plife()>=accretor->plife() ? 0 : 1;
    else
        get_id_more_evolved_star = donor->getp(Phase::ID)>=accretor->getp(Phase::ID) ? 0 : 1;
    if(get_id_more_evolved_star==0)
        //Now the star 0 becomes the accretor and the star 1 the donor
        utilities::swap_stars(donor,accretor);

    //Now handle a special case: If the accretor is a naked helium and the donor not and they have the same innermost core
    //we swap so that the accretor will be the normal star.  This is done because jumping to a normal
    //track from a pureHe is more dangerous since if there are problems on finding a good match the star
    //does not jump and the code raises an error.  In case of a normal track, if we don't find a match
    //we can just continue following the original track.
    //GI 18/02/21: bug fix, to really swap the donor have to have at least a Helium core (the CO core is handled inside).
    //We check >1E-3 instead of >0 to avoid to select the accretor star as the one that is just starting to grow its HE core.
    //GI 02/09/22: We further extend the bug fix to disable the donor swab if the H-star is growing its core from 0, i.e. it is the TMS phase
    //
    //if (accretor->aminakedhelium() and (!donor->aminakedhelium() and donor->getp(MHE::ID)>1E-3)){
    if (accretor->aminakedhelium() and (!donor->aminakedhelium() and
        donor->getp(Phase::ID)>Lookup::TerminalMainSequence and
            donor->getp(Phase::ID)!=Lookup::TerminalCoreHeBurning)){
        unsigned int id_inner_core_accretor = accretor->getp(MCO::ID)!=0 ? MCO::ID : MHE::ID;
        unsigned int id_inner_core_donor =  donor->getp(MCO::ID)!=0 ? MCO::ID : MHE::ID;
        if (id_inner_core_accretor==id_inner_core_donor)
            utilities::swap_stars(donor,accretor);
    }



    ///Let's mix
    //Case 1 - Mix between not remnant stars
    if (!accretor->amiremnant() and !donor->amiremnant()){

        //Set new maximumCO and minmum HE
        accretor->set_MCO_max_aftermerge(donor);
        accretor->set_MHE_min_aftermerge(donor);

        //Mass from the donor
        double DM_new     = donor->getp(Mass::ID);
        double DMHE_new   = donor->getp(MHE::ID);
        double DMCO_new   = donor->getp(MCO::ID);
        double MCO_old    = accretor->getp(MCO::ID);
        double MHE_old    = accretor->getp(MHE::ID);

        //Update masses of the accretor
        accretor->update_from_binary(Mass::ID, DM_new);
        accretor->update_from_binary(dMcumul_binary::ID, DM_new);
        accretor->update_from_binary(MHE::ID, DMHE_new);
        accretor->update_from_binary(MCO::ID, DMCO_new);

        ///Jump to a new track
        //Sanity check on MCO
        if (MCO_old==0 and DMCO_new>0){
            throw std::runtime_error("This is an embarrassing situation, in Mix the donor has a CO core and the accretor not");
            // svlog.critical("This is an embarrassing situation, in Mix the donor has a CO core and the accretor not",
            //                 __FILE__,__LINE__,sevnstd::ce_error());
        }
        // Sanity check on MHE
        else if (MHE_old==0 and DMHE_new>0){
            throw std::runtime_error("This is an embarrassing situation, in Mix the donor has a HE core and the accretor not");
            // svlog.critical("This is an embarrassing situation, in Mix the donor has a HE core and the accretor not",
            //                 __FILE__,__LINE__,sevnstd::ce_error());
        }
        // The accretor is a nakedHelium and the donor not, so we have to trigger a jump to a normal track
        // There is no possibility that the donor has a CO core and the naked helium not because we already check this
        // at the beginning
        else if (accretor->aminakedhelium() and !donor->aminakedhelium()){
            accretor->jump_to_normal_tracks();
        }
        else{
            accretor->find_new_track_after_merger();
        }
    }
    //Case 2 - Mix between WDs
    //TODO Notice that here mixing with a WD is treated as mixing with a NS/BH, in Hurley we have a more complicate outcomes
    else if (accretor->amiWD() and donor->amiWD()) {

        //Check if we have two HeWD
        bool check_double_HeWD = donor->getp(RemnantType::ID)==Lookup::Remnants::HeWD && accretor->getp(RemnantType::ID)==Lookup::Remnants::HeWD;

        //Case 2a - SNI explosion
        if (check_double_HeWD){
            accretor->explode_as_SNI(); // Eunwoo: this should be treated carefully!
            // if (binstar->onesurvived) binstar->set_onesurvived(); //Why this condition? I forgot, we have to check
        }
        else{
            //Update the Mass
            accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
            // binstar->set_onesurvived();
        }

    }
    //Case 3 -  NS, BH or WD mixing with a NS/BH
    else if (accretor->amiCompact() and donor->amiremnant()){
        //Update the Mass
        accretor->update_from_binary(Mass::ID, donor->getp(Mass::ID));
    }
    //Case 4 - A not remnant donor over a NS, BH, WD
    //-> Do not accreate nothing, just set the donor to empty, see below
    //From commonenvelope::main_stellar_collision in common_envelope.cpp in SEVN1


    ///LOGPRINT AND EVENT SET
    // std::string w =Mix::log_message(binstar,accretor,donor);
    // binstar->print_to_log(w);
    //utilities::hardwait(" Event before ",get_event());
    // set_event((double)Lookup::EventsList::Merger);
    //utilities::hardwait(" Event After ",get_event());

    ///LAST STUFF
    //In any case set the binary to broken and put the donor to empty
    //donor->set_empty_in_bse(binstar);
    donor->set_empty();

}

// Use this function when merger happened
void SetRadius(Particle* ptcl) {

    if (ptcl->particle_type < REMNANT) {
        ptcl->radius = ptcl->stellar_evolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->mass*mass_unit > ptcl->stellar_evolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->mass*mass_unit/ptcl->stellar_evolution->get_max_zams(), 1./3);
    }
    else if (ptcl->particle_type <= WHITE_DWARF_ONE) {
        double RNS = 11/(velocity_unit/yr*pc/1e5); // 11 km/s in code unit
        double Mch = 1.41/mass_unit;
        double RWD = 0.0115*std::sqrt(pow(Mch/ptcl->mass,0.6666666667) -  pow(ptcl->mass/Mch,0.6666666667));
        
        ptcl->radius = std::max(RNS,RWD);
    }
    else if (ptcl->particle_type <= NEUTRON_STAR_CCSN)
        ptcl->radius = 11/(velocity_unit/yr*pc/1e5); // 11 km/s in code unit
    else
        ptcl->radius = 2*ptcl->mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzschild radius in code unit
}

#ifdef SEVN_BINARY
/* by EW 2025.6.27
    Important note: In current version, binary members are DELETED in the SEVNList.
    If this takes too much time, we should consider to use a different approach.
    How about not using SEVNList at all, and applying stellar evolution at every irregular timestep?
*/
// Return false if a member is empty or kicked
bool makeSEVNBinary(Particle* ptclCM) {

    assert(!ptclCM->is_active);

    // Activate binary stellar evolution iff there are two members in the binary
	int NumberOfMembers=0;
    Particle* ptcl1 = &particles[ptclCM->new_members[0]];
	for (int i = 0; i < ptclCM->new_num_members; ++i) {
        Particle* members = &particles[ptclCM->new_members[i]];
        if (members->current_time_irr > ptcl1->current_time_irr) {
        	ptcl1 = members;
    	}
		if (!members->is_cm_particle)
			NumberOfMembers++;
		else {
			NumberOfMembers += members->num_members;
            if (members->binary_evolution != nullptr) {
                fprintf(sevn_output_file, "SEVN BSE... BSE can't be applied to many-body case. Binary object (PID: %d) should be deleted!!!\n", members->pid);
                convertBinaryToSingle(members);
            }
        }
    }
    if (NumberOfMembers != 2)
        return true;

    double BinaryFormationTime = ptcl1->current_time_irr; // in code unit
    
    ptcl1 = &particles[ptclCM->new_members[0]];
    StarSEVN* star1 = ptcl1->stellar_evolution;

    Particle* ptcl2 = &particles[ptclCM->new_members[1]];
    StarSEVN* star2 = ptcl2->stellar_evolution;

    // Make BSE object iff both members are not remnants
    // If they are remnants, let's consider orbit-shrinking, TDE, GW_merge in SDAR
    if ((ptcl1->particle_type > REMNANT) || (ptcl2->particle_type > REMNANT))
        return true;

    // Make BSE object iff both members have StellarEvolution objects
    if (star1 == nullptr || star2 == nullptr)
        return true;

    double pos1[3], pos2[3], vel1[3], vel2[3]; // position and velocity vectors at BinaryFormationTime
    double dr[3], dv[3], m_tot; // relative position and velocity vectors, total mass
    double semi, ecc; // semi-major axis and eccentricity
    double rv; // inner product of relative position and velocity vectors

    ptcl1->predict_particle_second_order(BinaryFormationTime - ptcl1->current_time_irr, pos1, vel1);
    ptcl2->predict_particle_second_order(BinaryFormationTime - ptcl2->current_time_irr, pos2, vel2);

    for (int dim=0; dim<DIM; dim++) {
        dr[dim] = pos2[dim] - pos1[dim];
        dv[dim] = vel2[dim] - vel1[dim];
    }
    double r = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    rv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];

    m_tot = ptcl1->mass + ptcl2->mass;
    double v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
    semi = 1.0 / (2.0 / r - v2 / m_tot); // semi-major axis in code unit

    double p = 1.0 - r/semi;

    ecc = sqrt(p*p + rv*rv/semi/m_tot);

    // Apply binary stellar evolution iff the binary orbit is elliptical
    if (ecc >= 1.0)
        return true;

    // for debugging by EW 2025.6.26
    // fprintf(sevn_output_file, "SEVN BSE... CM PID: %d. semi: %e pc, ecc: %e\n", ptclCM->pid, semi * position_unit, ecc);

    semi = semi * position_unit * utilities::parsec_to_Rsun; // semi-major axis in Rsun

    // fprintf(sevn_output_file, "SEVN BSE... Binary candidate members should be deleted from SEVNList\n");
    int num_del = 0;
    auto it = SEVNList.begin();
    while (it != SEVNList.end()) {
        if (it->second == ptcl1->particle_index || it->second == ptcl2->particle_index) {
            Particle* ptcl = &particles[it->second];
            it = SEVNList.erase(it);
            num_del++;
            fprintf(sevn_output_file, "BSE candidate member (PID: %d) is deleted from SEVNList\n", ptcl->pid);
            if (num_del == 2)
                break;
        }
        else
            it++;
    }
    assert(num_del == 2); // we should delete two particles from the SEVNList

    BinaryFormationTime *= enzo_time_step*1e4; // in Myr unit
    /*
    bool evolve = false;
    while (ptcl1->world_time + ptcl1->stellar_evolution->getp(Timestep::ID) <= BinaryFormationTime) {
        ptcl1->world_time += ptcl1->stellar_evolution->getp(Timestep::ID);
        ptcl1->stellar_evolution->evolve();
        evolve = true;
    }
    if (evolve) {
        UpdateEvolution(ptcl1, true);
        if (ptcl1->particle_type > REMNANT) {
            if (star1->amiempty() || (star1->vkick[3] > 0.0)) {
                fprintf(sevn_output_file, "SEVN BSE... No binary is created!!!\n");
                ptclCM->new_num_members = 0;
                return false;
            }
            return true;
        }
    }
    ptcl1->world_time = BinaryFormationTime;

    evolve = false;
    while (ptcl2->world_time + ptcl2->stellar_evolution->getp(Timestep::ID) <= BinaryFormationTime) {
        ptcl2->world_time += ptcl2->stellar_evolution->getp(Timestep::ID);
        ptcl2->stellar_evolution->evolve();
        evolve = true;
    }
    if (evolve) {
        UpdateEvolution(ptcl2, true);
        if (ptcl2->particle_type > REMNANT) {
            if (star2->amiempty() || (star2->vkick[3] > 0.0)) {
                fprintf(sevn_output_file, "SEVN BSE... No binary is created!!!\n");
                ptclCM->new_num_members = 0;
                return false;
            }
            return true;
        }
    }
    ptcl2->world_time = BinaryFormationTime;
    */

    std::vector<std::string> init1 = getCustomInitParams(ptcl1);
    std::vector<std::string> init2 = getCustomInitParams(ptcl2);

    std::vector<std::string> init_binstar_params(init1.begin(), init1.begin()+5);
    init_binstar_params.insert(init_binstar_params.end(), init2.begin(), init2.begin() + 5);

    std::stringstream semi_str;
    semi_str << std::setprecision(17) << semi; // semi-major axis in Rsun

    std::stringstream ecc_str;
    ecc_str << std::setprecision(17) << ecc; // eccentricity

    init_binstar_params.insert(init_binstar_params.end(), {semi_str.str(), ecc_str.str(), "broken", "all"});

    size_t id = ptclCM->pid;
    Binstar* binstar = new Binstar(sevnio, init_binstar_params, id);
    fprintf(sevn_output_file, "SEVN BSE... New Binstar is successfully created. CM PID: %d at %e Myr (semi: %e pc, ecc: %e, GWtime: %e Myr)\n", 
            ptclCM->pid, BinaryFormationTime, binstar->getp(Semimajor::ID) / utilities::parsec_to_Rsun, binstar->getp(Eccentricity::ID), binstar->getp(GWtime::ID));

    ptclCM->binary_evolution = binstar;
    ptclCM->formation_time   = BinaryFormationTime;    // in Myr unit
    ptclCM->world_time       = BinaryFormationTime;    // in Myr unit

    return true;
}

void convertBinaryToSingle(Particle* ptclCM) {

    assert(ptclCM->num_members == 2);

    Particle* ptcl1 = &particles[ptclCM->members[0]];
    Particle* ptcl2 = &particles[ptclCM->members[1]];
    
    assert(ptcl1->stellar_evolution == nullptr);
    assert(ptcl2->stellar_evolution == nullptr);

    StarSEVN* star1 = ptclCM->binary_evolution->getstar(0);
    StarSEVN* star2 = ptclCM->binary_evolution->getstar(1);
    if ((int)star1->get_ID() != 0) {
        fprintf(sevn_output_file, "In convertBinaryToSingle.. star1->get_ID(): %d, star2->get_ID(): %d\n", (int)star1->get_ID(), (int)star2->get_ID());
        assert((int)star1->get_ID() == 1);
        std::swap(star1, star2);
    }

    ptcl1->stellar_evolution = star1;
    ptcl2->stellar_evolution = star2;

    ptclCM->binary_evolution->custom_destructor();
    ptclCM->binary_evolution = nullptr;
    fprintf(sevn_output_file, "In convertBinaryToSingle... Binary object (PID: %d) SEVN memory is free now\n", ptclCM->pid);

    if (!star1->amiremnant())
        SEVNList.insert({ptcl1->world_time + ptcl1->stellar_evolution->getp(Timestep::ID), ptcl1->particle_index});
    else if (star1->amiempty()) {
        delete star1;
        ptcl1->stellar_evolution = nullptr;
    }
    if (!star2->amiremnant())
        SEVNList.insert({ptcl2->world_time + ptcl2->stellar_evolution->getp(Timestep::ID), ptcl2->particle_index});
    else if (star2->amiempty()) {
        delete star2;
        ptcl2->stellar_evolution = nullptr;
    }

}

std::vector<std::string> getCustomInitParams(Particle* ptcl) {

    StarSEVN* star = ptcl->stellar_evolution;
    
    std::stringstream Mass;
    if (star->amiremnant()) {
        int remnant_type = int(star->getp(RemnantType::ID));
        std::string remnant_suffix;

        if (remnant_type == 1)
            remnant_suffix = "HEWD";
        else if (remnant_type == 2)
            remnant_suffix = "COWD";
        else if (remnant_type == 3)
            remnant_suffix = "ONEWD";
        else if (remnant_type == 4)
            remnant_suffix = "NSEC";
        else if (remnant_type == 5)
            remnant_suffix = "NS";
        else if (remnant_type == 6)
            remnant_suffix = "BH";

        Mass << std::setprecision(17) << star->getp(Mass::ID) << remnant_suffix;
    }
    else if (star->aminakedhelium()) { // This includes naked helium stars and naked CO stars by EW 2025.7.3 // MHE always includes MCO in SEVN
        Mass << std::setprecision(17) << "(" << star->get_zams() << "," 
            << star->getp(Mass::ID) << "," 
            << star->getp(MCO::ID) << ")HE";
    }
    else {
        Mass << std::setprecision(17) << "(" << star->get_zams() << "," 
            << star->getp(Mass::ID) << "," 
            << star->getp(MHE::ID) << "," 
            << star->getp(MCO::ID) << ")";
    }
    

    std::stringstream Z;
    Z << std::setprecision(17) << star->get_Z();

    std::stringstream tini;
    if (star->amiremnant())
        tini << "zams";
    else {
        tini << std::setprecision(17) << "%" << star->plife() * 100
            << ":" << static_cast<int>(star->getp(Phase::ID));
    }
    
    std::vector<std::string> init_params{Mass.str(), Z.str(), "0.0", "delayed", tini.str(), "end", "all"};

    delete star;
    ptcl->stellar_evolution = nullptr;

    return init_params;
}

void deleteSEVNBinary(Particle* ptclCM) {

    convertBinaryToSingle(ptclCM);
    return;
}

void BinaryEvolution(Particle* ptclCM) {

    Binstar* binary = ptclCM->binary_evolution;
    
    if (ptclCM->world_time + binary->getp(BTimestep::ID) > (ptclCM->current_time_irr + ptclCM->time_step_irr) * enzo_time_step * 1e4)
        return;

    if (ptclCM->get_binary_interrupt_state() == BinaryInterruptState::merger) {
        convertBinaryToSingle(ptclCM);
        return;
    }

    if (ptclCM->spin_param[1] >= 1.0) {
        fprintf(sevn_output_file, "In BinaryEvolution... Initial binary orbit became hyperbolic in SDAR!!!\n");
        convertBinaryToSingle(ptclCM);
        return;
    }

    while ((ptclCM->world_time + binary->getp(BTimestep::ID) <= (ptclCM->current_time_irr + ptclCM->time_step_irr) * enzo_time_step * 1e4) &&
            !binary->getstar(0)->amiremnant() &&
            !binary->getstar(1)->amiremnant()) {

        ptclCM->world_time += binary->getp(BTimestep::ID);
        for (int i=0; i<ptclCM->num_members; i++) {
            Particle* member = &particles[ptclCM->members[i]];
            member->world_time += binary->getp(BTimestep::ID);
        }
        binary->evolve();
    }
    UpdateBinaryEvolution(ptclCM);
}

void UpdateBinaryEvolution(Particle* ptclCM) {
    
    /*
    There are 20 binary events...
    Let's check what happens to Binstar if PISN happen for a member -> broken or not?
    How about merger/CE driven merger/RLOF driven merger -> broken or not / what happens to donor? 

    If merger-like things happen, we should delete Binstar, and terminate SDAR too! -> FBTermination?
    After binary evolution, let's correct semi-major axis and eccentricity by calculating orbital parameters

    Let's store newly calculated semi-major axis and eccentricity in ptclCM->spin_param (+ CM ptcl mass?)
    */

    Binstar* binary = ptclCM->binary_evolution;

    StarSEVN* star1 = binary->getstar(0);
    Particle* ptcl1 = &particles[ptclCM->members[0]];

    StarSEVN* star2 = binary->getstar(1);
    Particle* ptcl2 = &particles[ptclCM->members[1]];

    if ((int)star1->get_ID() != 0){
        fprintf(sevn_output_file, "In UpdateBinaryEvolution... star1->get_ID(): %d, star2->get_ID(): %d\n", (int)star1->get_ID(), (int)star2->get_ID());
        assert((int)star1->get_ID() == 1);
        std::swap(star1, star2);
    }

    double dm = 0.0; // ejected mass by stellar wind to nearby gas cells

    fprintf(sevn_output_file, "BSE... CM PID: %d, BEvent: %d, WorldTime: %e Myr, GW merger time: %e Myr, semi: %e pc, ecc: %e\n",
            ptclCM->pid, binary->getp(BEvent::ID), ptclCM->world_time, 
            binary->getp(GWtime::ID), binary->getp(Semimajor::ID)/utilities::parsec_to_Rsun, binary->getp(Eccentricity::ID));
    fprintf(sevn_output_file, "\tPhase1: %d, Mass1: %e Msol, Radius1: %e pc, T_eff1: %e K\n", 
            star1->getp(Phase::ID), star1->getp(Mass::ID), star1->getp(Radius::ID)/(utilities::parsec_to_Rsun), star1->getp(Temperature::ID));
    fprintf(sevn_output_file, "\tPhase2: %d, Mass2: %e Msol, Radius2: %e pc, T_eff2: %e K\n", 
            star2->getp(Phase::ID), star2->getp(Mass::ID), star2->getp(Radius::ID)/(utilities::parsec_to_Rsun), star2->getp(Temperature::ID));

    // Let's update ParticleType
    if (star1->amiremnant())
        ptcl1->particle_type = REMNANT + (int)star1->getp(RemnantType::ID);
    else
        ptcl1->particle_type = (int)star1->getp(Phase::ID);

    if (star2->amiremnant())
        ptcl2->particle_type = REMNANT + (int)star2->getp(RemnantType::ID);
    else
        ptcl2->particle_type = (int)star2->getp(Phase::ID);

    // Let's store the semi-major axis and eccentricity in ptclCM->spin_param if RLOF is triggered
    // a_spin[2] > 0.0: RLOF is triggered, and applied to SDAR
    // a_spin[2] = 0.0: RLOF is not triggered yet
    // a_spin[2] < 0.0: RLOF is triggered, but not applied to SDAR yet
    if (binary->getp(Eccentricity::ID) == 0.0 && ptclCM->spin_param[2] == 0.0) {
        ptclCM->spin_param[0] = binary->getp(Semimajor::ID)/utilities::parsec_to_Rsun/position_unit; // in code unit
        ptclCM->spin_param[1] = binary->getp(Eccentricity::ID);
        ptclCM->spin_param[2] = -1.0;
    }

    if (star1->amiempty()) {
        fprintf(sevn_output_file, "\tPID: %d. Empty!\n", ptcl1->pid);
        ptcl1->delta_mass += ptcl1->mass;
        ptcl1->mass = -1.0;
        ptcl1->set_binary_interrupt_state(BinaryInterruptState::kicked);
        ptclCM->set_binary_interrupt_state(BinaryInterruptState::terminated);
    } else {

        ptcl1->mass = star1->getp(Mass::ID)/mass_unit;
        ptcl1->radius = star1->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;

        if (!star2->amiempty()) {
            dm = ptclCM->mass - (ptcl1->mass + star2->getp(Mass::ID)/mass_unit);
            ptclCM->mass = ptcl1->mass + star2->getp(Mass::ID)/mass_unit;
        }

        if (star1->vkick[3] > 0.0) {
            fprintf(sevn_output_file, "\tPID: %d. Kicked velocity1: (%e, %e, %e) [km/s]\n", ptcl1->pid, star1->vkick[0], star1->vkick[1], star1->vkick[2]);
            ptcl1->set_binary_interrupt_state(BinaryInterruptState::kicked);
            ptclCM->set_binary_interrupt_state(BinaryInterruptState::terminated);
            for (int dim = 0; dim < DIM; dim++)
                ptcl1->velocity[dim] += star1->vkick[dim]/(velocity_unit/yr*pc/1e5);

            ptcl2->delta_mass += dm;
        } else
            ptcl1->delta_mass += dm;
    }

    if (star2->amiempty()) {
        fprintf(sevn_output_file, "\tPID: %d. Empty!\n", ptcl2->pid);
        ptcl2->delta_mass += ptcl2->mass;
        ptcl2->mass = -1.0;
        ptcl2->set_binary_interrupt_state(BinaryInterruptState::kicked);
        ptclCM->set_binary_interrupt_state(BinaryInterruptState::terminated);
    } else {

        ptcl2->mass = star2->getp(Mass::ID)/mass_unit;
        ptcl2->radius = star2->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;

        if (star2->vkick[3] > 0.0) {
            fprintf(sevn_output_file, "\tPID: %d. Kicked velocity2: (%e, %e, %e) [km/s]\n", ptcl2->pid, star2->vkick[0], star2->vkick[1], star2->vkick[2]);
            ptcl2->set_binary_interrupt_state(BinaryInterruptState::kicked);
            ptclCM->set_binary_interrupt_state(BinaryInterruptState::terminated);
            for (int dim = 0; dim < DIM; dim++)
                ptcl2->velocity[dim] += star2->vkick[dim]/(velocity_unit/yr*pc/1e5);
        }
    }

    if (star1->amiremnant() || star2->amiremnant())
        convertBinaryToSingle(ptclCM);

}
#endif

#endif