#ifdef SEVN

#include "global.h"
#include <random>

void UpdateEvolution(Particle* ptcl);

// Start SEVN stellar evolution
// Set stellar radii, BH spin, etc
void initializeStellarEvolution() {
/*
    std::vector<std::string> args = {"empty", // Not used
                                // "-myself", "/data/vinicius/NbodyPlus/SEVN",
                                "-tables", "/data/vinicius/sevn/tables/SEVNtracks_parsec_ov04_AGB", 
                                //  "-tables", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_MIST_AGB",
                                // "-tables_HE", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_parsec_pureHe36",
                                // "-turn_WR_to_pureHe", "false",
                                "-snmode", "delayed",
                                "-Z", "0.0002",
                                "-spin", "0.0",
                                "-tini", "zams", 
                                "-tf", "end",
                                // "-tf", "0.000122",
                                "-dtout", "events",
                                "-xspinmode", "geneva"};
    std::vector<char*> c_args;
    for (auto& arg : args) {
        c_args.push_back(&arg[0]);
    }

    IO* sevnio;
    sevnio = new IO;
    sevnio->load(c_args.size(), c_args.data());
*/

    std::ostringstream oss;
	oss << sevnio->svpar;
	fprintf(SEVNout, "%s\n\n\n\n\n", oss.str().c_str());
    // fprintf(SEVNout, "PID\tPhase\tMass (Msun)\tRadius (pc)\tTime (Myr)\tWorldtime (Myr)\n")
	fflush(SEVNout);

    for (int i=0; i<NumberOfParticle; i++) {

        Particle* ptcl = &particles[i];

        ptcl->ParticleType = NormalStar+SingleStar;
        ptcl->FormationTime = 0.0;
        ptcl->WorldTime = 0.0;

        if (ptcl->Mass*mass_unit < 2.2) {
			ptcl->radius = 2.25461e-8/position_unit*pow(ptcl->Mass*mass_unit, 1/3);
            // stellar radius in code unit
            // extrapolated from the solar radius (it is assumed that low-mass stars have the same stellar density to the sun)

			// fprintf(SEVNout, "PID: %d. Mass: %e Msol, Radius: %e pc\n", ptcl->PID, ptcl->Mass*mass_unit, ptcl->radius*position_unit);
			continue;
		}

        // Mass, metallicity, spin, sn model, tini, tf, dtout, random seed(optional)
		std::vector<std::string> init_params{std::to_string(double(ptcl->Mass*mass_unit)), "0.0002", "0.0", "delayed", "zams", "end", "events"};

        size_t id = ptcl->PID;
        fprintf(stdout, "PID: %d, Mass: %e Msol\n", ptcl->PID, ptcl->Mass*mass_unit);
        fflush(stdout);
        ptcl->StellarEvolution = new Star(sevnio, init_params, id, false);

		ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // stellar radius in code unit
    }
}

void setBHspin(Particle* ptcl) {
    std::random_device rd; // Obtain a random number from hardware
    std::mt19937 mt(rd()); // Seed the generator
    std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
    double phi = 2 * M_PI * distr(mt);
    double theta = M_PI * distr(mt);

    ptcl->a_spin[0] = ptcl->StellarEvolution->getp(Xspin::ID) * sin(theta)*cos(phi);
    ptcl->a_spin[1] = ptcl->StellarEvolution->getp(Xspin::ID) * sin(theta)*sin(phi);
    ptcl->a_spin[2] = ptcl->StellarEvolution->getp(Xspin::ID) * cos(theta);
}

void StellarEvolution() {

    bool evolved;
    int OriginalNumberOfParticle = NumberOfParticle;
    Particle* ptcl;

    for (int i=0; i<=LastParticleIndex; i++) {

        ptcl = &particles[i];

        evolved = false;
    
        if (ptcl->StellarEvolution == nullptr || ptcl->StellarEvolution->amiremnant()) // CM particle & Stars with mass < 2.2 Msol don't have star class.
            continue;

        while (ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID) < global_time*EnzoTimeStep*1e4) {
            ptcl->WorldTime += ptcl->StellarEvolution->getp(Timestep::ID);
            ptcl->StellarEvolution->evolve();
            evolved = true;
        }

        if (evolved) {
            UpdateEvolution(ptcl);
        }          
    }

    if (OriginalNumberOfParticle != NumberOfParticle) {
        for (int i=0; i<=LastParticleIndex; i++) {

            ptcl = &particles[i];
            if (!ptcl->isActive)
                continue;

            auto newEnd = std::remove_if(
                ptcl->Neighbors, 
                ptcl->Neighbors + ptcl->NumberOfNeighbor, 
                [&particles](int j) {
                    return !particles[j].isActive;
                }
            );

            if (newEnd != ptcl->Neighbors + ptcl->NumberOfNeighbor)
                ptcl->NumberOfNeighbor = newEnd - ptcl->Neighbors;

        }
    }
}

void UpdateEvolution(Particle* ptcl) {

    if (!ptcl->StellarEvolution->amiremnant()) {
        fprintf(SEVNout, "PID: %d, Phase: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, int(ptcl->StellarEvolution->getp(Phase::ID)), ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->Mass*mass_unit > ptcl->StellarEvolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->Mass*1e9/ptcl->StellarEvolution->get_max_zams(), 1./3);
    }
    else if (ptcl->StellarEvolution->amiWD()) {
        fprintf(SEVNout, "WD. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->StellarEvolution->amiNS()) {
        fprintf(SEVNout, "NS. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->StellarEvolution->amiBH()) {
        fprintf(SEVNout, "BH. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        setBHspin(ptcl);
        fprintf(SEVNout, "\tDimless spin. mag: %e, (%e, %e, %e)\n", ptcl->StellarEvolution->getp(Xspin::ID), ptcl->a_spin[0], ptcl->a_spin[1], ptcl->a_spin[2]);
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        ptcl->ParticleType = Blackhole+SingleStar;
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->StellarEvolution->amiempty()) {
        fprintf(SEVNout, "Empty. PID: %d, Mass: %e Msol, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = 0.0;
        if (ptcl->CMPtclIndex != -1) {
            ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
        else {
            assert(ptcl->isActive); // for debugging by EW 2025.1.20
            ptcl->isActive = false;
            NumberOfParticle--;
        }
    }
    fflush(SEVNout);
}

// Reference: int Mix::special_evolve(Binstar *binstar) in Processes.cpp of SEVN
void Mix(Star* star1, Star* star2) {

    // utilities::wait("Hey I have to mix",binstar->getp(BWorldtime::ID),__FILE__,__LINE__);

    Star *donor = star1; //Donor is the star that will set to empty
    Star *accretor = star2; //Accretor is the star that will remain as results of the mix

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

    if (!ptcl->StellarEvolution->amiremnant()) {
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->Mass*mass_unit > ptcl->StellarEvolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->Mass*1e9/ptcl->StellarEvolution->get_max_zams(), 1./3);
    }
    else if (ptcl->StellarEvolution->amiWD()) {
        double RNS = 11/(velocity_unit/yr*pc/1e5); // 11 km/s in code unit
        double Mch = 1.41/mass_unit;
        double RWD = 0.0115*std::sqrt(pow(Mch/ptcl->Mass,0.6666666667) -  pow(ptcl->Mass/Mch,0.6666666667));
        
        ptcl->radius = std::max(RNS,RWD);
    }
    else if (ptcl->StellarEvolution->amiNS())
        ptcl->radius = 11/(velocity_unit/yr*pc/1e5); // 11 km/s in code unit
    else if (ptcl->StellarEvolution->amiBH())
        ptcl->radius = 2*ptcl->Mass/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzschild radius in code unit
}

#endif