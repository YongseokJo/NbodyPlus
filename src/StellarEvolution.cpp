#ifdef SEVN

#include "global.h"
#include <random>
#include <map>

#ifdef SEVN_BINARY
void UpdateEvolution(Particle* ptcl, bool from_BSE);
void UpdateBinaryEvolution(Particle* ptcl);
#else
void UpdateEvolution(Particle* ptcl);
#endif

// Start SEVN stellar evolution
// Set stellar radii, BH spin, etc
void initializeStellarEvolution() {

    std::ostringstream oss;
	oss << sevnio->svpar;
	fprintf(SEVNout, "%s\n\n\n\n\n", oss.str().c_str());
    // fprintf(SEVNout, "PID\tPhase\tMass (Msun)\tRadius (pc)\tTime (Myr)\tWorldtime (Myr)\n")
	fflush(SEVNout);

    assert(SEVNList.empty());

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
        ptcl->StellarEvolution = new StarSEVN(sevnio, init_params, id, false);
        SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});

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

    Particle* ptcl;
    while (!SEVNList.empty()) {
         
        auto it = SEVNList.begin();
        ptcl = &particles[it->second];
#ifdef SEVN_BINARY
        if (ptcl->isCMptcl) {
            ptcl->WorldTime += ptcl->BinaryEvolution->getp(BTimestep::ID);
            for (int i=0; i<ptcl->NumberOfMember; i++) {
                Particle* member = &particles[ptcl->Members[i]];
                member->WorldTime += ptcl->BinaryEvolution->getp(BTimestep::ID);
            }
            ptcl->BinaryEvolution->evolve();

            while (ptcl->WorldTime + ptcl->BinaryEvolution->getp(BTimestep::ID) <= global_time * EnzoTimeStep * 1e4) {
                ptcl->WorldTime += ptcl->BinaryEvolution->getp(BTimestep::ID);
                for (int i=0; i<ptcl->NumberOfMember; i++) {
                    Particle* member = &particles[ptcl->Members[i]];
                    member->WorldTime += ptcl->BinaryEvolution->getp(BTimestep::ID);
                }
                ptcl->BinaryEvolution->evolve();
            }

            it = SEVNList.erase(it);
            UpdateBinaryEvolution(ptcl);

            if (SEVNList.empty() || SEVNList.begin()->first > global_time * EnzoTimeStep * 1e4)
                break;
        }
#endif
        ptcl->WorldTime += ptcl->StellarEvolution->getp(Timestep::ID);
        ptcl->StellarEvolution->evolve();

        while (ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID) <= global_time * EnzoTimeStep * 1e4) {
            ptcl->WorldTime += ptcl->StellarEvolution->getp(Timestep::ID);
            ptcl->StellarEvolution->evolve();
        }

        it = SEVNList.erase(it);
#ifdef SEVN_BINARY
        UpdateEvolution(ptcl, false);
#else
        UpdateEvolution(ptcl);
#endif

        if (SEVNList.empty() || SEVNList.begin()->first > global_time * EnzoTimeStep * 1e4)
            break;
    }
    fflush(SEVNout);
}

#ifdef SEVN_BINARY
void UpdateEvolution(Particle* ptcl, bool from_BSE) {
#else
void UpdateEvolution(Particle* ptcl) {
#endif

    if (!ptcl->StellarEvolution->amiremnant()) {
#ifdef SEVN_BINARY
        if (!from_BSE)
            SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});
#else
        SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});
#endif
        // fprintf(SEVNout, "PID: %d, Phase: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, int(ptcl->StellarEvolution->getp(Phase::ID)), ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->Mass*mass_unit > ptcl->StellarEvolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->Mass*mass_unit/ptcl->StellarEvolution->get_max_zams(), 1./3);
        fprintf(SEVNout, "PID: %d, Phase: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, int(ptcl->StellarEvolution->getp(Phase::ID)), ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->radius*position_unit, ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
    }
    else if (ptcl->StellarEvolution->amiWD()) {
        // SEVNList.erase(ptcl->ParticleIndex);
        // fprintf(SEVNout, "WD. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        // ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        fprintf(SEVNout, "WD. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->radius*position_unit, ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->StellarEvolution->amiNS()) {
        // SEVNList.erase(ptcl->ParticleIndex);
        // fprintf(SEVNout, "NS. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        // ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        fprintf(SEVNout, "NS. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->radius*position_unit, ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        if (ptcl->StellarEvolution->vkick[3] > 0.0) {
            fprintf(SEVNout, "\tKicked velocity: (%e, %e, %e) [km/s]\n", ptcl->StellarEvolution->vkick[0], ptcl->StellarEvolution->vkick[1], ptcl->StellarEvolution->vkick[2]);
            for(int i=0; i<Dim; i++)
                ptcl->Velocity[i] += ptcl->StellarEvolution->vkick[i]/(velocity_unit/yr*pc/1e5);
            if (ptcl->CMPtclIndex != -1)
                ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
    }
    else if (ptcl->StellarEvolution->amiBH()) {
        // SEVNList.erase(ptcl->ParticleIndex);
        // fprintf(SEVNout, "BH. PID: %d, Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        setBHspin(ptcl);
        ptcl->dm += ptcl->Mass - ptcl->StellarEvolution->getp(Mass::ID)/mass_unit; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = ptcl->StellarEvolution->getp(Mass::ID)/mass_unit;
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit; // this might be wrong!
        // ptcl->WorldTime = NUMERIC_FLOAT_MAX;
        fprintf(SEVNout, "BH. PID: %d, Mass: %e Msol, ZAMS Mass: %e Msol, Radius: %e pc, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->Mass*mass_unit, ptcl->StellarEvolution->get_zams(), ptcl->radius*position_unit, ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        fprintf(SEVNout, "\tDimless spin. mag: %e, (%e, %e, %e)\n", ptcl->StellarEvolution->getp(Xspin::ID), ptcl->a_spin[0], ptcl->a_spin[1], ptcl->a_spin[2]);
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
        // SEVNList.erase(ptcl->ParticleIndex);
        // fprintf(SEVNout, "Empty. PID: %d, Mass: %e Msol, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->getp(Mass::ID), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        ptcl->dm += ptcl->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
        ptcl->Mass = -1.0;
        fprintf(SEVNout, "Empty. PID: %d, ZAMS Mass: %e Msol, Time: %e Myr, Worldtime: %e Myr\n", ptcl->PID, ptcl->StellarEvolution->get_zams(), ptcl->WorldTime, ptcl->StellarEvolution->getp(Worldtime::ID));
        if (ptcl->CMPtclIndex != -1) {
            ptcl->setBinaryInterruptState(BinaryInterruptState::kicked);
        }
        else {
            if (!ptcl->isActive) {
                fprintf(stderr, "Particle %d is not active\n", ptcl->PID);
                fflush(stderr);
                assert(ptcl->isActive);
            }
            // assert(ptcl->isActive); // for debugging by EW 2025.1.20
            ptcl->isActive = false;
            NumberOfParticle--;
        }
        delete ptcl->StellarEvolution;
        ptcl->StellarEvolution = nullptr;
    }
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

    if (!ptcl->StellarEvolution->amiremnant()) {
        ptcl->radius = ptcl->StellarEvolution->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
        if (ptcl->Mass*mass_unit > ptcl->StellarEvolution->get_max_zams()) // VMS correction; constant stellar density is assumed
            ptcl->radius *= pow(ptcl->Mass*mass_unit/ptcl->StellarEvolution->get_max_zams(), 1./3);
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

#ifdef SEVN_BINARY
/* by EW 2025.6.27
    Important note: In current version, binary members are DELETED in the SEVNList.
    + CM ptcl will be added in the SEVNList.
    If this takes too much time, we should consider to use a different approach.
    How about not using SEVNList at all, and applying stellar evolution at every irregular timestep?
*/
// Return false if (P)PISN happens
bool makeSEVNBinary(Particle* ptclCM) {

    assert(!ptclCM->isActive);

    // Activate binary stellar evolution iff there are two members in the binary
	int NumberOfMembers=0;
    Particle* ptcl1 = &particles[ptclCM->NewNeighbors[0]];
	for (int i = 0; i < ptclCM->NewNumberOfNeighbor; ++i) {
		Particle* members = &particles[ptclCM->NewNeighbors[i]];
        if (members->CurrentTimeIrr > ptcl1->CurrentTimeIrr) {
        	ptcl1 = members;
    	}
		if (!members->isCMptcl)
			NumberOfMembers++;
		else {
			NumberOfMembers += members->NewNumberOfNeighbor;
            if (members->BinaryEvolution != nullptr) {

                assert(members->NewNumberOfNeighbor == 2);

                fprintf(SEVNout, "SEVN BSE... BSE can't be applied to many-body case. Binary object (PID: %d) should be deleted!!!\n", members->PID);

                auto it = SEVNList.begin();
                while (it != SEVNList.end()) {
                    if (it->second == members->ParticleIndex) {
                        it = SEVNList.erase(it);
                        fprintf(SEVNout, "Binary object (PID: %d) is deleted from SEVNList\n", members->PID);
                        break;
                    }
                    else
                        it++;
                }
                members->BinaryEvolution->custom_destructor();
                fprintf(SEVNout, "Binary object (PID: %d) SEVN memory is free now\n", members->PID);

                for (int j=0; j<members->NewNumberOfNeighbor; j++) {
                    Particle* ptcl = &particles[members->NewNeighbors[j]];
                    if (!ptcl->StellarEvolution->amiremnant())
                        SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});
                }
                return true;
            }
        }
    }
    if (NumberOfMembers != 2)
        return true;

    double BinaryFormationTime = ptcl1->CurrentTimeIrr; // in code unit
    
    ptcl1 = &particles[ptclCM->NewNeighbors[0]];
    StarSEVN* star1 = ptcl1->StellarEvolution;

    Particle* ptcl2 = &particles[ptclCM->NewNeighbors[1]];
    StarSEVN* star2 = ptcl2->StellarEvolution;

    // Make BSE object iff both members have StellarEvolution objects
    if (star1 == nullptr || star2 == nullptr)
        return true;

    double pos1[3], pos2[3], vel1[3], vel2[3]; // position and velocity vectors at BinaryFormationTime
    double dr[3], dv[3], M; // relative position and velocity vectors, total mass
    double energy; // specific orbital energy
    double semi, ecc; // semi-major axis and eccentricity
    double h[3]; // specific angular momentum vector

    ptcl1->predictParticleSecondOrder(BinaryFormationTime, pos1, vel1);
    ptcl2->predictParticleSecondOrder(BinaryFormationTime, pos2, vel2);

    for (int dim=0; dim<Dim; dim++) {
        dr[dim] = pos2[dim] - pos1[dim];
        dv[dim] = vel2[dim] - vel1[dim];
    }
    M = ptcl1->Mass + ptcl2->Mass;
    energy = 0.5 * mag(dv) - M / sqrt(mag(dr));
    semi = - M / (2.0 * energy); // semi-major axis in code unit

    h[0] = dr[1] * dv[2] - dr[2] * dv[1];
    h[1] = dr[2] * dv[0] - dr[0] * dv[2];
    h[2] = dr[0] * dv[1] - dr[1] * dv[0];

    ecc = sqrt(1 + (2 * energy * mag(h)) / (M * M));

    // Apply binary stellar evolution iff the binary orbit is bound
    if (ecc < 0.0)
        return true;

    // for debugging by EW 2025.6.26
    fprintf(SEVNout, "SEVN BSE... CM PID: %d. semi: %e pc, ecc: %e\n", ptclCM->PID, semi * position_unit, ecc);

    semi = semi * position_unit * utilities::parsec_to_Rsun; // semi-major axis in Rsun

    BinaryFormationTime *= EnzoTimeStep*1e4; // in Myr unit
    double dt; // timestep in Myr

    if (BinaryFormationTime > ptcl1->WorldTime) {

        dt = BinaryFormationTime - ptcl1->WorldTime;
        star1->sync_with(dt);
        star1->evolve();
        ptcl1->WorldTime = BinaryFormationTime;
        UpdateEvolution(ptcl1, true);

        if (star1->amiremnant()) {
            auto it = SEVNList.begin();
            while (it != SEVNList.end()) {
                if (it->second == ptcl1->ParticleIndex) {
                    it = SEVNList.erase(it);
                    fprintf(SEVNout, "Remnant particle (PID: %d) is deleted from SEVNList\n", ptcl1->PID);
                    break;
                }
                else
                    it++;
            }
            if (star1->amiempty() || (star1->vkick[3] > 0.0)) {
                fprintf(SEVNout, "SEVN BSE... No binary is created!!!\n");
                ptclCM->NewNumberOfNeighbor = 0;
                return false;
            }
        }
    }
    if (BinaryFormationTime > ptcl2->WorldTime) {

        dt = BinaryFormationTime - ptcl2->WorldTime;
        star2->sync_with(dt);
        star2->evolve();
        ptcl2->WorldTime = BinaryFormationTime;
        UpdateEvolution(ptcl2, true);

        if (star2->amiremnant()) {
            auto it = SEVNList.begin();
            while (it != SEVNList.end()) {
                if (it->second == ptcl2->ParticleIndex) {
                    it = SEVNList.erase(it);
                    fprintf(SEVNout, "Remnant particle (PID: %d) is deleted from SEVNList\n", ptcl2->PID);
                    break;
                }
                else
                    it++;
            }
            if (star2->amiempty() || (star2->vkick[3] > 0.0)) {
                fprintf(SEVNout, "SEVN BSE... No binary is created!!!\n");
                ptclCM->NewNumberOfNeighbor = 0;
                return false;
            }
        }
    }

    std::stringstream star1_percent;
    star1_percent << "%" << static_cast<int>(star1->plife() * 100) 
        << ":" << static_cast<int>(star1->getp(Phase::ID));

    std::stringstream star2_percent;
    star2_percent << "%" << static_cast<int>(star2->plife() * 100)
        << ":" << static_cast<int>(star2->getp(Phase::ID));

    // Mass, metallicity, spin, sn model, tini for star1
    // Mass, metallicity, spin, sn model, tini for star2
    // Semi-major axis (Rsun), eccentricity, tf, dtout, random seed(optional)
    std::vector<std::string> init_params_bin{
        std::to_string(star1->get_zams()), std::to_string(star1->get_Z()), "0.0", "delayed", star1_percent.str(), 
        std::to_string(star2->get_zams()), std::to_string(star2->get_Z()), "0.0", "delayed", star2_percent.str(), 
        std::to_string(semi), std::to_string(ecc), "broken", "all"
    };

    size_t id = ptclCM->PID;
    Binstar* binstar = new Binstar(sevnio, init_params_bin, star1, star2, id);
    fprintf(SEVNout, "SEVN BSE... New Binstar is successfully created. CM PID: %d\n");

    ptclCM->BinaryEvolution = binstar;
    ptclCM->FormationTime = BinaryFormationTime; // in Myr unit
    ptclCM->WorldTime = BinaryFormationTime; // in Myr unit

    int total = 0;
    if (!star1->amiremnant())
        total++;
    if (!star2->amiremnant())
        total++;
    if (total > 0) {
        fprintf(SEVNout, "SEVN BSE... %d objects should be deleted from SEVNList\n", total);

        int num_del = 0;
        auto it = SEVNList.begin();
        while (it != SEVNList.end()) {
            if (it->second == ptcl1->ParticleIndex || it->second == ptcl2->ParticleIndex) {
                Particle* ptcl = &particles[it->second];
                it = SEVNList.erase(it);
                num_del++;
                fprintf(SEVNout, "PISN induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl->PID);
                if (num_del == total)
                    break;
            }
            else
                it++;
        }
        assert(num_del == total); // we should delete two particles from the SEVNList
    }

    SEVNList.insert({ptclCM->WorldTime + binstar->getp(Timestep::ID), ptclCM->ParticleIndex}); // (BSE Query) This should be treated very carefully!!! by EW 2025.6.27
    fprintf(SEVNout, "SEVN BSE... BSE object (PID: %d) is inserted into SEVNList\n", ptclCM->PID);

    return true;
}

void deleteSEVNBinary(Particle* ptclCM) {

    if (ptclCM->BinaryEvolution == nullptr)
        return;

    assert(ptclCM->NumberOfMember == 2);

    fprintf(SEVNout, "SEVN BSE... Binary object (PID: %d) should be deleted\n", ptclCM->PID);

    auto it = SEVNList.begin();
    while (it != SEVNList.end()) {
        if (it->second == ptclCM->ParticleIndex) {
            it = SEVNList.erase(it);
            fprintf(SEVNout, "Binary object (PID: %d) is deleted from SEVNList\n", ptclCM->PID);
            break;
        }
        else
            it++;
    }
    ptclCM->BinaryEvolution->custom_destructor();
    fprintf(SEVNout, "Binary object (PID: %d) SEVN memory is free now\n", ptclCM->PID);

    for (int j=0; j<ptclCM->NumberOfMember; j++) {
        Particle* ptcl = &particles[ptclCM->Members[j]];
        assert(ptcl->StellarEvolution != nullptr);
        if (!ptcl->StellarEvolution->amiremnant())
            SEVNList.insert({ptcl->WorldTime + ptcl->StellarEvolution->getp(Timestep::ID), ptcl->ParticleIndex});
    }
    return;
}

void UpdateBinaryEvolution(Particle* ptclCM) {
    
    /*
    There are 20 binary events...
    Let's check what happens to Binstar if PISN happen for a member -> broken or not?
    How about merger/CE driven merger/RLOF driven merger -> broken or not / what happens to donor? 

    If merger-like things happen, we should delete Binstar, and terminate SDAR too! -> FBTermination?
    After binary evolution, let's correct semi-major axis and eccentricity by calculating orbital parameters
    */

    
}
#endif

#endif