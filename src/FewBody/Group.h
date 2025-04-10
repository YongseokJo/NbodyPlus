#ifndef GROUP_H
#define GROUP_H
#ifdef FEWBODY

#include <vector>
#include <iostream>
#include "../def.h"
#include <cassert>

#include "Common/binary_tree.h"
#include "Common/Float.h"
#include "Common/io.h"
#include "AR/symplectic_integrator.h"
#include "AR/information.h"
#include "ar_interaction.hpp"
#include "ar_perturber.hpp"

struct Particle;
struct Group
{

	Particle* groupCM;
	bool isTerminate; // For later use: I will use this when Binary Interrupt State is being used
	bool isMerger; // True if merger happened

	double CurrentTime;  // this show how much the binary system has evolved

	AR::TimeTransformedSymplecticIntegrator<Particle, Particle, Perturber, Interaction, AR::Information<Particle,Particle>> sym_int;
	AR::TimeTransformedSymplecticManager<Interaction> manager;


	// Constructor
	Group(void) 
		: groupCM(nullptr),
		isTerminate(false),
		isMerger(false),
		CurrentTime(0.0),
		sym_int(),
		manager()
	{}

	~Group() {
		sym_int.clear();
		sym_int.particles.clear(); // This might be unnecessary by EW 2025.1.4
		manager.step.clear();
		groupCM = nullptr;
	}

	void ARIntegration(double next_time);
	bool CheckBreak();
	void initialManager();
	void initialIntegrator(int NumMembers);

};
#endif
#endif