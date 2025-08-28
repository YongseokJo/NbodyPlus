#include <iostream>
#include <stdio.h>
#include "global.h"
#ifdef SEVN
#include <map>
#endif

Particle *particles_original;
Particle *particles;

MPI_Win win;
MPI_Win win2;

MPI_Comm shared_comm;
int MyRank;
int NumberOfProcessor;
int NumberOfWorker;

MPI_Datatype QueueType;
MPI_Datatype IparticleType;
MPI_Datatype JparticleType;

GlobalVariable *global_variable;
GlobalVariable *global_variable_original;
int LastParticleIndex; // The last index of particle array
int NumberOfParticle; // The number of active particles
int NewPID;

// Time
double global_time;
double global_time_irr;
ULL NextRegTimeBlock;
int time_block;
double time_step;
ULL block_max;
double outputTimeStep;
double endTime;

double eta;

// Enzo to Nbody
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoTimeStep;


int FixNumNeighbor;
double InitialNeighborRadius;

// i/o
char* fname;
bool restart;
char* foutput;
double outputTime;
int outNum;
char *config_file;

FILE* binout;
FILE* mergerout;
#ifdef SEVN
FILE* SEVNout;
IO* sevnio = nullptr;
std::multimap<double, int> SEVNList;
#endif
FILE* workerout;

#ifdef PERFORMANCETRACE
Performance performance;
#endif

double E_binary;	// Total binary energy in a processor
double E_binary_SD; // Total slowdown binary energy in a processor
double E_merger; 	// Total merger energy in a processor
double E_PN;		// Total post-Newtonian energy in a processor // Experimental one.. It seems not working well by EW 2025.8.19

void DefaultGlobal() {

	/* Timesteps */
	endTime = 1;
	EnzoTimeStep   = endTime/1e10; // endTime should be Myr
	outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

	time_block = -30;
	block_max = static_cast<ULL>(pow(2, -time_block));
	time_step = std::pow(2,time_block);

	endTime = 0.0;
	outputTimeStep = 0.;

	global_time = 0.;
	outputTime = 0.;

	eta = 0.01;
	FixNumNeighbor = 100;
	InitialNeighborRadius = 0.011;

	E_binary = 0.0;
	E_binary_SD = 0.0;
	E_merger = 0.0;
	E_PN = 0.0;

#ifdef SEVN
	std::vector<std::string> args = {"empty", // Not used
		// "-myself", "/data/vinicius/NbodyPlus/SEVN",
		"-tables", "/home/vinicius/install/sevn_custom/tables/SEVNtracks_parsec_ov04_AGB", 
		//  "-tables", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_MIST_AGB",
		// "-tables_HE", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_parsec_pureHe36",
		// "-turn_WR_to_pureHe", "false",
		// "-snmode", "delayed",
		// "-Z", "0.0002",
		// "-spin", "0.0",
		// "-tini", "zams", 
		// "-tf", "end",
		// "-dtout", "events",
		"-xspinmode", "geneva",
		"-hardmode", "disabled", // binary hardening due to external perturbers
		// "-collmode", "disabled", // collision at periastron // If this is disabled, no common envelope or stellar merger in SEVN
		// "-circmode", "disabled", // orbital circularisation; default: circualise conseving the binary angular momentum
		"-tmode", "disabled" // equilibrium tides // Currently, we're not considering stellar rotation, so let's disable this
		// "-gwmode", "disabled"}; // Peters formula; this is considered in SDAR but let's turn this on
	};
	std::vector<char*> c_args;
	for (auto& arg : args) {
		c_args.push_back(&arg[0]);
	}

	sevnio = new IO;
	sevnio->load(c_args.size(), c_args.data());
#endif

}
