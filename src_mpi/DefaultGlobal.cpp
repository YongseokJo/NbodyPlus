#include <iostream>
#include <stdio.h>
#include "global.h"
#ifdef SEVN
#include <map>
#endif

Particle *particles_original;
Particle *particles;
int *ActiveIndexToOriginalIndex;
int *ActiveIndexToOriginalIndex_orginal;

MPI_Win win;
MPI_Win win2;
MPI_Win win3;

MPI_Comm shared_comm;
int MyRank;
int NumberOfProcessor;
int NumberOfWorker;
int NumberOfCommunication;


GlobalVariable *global_variable;
GlobalVariable *global_variable_original;
int LastParticleIndex; // The last index of particle array
int NumberOfParticle; // The number of active particles
int NewPID;

// Task
int Task[NumberOfTask];

// Time
double global_time;
double global_time_irr;
ULL NextRegTimeBlock;
int time_block;
double time_step;
ULL block_max;
double outputTimeStep;
double endTime;

double binary_time;
double binary_time_prev;
ULL binary_block;

// Enzo to Nbody
Particle* FirstEnzoParticle;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoTimeStep;



// i/o
char* fname;
double inputTime;
bool restart;
char* foutput;
bool IsOutput;
double outputTime;
int outNum;

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

void DefaultGlobal() {

	/* Task initialization */
	//int Task[NumberOfTask];
	for (int i=0;i<NumberOfTask; i++) {
		Task[i] = i;
	}

	NumberOfCommunication = 0;

	/* Timesteps */
	endTime = 1;
	EnzoTimeStep   = endTime/1e10; // endTime should be Myr
	outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

	time_block = -30;
	block_max = static_cast<ULL>(pow(2, -time_block));
	time_step = std::pow(2,time_block);

	inputTime = 0.0;
	endTime = 0.0;
	outputTimeStep = 0.;

	global_time = 0.;
	outputTime = 0.;

#ifdef SEVN
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

	sevnio = new IO; // (SEVN Query) We should initialize sevnio only once. Here in Abyss, and somewhere else in Enzo by EW 2025.3.27
	sevnio->load(c_args.size(), c_args.data());
#endif

}
