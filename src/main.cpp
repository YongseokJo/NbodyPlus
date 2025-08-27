#ifdef SEVN
#include "sevn.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "global.h"
#include <mpi.h>
#ifdef CUDA
#include <cuda_runtime.h>
#include "cuda/cuda_functions.h"
#endif


void broadcastFromRoot(int &data);
void DefaultGlobal();
void initializeMPI(int argc, char *argv[]);
void WorkerRoutines();
void RootRoutines();
void Parser(int argc, char *argv[]);
bool readData();
void readParameterFile();


int main(int argc, char *argv[]) {

	/* Initialize global variables */
	DefaultGlobal();

	binout = fopen("binary_output.txt", "w");
	fprintf(binout, "Starting nbody - Binary OUTPUT\n");
	fflush(binout);
	mergerout = fopen("merger_output.txt", "w");
	fprintf(mergerout, "Starting nbody - Merger OUTPUT\n");
	fflush(mergerout);
#ifdef SEVN
	SEVNout = fopen("SEVN_output.txt", "w");
	fprintf(SEVNout, "Starting nbody - SEVN OUTPUT\n");
	fflush(SEVNout);
#endif

	/* MPI Initialization */
	initializeMPI(argc, argv);

#ifdef CUDA
	if (MyRank == ROOT) {
		OpenDevice();
		cudaDeviceSynchronize();
	}
#endif

	/*
	// Insert this function definition at the top of your code after the include directives.
	char hostname[256];
	gethostname(hostname, sizeof(hostname));

	// Insert this code right after the  MPI initialization routines (though not a mandatory requirement 
	// to add there only). Please make a judgement based on your code.
	// Retrieve process ID and hostname
	pid_t pid = getpid();

	volatile int i = 0;
	while (0 == i)
	{
		std::cout << "My rank = " << MyRank << " PID = " << pid << " running on Host = " << hostname << " in sleep " << std::endl;
		sleep(5);
	}
	*/

	/* Input options */
	Parser(argc, argv);
	readParameterFile();

	// Write Particles
	if (MyRank == ROOT && !readData())
		fprintf(stderr, "Read Data Failed!\n");
	

	if (MyRank == ROOT) {
		global_variable->LastParticleIndex = LastParticleIndex;

		RootRoutines();
	} else {
		// /* // by EW 2025.1.27
		std::string filename = "worker_output_" + std::to_string(MyRank) + ".txt";
		workerout = fopen(filename.c_str(), "w");
		fprintf(workerout, "Starting nbody - WORKER OUTPUT\n");
		fflush(workerout);
		// */
		
		WorkerRoutines();
	}

	// Finalize the window and MPI environment
	MPI_Win_free(&win);
	MPI_Finalize();
	return 0;
}

