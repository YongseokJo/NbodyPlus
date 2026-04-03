#ifdef SEVN
#include "sevn.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include "def.h"
#include "particle.h"
#include "global_state.h"
#include "global.h"
#include "mcluster_runner.h"
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

	bin_output_file = fopen("binary_output.txt", "w");
	fprintf(bin_output_file, "Starting ABYSS - Binary OUTPUT\n");
	fflush(bin_output_file);
	merger_output_file = fopen("merger_output.txt", "w");
	fprintf(merger_output_file, "Starting ABYSS - Merger OUTPUT\n");
	fflush(merger_output_file);
#ifdef SEVN
	sevn_output_file = fopen("SEVN_output.txt", "w");
	fprintf(sevn_output_file, "Starting ABYSS - SEVN OUTPUT\n");
	fflush(sevn_output_file);
#endif

	/* MPI Initialization */
	initializeMPI(argc, argv);

#ifdef CUDA
	if (my_rank == ROOT) {
		OpenDevice();
		cudaDeviceSynchronize();
	}
#endif

	/* Input options */
	Parser(argc, argv);
	readParameterFile();

	// Phase 27: McLuster IC generation (if configured)
	if (mcluster_config.has_mcluster_section) {
		bool mcluster_success = true;
		std::string abyss_ic = "mcluster_abyss.dat";

		if (my_rank == ROOT) {
			std::cout << "\n=== McLuster IC Generation ===" << std::endl;

			// Check if IC file was pre-generated (by separate McLuster job)
			struct stat ic_stat;
			bool ic_exists = (stat(abyss_ic.c_str(), &ic_stat) == 0 && ic_stat.st_size > 0);

			if (ic_exists) {
				// Use pre-generated IC file
				std::cout << "Found pre-generated IC file: " << abyss_ic << std::endl;
				std::cout << "Skipping McLuster execution" << std::endl;

				// Update fname to use the pre-generated IC file
				static std::string generated_ic_path = abyss_ic;
				fname = const_cast<char*>(generated_ic_path.c_str());
				std::cout << "IC file ready: " << fname << std::endl;
			} else {
				// Run McLuster to generate IC
				std::vector<std::string> args = buildMclusterArgs(mcluster_config);

				std::cout << "Running: " << MCLUSTER_BINARY;
				for (const auto& arg : args) {
					std::cout << " " << arg;
				}
				std::cout << std::endl;

				RunResult result = runMclusterSubprocess(MCLUSTER_BINARY, args);

				if (!result.success) {
					std::cerr << "McLuster failed with exit code: " << result.exit_code << std::endl;
					if (!result.stderr_content.empty()) {
						std::cerr << "Error output:\n" << result.stderr_content << std::endl;
					}
					mcluster_success = false;
				} else {
					// Validate output
					std::string mcluster_output = MCLUSTER_OUTPUT_BASE + ".txt";
					if (!validateMclusterOutput(mcluster_output, 0)) {
						std::cerr << "McLuster output validation failed" << std::endl;
						mcluster_success = false;
					} else {
						// Transform to ABYSS format
						if (!transformMclusterOutput(mcluster_output, abyss_ic)) {
							std::cerr << "Failed to transform McLuster output" << std::endl;
							mcluster_success = false;
						} else {
							// Update fname to use generated IC file
							static std::string generated_ic_path = abyss_ic;
							fname = const_cast<char*>(generated_ic_path.c_str());
							std::cout << "IC file ready: " << fname << std::endl;
						}
					}
				}
			}

			std::cout << "================================\n" << std::endl;
		}

		// Broadcast success/failure to all ranks
		int success_flag = mcluster_success ? 1 : 0;
		MPI_Bcast(&success_flag, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

		if (!success_flag) {
			if (my_rank == ROOT) {
				std::cerr << "Aborting due to McLuster failure" << std::endl;
			}
			MPI_Finalize();
			return 1;
		}

		// Handle generate_only mode
		if (mcluster_config.generate_only) {
			if (my_rank == ROOT) {
				std::cout << "generate_only=true: IC generation complete, exiting." << std::endl;
			}

			// Clean shutdown - free MPI resources that were allocated
			MPI_Win_free(&win);
			MPI_Win_free(&win2);
			MPI_Win_free(&win3);
			MPI_Win_free(&win4);
			particle_data.deallocate_shared();
			MPI_Comm_free(&shared_comm);
			MPI_Type_free(&queue_type_mpi);
			MPI_Type_free(&iparticle_type_mpi);
			MPI_Type_free(&jparticle_type_mpi);
			MPI_Finalize();
			return 0;
		}
	}

	// Write Particles
	if (my_rank == ROOT && !readData())
		fprintf(stderr, "Read Data Failed!\n");
	

	if (my_rank == ROOT) {
		RootRoutines();
	} else {
		std::string filename = "worker_output_" + std::to_string(my_rank) + ".txt";
		worker_output_file = fopen(filename.c_str(), "w");
		fprintf(worker_output_file, "Starting ABYSS - WORKER OUTPUT\n");
		fflush(worker_output_file);
		
		WorkerRoutines();
	}

	// Finalize the window and MPI environment
	MPI_Win_free(&win);
	MPI_Win_free(&win2);
	MPI_Win_free(&win3);
	MPI_Win_free(&win4);

	// Deallocate SoA particle data
	particle_data.deallocate_shared();

	MPI_Comm_free(&shared_comm);

	MPI_Type_free(&queue_type_mpi);
    MPI_Type_free(&iparticle_type_mpi);
    MPI_Type_free(&jparticle_type_mpi);

	MPI_Finalize();

	return 0;
}

