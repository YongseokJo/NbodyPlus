#include <iostream>
#include "global.h"


void Parser(int argc, char *argv[]) {
	for (int i = 1; i < argc; ++i) {

		std::string arg = argv[i];

		if (arg == "-h" || arg == "--help") {
			// Display help message
			std::cout << "Help message: Your program description here" << std::endl;

		} else if (arg == "-c" || arg == "--config") {
			// Check if the next argument exists
			if (i + 1 < argc) {
				config_file = argv[i + 1];
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "TASK_ERROR: Missing argument for -c/--config option" << std::endl;
				std::exit(EXIT_FAILURE); // TASK_ERROR: Missing argument
			}

		} else {
			// Handle unrecognized arguments
			std::cerr << "TASK_ERROR: Unrecognized argument: " << arg << std::endl;
			std::exit(EXIT_FAILURE); // TASK_ERROR: Unrecognized argument
		}
	}

	enzo_time_step   = end_time/1e10; // end_time should be Myr
	output_time_step = output_time_step/end_time; // end_time should be Myr

	if (my_rank == ROOT) {
		std::cout << "Configuration file: " << config_file <<  std::endl;
	}
}


/*
		} else if (arg == "-dtdump" || arg == "--timestep") {
			if (i + 1 < argc) {
				output_time_step = std::atof(argv[i + 1]);
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "TASK_ERROR: Missing argument for -t/--double option" << std::endl;
				return -2; // TASK_ERROR: Missing argument
			}

		} else if (arg == "-tend" || arg == "--end") {
			if (i + 1 < argc) {
				end_time = std::atof(argv[i + 1]);
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "TASK_ERROR: Missing argument for -tend/--double option" << std::endl;
				return -2; // TASK_ERROR: Missing argument
			}
		} else if (arg == "-d" || arg == "--dir") {
			if (i + 1 < argc) {
				foutput = argv[i + 1];
				i++; // Skip the next argument as it's already processed
			}
*/
