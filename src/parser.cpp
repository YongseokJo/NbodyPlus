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
				std::cerr << "Error: Missing argument for -c/--config option" << std::endl;
				std::exit(EXIT_FAILURE); // Error: Missing argument
			}

		} else {
			// Handle unrecognized arguments
			std::cerr << "Error: Unrecognized argument: " << arg << std::endl;
			std::exit(EXIT_FAILURE); // Error: Unrecognized argument
		}
	}

	EnzoTimeStep   = endTime/1e10; // endTime should be Myr
	outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

	if (MyRank == ROOT) {
		std::cout << "Configuration file: " << config_file <<  std::endl;
	}
}


/*
		} else if (arg == "-dtdump" || arg == "--timestep") {
			if (i + 1 < argc) {
				outputTimeStep = std::atof(argv[i + 1]);
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -t/--double option" << std::endl;
				return -2; // Error: Missing argument
			}

		} else if (arg == "-tend" || arg == "--end") {
			if (i + 1 < argc) {
				endTime = std::atof(argv[i + 1]);
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -tend/--double option" << std::endl;
				return -2; // Error: Missing argument
			}
		} else if (arg == "-d" || arg == "--dir") {
			if (i + 1 < argc) {
				foutput = argv[i + 1];
				i++; // Skip the next argument as it's already processed
			}
*/
