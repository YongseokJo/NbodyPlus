#include <iostream>
#include "global.h"

char* fname;
char* foutput;
REAL inputTime = 0.0;
REAL endTime = 0.0;
REAL outputTimeStep = 0.;

int Parser(int argc, char *argv[]) {
	for (int i = 1; i < argc; ++i) {

		std::string arg = argv[i];

		if (arg == "-h" || arg == "--help") {
			// Display help message
			std::cout << "Help message: Your program description here" << std::endl;

		} else if (arg == "-f" || arg == "--file") {
			// Check if the next argument exists
			if (i + 1 < argc) {
				fname = argv[i + 1];
				std::cout << "File name: " << fname << std::endl;
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -f/--file option" << std::endl;
				return -1; // Error: Missing argument
			}

		} else if (arg == "-dtdump" || arg == "--timestep") {
			if (i + 1 < argc) {
				outputTimeStep = std::atof(argv[i + 1]);
				std::cout << "Double value: " << outputTimeStep << std::endl;
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -t/--REAL option" << std::endl;
				return -2; // Error: Missing argument
			}

		} else if (arg == "-tend" || arg == "--end") {
			if (i + 1 < argc) {
				endTime = std::atof(argv[i + 1]);
				std::cout << "Double value: " << endTime << std::endl;
				i++; // Skip the next argument as it's already processed
			} else {
				std::cerr << "Error: Missing argument for -tend/--REAL option" << std::endl;
				return -2; // Error: Missing argument
			}
		} else if (arg == "-d" || arg == "--dir") {
			if (i + 1 < argc) {
				foutput = argv[i + 1];
				std::cout << "Output file name: " << fname << std::endl;
				i++; // Skip the next argument as it's already processed
			}
		} else {
			// Handle unrecognized arguments
			std::cerr << "Error: Unrecognized argument: " << arg << std::endl;
			return -11; // Error: Unrecognized argument
		}
	}
	return 1;
}
