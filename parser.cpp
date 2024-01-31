#include <iostream>
#include "global.h"

char* fname;
double inputTime = 0.0;

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

		} else if (arg == "-t" || arg == "--time") {
            // Check if the next argument exists
            if (i + 1 < argc) {
                inputTime = std::atof(argv[i + 1]);
                std::cout << "Double value: " << inputTime << std::endl;
                i++; // Skip the next argument as it's already processed
            } else {
                std::cerr << "Error: Missing argument for -d/--double option" << std::endl;
                return -2; // Error: Missing argument
			}

		} else if (arg == "-dt" || arg == "--time") {
            // Check if the next argument exists
            if (i + 1 < argc) {
                inputTime = std::atof(argv[i + 1]);
                std::cout << "Double value: " << outputTimeStep << std::endl;
                i++; // Skip the next argument as it's already processed
            } else {
                std::cerr << "Error: Missing argument for -d/--double option" << std::endl;
                return -2; // Error: Missing argument
			}
		} else {
			// Handle unrecognized arguments
			std::cerr << "Error: Unrecognized argument: " << arg << std::endl;
			return -11; // Error: Unrecognized argument
		}
	}
	return 1;
}
