#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "Particle.h"
#include "global.h"


int getLineNumber();
const int NUM_COLUMNS = 7; // Define the number of columns


int readData(std::vector<Particle*> &particle) {
	std::ifstream inputFile("data.dat");

	if (!inputFile) {
		std::cerr << "Error: Could not open the file." << std::endl;
		return FAIL;
	}


	int NumParticle;
	NumParticle = getLineNumber();
	NNB = NumParticle;

	// Declaration
	Particle *particle_temp;	
	particle_temp = new Particle[NumParticle];
	double** data = new double*[NumParticle];

	for (int i = 0; i < NumParticle; ++i) {
		data[i] = new double[NUM_COLUMNS];
	}

	// Initialization
	for (int i = 0; i < NumParticle; ++i) {
		for (int j = 0; j < NUM_COLUMNS; ++j) {
			data[i][j] = 0;
		}
	}


	int row = 0;

	std::string line;
	while (std::getline(inputFile, line) && row < NumParticle) { // Read lines from the file
		std::istringstream iss(line); // Create a stringstream for each line

		double value;
		int col = 0;
		while (iss >> value && col < NUM_COLUMNS) { // Read values from the stringstream
			data[row][col] = value;
			++col;
		}
		particle_temp[row].setParticleInfo(data[row], row);
		particle.push_back(&particle_temp[row]);
		++row;
	}




	/***
	for (int col = 0; col < NUM_COLUMNS; ++col) {
		std::cout << "Column " << col + 1 << " values: ";
		for (int r = 0; r < row; ++r) {
			std::cout << columnData[col][r] << ' ';
		}
		std::cout << std::endl;
	}
	***/

	inputFile.close();




	// Deallocate memory 
	for (int i = 0; i < NUM_COLUMNS; ++i) {
		delete[] data[i];
	}
	delete[] data;



	return SUCCESS;
}


int getLineNumber() {
	    std::ifstream inputFile("data.dat"); // Open the file

			if (!inputFile) {
				std::cerr << "Error: Could not open the file." << std::endl;
				return 1;
			}

			int lineCount = 0;
			std::string line;
			while (std::getline(inputFile, line)) { // Read lines from the file
				lineCount++;
			}

			std::cout << "Number of lines in the file: " << lineCount << std::endl;

			inputFile.close(); // Close the file

			return lineCount;
}

int WriteData() {
	return SUCCESS;
}


