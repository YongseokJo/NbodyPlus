#include <iostream>

using namespace std;
void start(Particle *particle);


int main() {
	cout << "This Nbody+." << endl;



    Particle *particle;
    particle = new Particle[NNB];

	start(particle);
	return 0;



}
