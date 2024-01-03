#include "Particle.h"
#include <cmath>

double Particle::getDistanceTo(Particle *particle) {
	double d=0; 
	for (int i=0; i<Dim; i++)
		d += std::pow(this->Position[i]-particle->Position[i],2);
	return std::sqrt(d); 
}

void Particle::setParticleInfo(double *data, int PID) {
	this->PID          = PID;
	this->Position[0]  = data[0];
	this->Position[1]  = data[1];
	this->Position[2]  = data[2];
	this->Velocity[0]  = data[3];
	this->Velocity[1]  = data[4];
	this->Velocity[2]  = data[5];
	this->Mass         = data[6];
	this->ParticleType = Star;
}
