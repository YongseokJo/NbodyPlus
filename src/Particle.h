#include "parameters.h"


class Particle
{
	private:
		int PID;
		float Mass;
		float Velocity[Dim];
		float Position[Dim];
		float CurrentTime[2]; // 0 for irr and 1 for reg
		float TimeStep[2]; // 0 for irr and 1 for reg
		float PositionPred[Dim];
		float VelocityPred[Dim];
		float Force[Dim];
		float ForceDot[Dim];
		float FIrr[Dim][HERMITE_ORDER];
		float FReg[Dim][HERMITE_ORDER];
		int *ACList;     // list of AC neighbor 
		int NumberOfAC; // number of neighbors
		float RadiusOfAC;	 
		int NextParticle;
		int ParticleType;

	// Variables for KS Regularization (117p in GNS by Aarseth, will be added)	
	//
	//
	//
	

}
