#include <iostream>
#include "defs.h"
#include "global.h"


class Particle
{
	private:
		int PID;
		int ParticleType;

		// Variables for KS Regularization (117p in GNS by Aarseth, will be added)	
		//
		//
		//
		//


	public:
		double Mass;
		double Velocity[Dim];
		double Position[Dim];
		double PredTimeIrr;
		double PredTimeReg;
		double CurrentTimeIrr;
		double CurrentTimeReg;
		double TimeStepIrr;
		double TimeStepReg;
		double PredPosition[Dim];
		double PredVelocity[Dim];
		double Force[Dim];
		double ForceDot[Dim];
		double FIrr[Dim];
		double FReg[Dim];
		double FIrrDot[Dim];
		double FRegDot[Dim];
		double dFIrr[Dim][HERMITE_ORDER];
		double dFReg[Dim][HERMITE_ORDER];
		int NextParticle;
		std::vector<Particle*> ACList;     // list of AC neighbor 
		int NumberOfAC; // number of neighbors
		double RadiusOfAC;

		// Constructor
		Particle(void) {
			//std::cout << "Constructor called" << std::endl;
			Mass = 0;
			NumberOfAC = 0; // number of neighbors
			RadiusOfAC = 0.;
			NextParticle = 0;
			ParticleType = -9999;
			CurrentTimeIrr = 0;
			CurrentTimeReg = 0;
			PredTimeIrr = 0;
			PredTimeReg = 0;
			TimeStepIrr = 0;
			TimeStepReg = 0;
			for (int i=0; i<Dim; i++) {
				Velocity[i] = 0;
				Position[i] = 0;
				PredPosition[i] = 0;
				PredVelocity[i] = 0;
				Force[i] = 0;
				ForceDot[i] = 0;
				FIrr[i] = 0;
				FReg[i] = 0;
			  FIrrDot[i]=0;
		    FRegDot[i]=0;
				for (int j=0; j<HERMITE_ORDER; j++) {
					dFIrr[i][j] = 0;
					dFReg[i][j] = 0;
				}
			}
		}

		void updateParticle(double mass, double *vel, double pos[], int particletype) {

			Mass = mass;
			ParticleType = particletype;

			for (int i=0; i<Dim; i++) {
				Velocity[i] = vel[i];
				Position[i] = pos[i];
			}
		}
		double getDistanceTo(Particle *particle);
		void setParticleInfo(double *data, int PID);
		void initializeTimeStep();
		int getPID() {return PID;};
		void calculateIrrForce(double next_time);
		void predictPosAndVel(double next_time);
};

