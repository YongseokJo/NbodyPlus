/*
*  Purporse: Particle Class
*
*  Date    : 2024.01.02  by Yongseok Jo
*
*/


#include <vector>
#include <iostream>
#include "../defs.h"
//#include "../global.h"

class Binary;
class Particle
{
	private:

	// Variables for KS Regularization (117p in GNS by Aarseth, will be added)	
	//
	//

	public:
		int ParticleOrder;
		int PID;
		int ParticleType;
		REAL Mass;
		REAL InitialMass;
		REAL CreationTime;
		REAL DynamicalTime;
		REAL Velocity[Dim];
		REAL Position[Dim];
		REAL NewVelocity[Dim];
		REAL NewPosition[Dim];
		REAL PredTime;
		REAL PredTimeIrr;
		REAL PredTimeReg;
		REAL CurrentTimeIrr;
		REAL CurrentTimeReg;
		ULL CurrentBlockIrr;
		ULL CurrentBlockReg;
		REAL TimeStepIrr;
		REAL TimeStepReg;
		ULL TimeBlockIrr;
		ULL TimeBlockReg;
		int TimeLevelIrr;
		int TimeLevelReg;
		REAL PredPosition[Dim];
		REAL PredVelocity[Dim];
		REAL a_tot[Dim][HERMITE_ORDER];
		REAL a_reg[Dim][HERMITE_ORDER];
		REAL a_irr[Dim][HERMITE_ORDER];
		REAL BackgroundAcceleration[Dim];
		REAL LocalDensity;
		Particle* NextParticleInEnzo;
		Particle* NextParticleForComputation;
		Particle* BinaryPairParticle;
		Particle* BinaryParticleI;
		Particle* BinaryParticleJ;
		Binary* BinaryInfo;
		std::vector<Particle*> ACList;     // list of AC neighbor 
		int NumberOfAC; // number of neighbors
		REAL RadiusOfAC;
		bool isStarEvolution;
		bool isBinary; // check whether this is a member of the binary
		bool isCMptcl; // check if this particle is center-of-mass particle
		bool isErase;

		// Constructor
		Particle(void) {__initializer__();};
		Particle(REAL *data, int PID);
		void __initializer__(void) {
			//std::cout << "Constructor called" << std::endl;
			PID             = -1;
			ParticleOrder   = -1;
			Mass            = 0;
			InitialMass     = 0;
			NumberOfAC      = 0; // number of neighbors
			RadiusOfAC      = -1;
			ParticleType    = -9999;
			CurrentTimeIrr  = 0.; // consistent with actual current time
			CurrentTimeReg  = 0.;
			CurrentBlockIrr = 0; // consistent with actual current time
			CurrentBlockReg = 0;	
			PredTimeIrr     = 0;
			PredTimeReg     = 0;
			TimeStepIrr     = 0;
			TimeStepReg     = 0;
			TimeLevelIrr    = 0;
			TimeLevelReg    = 0;
			TimeBlockIrr    = 0;
			TimeBlockReg    = 0;
			LocalDensity    = 0;
			isStarEvolution = true;
			isBinary        = false;
			isCMptcl        = false;
			isErase         = false;
			for (int i=0; i<Dim; i++) {
				Velocity[i]     = 0.;
				Position[i]     = 0.;
				PredPosition[i] = 0.;
				PredVelocity[i] = 0.;
				BackgroundAcceleration[i] = 0.;
				for (int j=0; j<HERMITE_ORDER; j++) {
					a_tot[i][j] = 0.;
					a_reg[i][j] = 0.;
					a_irr[i][j] = 0.;
				}
			}
			NextParticleInEnzo         = nullptr;
			NextParticleForComputation = nullptr;
			BinaryPairParticle         = nullptr;
			BinaryParticleI            = nullptr;
			BinaryParticleJ            = nullptr;
			BinaryInfo                 = nullptr;
			ACList.clear();
		};

		void updateParticle(REAL mass, REAL *vel, REAL pos[], int particletype) {

			Mass = mass;
			ParticleType = particletype;

			for (int i=0; i<Dim; i++) {
				Velocity[i] = vel[i];
				Position[i] = pos[i];
			}
		};
		REAL getDistanceTo(Particle *particle);
		void setParticleInfo(REAL *data, int PID);
		void setParticleInfo(REAL *data, int PID, Particle* NextParticleInEnzo);
		void setParticleInfo(int *PID, REAL *Mass, REAL *CreationTime, REAL *DynamicalTime,
		REAL *Position[Dim], REAL *Velocity[Dim],
		REAL *BackgroundAcceleration[Dim], Particle* NextParticleInEnzo, int i);
		void setParticleInfo(int *PID, REAL *Mass, REAL *CreationTime, REAL *DynamicalTime,
		REAL *Position[Dim], REAL *Velocity[Dim],
		REAL *BackgroundAcceleration[Dim],  int ParticleType, Particle* NextParticleInEnzo, int i);
		void setParticleInfo(int *PID, REAL *BackgroundAcceleration[Dim],
		Particle* NextParticleInEnzo, int i);
		void setParticleInfo(REAL *Mass, REAL *BackgroundAcceleration[Dim],
		Particle* NextParticleInEnzo, int i);
		void initializeTimeStep();
		int getPID() {return PID;};
		void calculateIrrForce();
		void calculateRegAccelerationSecondOrder(std::vector<Particle*> &particle);
		void calculateRegAccelerationFourthOrder(std::vector<Particle*> &particle);

		void predictParticleSecondOrder(REAL time);
		void predictParticleSecondOrderIrr(REAL time);
		void correctParticleFourthOrder(REAL current_time, REAL next_time, REAL a[3][4]);

		void normalizeParticle();
		void calculateTimeStepIrr(REAL f[3][4], REAL df[3][4]);
		void calculateTimeStepReg();
		bool checkNeighborForEvolution();
		void updateEvolveParticle(std::vector<Particle*> &particle);
		void updateParticle();
		REAL evolveStarMass(REAL t1, REAL t2);
		void isKSCandidate();
		void convertBinaryCoordinatesToCartesian();
		void polynomialPrediction(REAL current_time);
		void UpdateRadius();
		void UpdateNeighbor(std::vector<Particle*> &particle);
		~Particle() {
			ACList.clear();
			// Deallocate memory
			ACList.shrink_to_fit();
			NextParticleInEnzo         = nullptr;
			NextParticleForComputation = nullptr;
			BinaryPairParticle         = nullptr;
			BinaryParticleI            = nullptr;
			BinaryParticleJ            = nullptr;
			BinaryInfo                 = nullptr;
			fprintf(stderr, "deleting particle, pid=%d\n", PID);
		};
};


