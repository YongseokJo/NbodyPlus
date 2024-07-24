#include <vector>
#include <iostream>
#include "../defs.h"

class Particle;
class Binary
{
	private:


	public:

		// corresponding single particles and neighbor particles
		// overlapping variables with ptclCM are deleted. 
		// int NumberOfAC;     // number of neighbors - not currently in use
		// REAL RadiusOfAC;
		// std::vector<Particle*> ACList;     // list of AC neighbor - not currently in use

		Particle* ptclCM;
		bool isTerminate;
		bool isErase;

		// information of binary particles in cartesian coordinates
		REAL Mass;

		// REAL Position[Dim];
		// REAL Velocity[Dim];
		// REAL PredPosition[Dim];
		// REAL PredVelocity[Dim];
		REAL PredTime;
		REAL CurrentTime;  // this show how much the binary system has evolved
		REAL TimeStep;
		int TimeLevel;

		// REAL a_tot[Dim][HERMITE_ORDER];
		// REAL a_reg[Dim][HERMITE_ORDER];
		// REAL a_irr[Dim][HERMITE_ORDER];
		// REAL BackgroundAcceleration[Dim];


		// information of binary particles after Levi-civita transformation
		REAL u[4];
		REAL udot[4];
		REAL u_pred[4];
		REAL udot_pred[4];
		REAL u2dot[4];
		REAL u3dot[4];
		REAL u4dot[4];
		REAL u5dot[4];

		// Q = (L^T)F_perturb : Levi-civita transformation of perturbing force

		REAL Q[4];
		REAL Qdot[4];

		// r: distance between two binaries

		REAL r;
		REAL r0; // initial separation between binaries
		REAL rdot;
		REAL r2dot;
		REAL r3dot;
		REAL r4dot;
		REAL r5dot;

		// t: physical time. time derivatives with respect to tau
		// tdot = r, t2dot = rdot and so on

		REAL tdot;
		REAL t2dot;
		REAL t3dot;
		REAL t4dot;
		REAL t5dot;
		REAL t6dot;



		// h: energy of binary

		REAL h;
		REAL hdot;
		REAL h2dot;
		REAL h3dot;
		REAL h4dot;


		// tau: time unit for new coordinate system, tau = r*t

		REAL PredTau;
		REAL CurrentTau;
		REAL dTau;
		// REAL TauStep;

		// other variables

		REAL gamma;       // gamma: relative energy ratio of binary
		REAL a;       // semi-major axis of orbit

		// stumpff coefficient related variables

		REAL cn[stumpffN+1];
		REAL cn_4z[stumpffN+1];


		// Constructor
		Binary(void) {

			//std::cout << "Constructor called" << std::endl;

			isTerminate = false;
			isErase     = false;
			// NumberOfAC     = 0; // number of neighbors
			// RadiusOfAC     = InitialRadiusOfAC;
			// Mass           = 0;
			PredTime       = 0;
			CurrentTime    = 0;
			TimeStep       = 0;
			TimeLevel      = 0;

			PredTau       = 0;
			CurrentTau    = 0;
			dTau          = 9999;
			// TauStep       = 9999;

			gamma = 0;
			a     = 0;

			r0    = 0; // initial separation between binaries
			r     = 0;
			rdot  = 0;
			r2dot = 0;
			r3dot = 0;
			r4dot = 0;
			r5dot = 0;

			tdot = 0;
			t2dot = 0;
			t3dot = 0;
			t4dot = 0;
			t5dot = 0;
			t6dot = 0;


			h     = 0;
			hdot  = 0;
			h2dot = 0;
			h3dot = 0;
			h4dot = 0;

			// for (int i=0; i<Dim; i++) {

			// 	Position[i]     = 0;
			// 	Velocity[i]     = 0;
			// 	PredPosition[i] = 0;
			// 	PredVelocity[i] = 0;

			// 	BackgroundAcceleration[i] = 0;

			// 	for (int j=0; j<HERMITE_ORDER; j++) {
			// 		a_reg[i][j] = 0;
			// 		a_irr[i][j] = 0;
			// 		a_tot[i][j] = 0;
			// 	}
			// }  -> it is all saved in the CM particle information

			for (int i=0; i<4; i++){

				u[i]          = 0;
				udot[i]       = 0;
				u_pred[i]     = 0;
				udot_pred[i]  = 0;

				u2dot[i]      = 0;
				u3dot[i]      = 0;
				u4dot[i]      = 0;
				u5dot[i]      = 0;

				Q[i]           = 0;
				Qdot[i]        = 0;
			}

			for (int i=0; i<stumpffN+1; i++) {
				cn[i]    = 1.0;
				cn_4z[i] = 1.0;
			}
		}

		void InitializeBinary(REAL current_time);
		void getStumpffCoefficients(REAL z);
		void KSIntegration(REAL next_time, int &calnum);
		void predictBinary(REAL next_time);
		void IntegrateBinary(REAL next_time);

		~Binary() {};

};
