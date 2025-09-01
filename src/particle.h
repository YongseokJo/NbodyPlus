#ifndef PARTICLE_H
#define PARTICLE_H

#include "def.h"
#include <cmath>
#include "cstring"
#include <vector>
#include "cuda/cuda_defs.h"
#include "GlobalVariable.h"

// SDAR
#include "Common/Float.h"
#include <iostream>
#include <iomanip>
#ifdef SEVN
#include "star.h" // Eunwoo added for SEVN
#ifdef SEVN_BINARY
#include "binstar.h"
#endif
#endif

extern double InitialNeighborRadius;
extern double EnzoTimeStep;
extern GlobalVariable *global_variable;
struct Particle;
extern Particle *particles;
extern double time_step;

enum class BinaryInterruptState:int {
	none = 0, 
	form = 1, 
	exchange = 2, 
	collisioncandidate = 3,
	collision = 4,
	threebody = 5,	// particles from 3body interaction // added by EW 2025.2.18
	manybody = 6,	// many-body (>3) group terminated // added by EW 2025.1.19
	kicked = 7,		// kicked or exploded as PISN during stellar evolution // added by EW 2025.1.19
	terminated = 8,	// CM particle only; terminated // added by EW 2025.1.20
	merger = 9		// CM particle only; merger inside // added by EW 2025.1.20
};
#define BINARY_STATE_ID_SHIFT 4 // Eunwoo
#define BINARY_INTERRUPT_STATE_MASKER 0xF // Eunwoo

struct Group;
struct Particle {
	int PID;
	int ParticleIndex; // (Query) please change it to index when you get a chance // (Answer) done!
	int ParticleType;
	double Position[3];      // Position in 3D space (x, y, z)
	double Velocity[3];      // Velocity in 3D space (vx, vy, vz)
	double Mass;             // Mass of the particle
	double a_tot[Dim][HERMITE_ORDER];
	double a_reg[Dim][HERMITE_ORDER];
	double a_irr[Dim][HERMITE_ORDER];
	int    NumberOfNeighbor;    // 
	int    NewNumberOfNeighbor;    // 

	double CurrentTimeIrr;
	double CurrentTimeReg;
	ULL CurrentBlockIrr;
	ULL NewCurrentBlockIrr;
	ULL CurrentBlockReg;
	ULL NextBlockIrr;
	double TimeStepIrr;
	double TimeStepReg;
	ULL TimeBlockIrr;
	ULL TimeBlockReg;
	int TimeLevelIrr;
	int TimeLevelReg;
	double NewPosition[Dim]; // this might be initialized in NewFBInitialization by EW 2025.1.7
	double NewVelocity[Dim]; // this might be initialized in NewFBInitialization by EW 2025.1.7
	double RadiusOfNeighbor; //this should be squared of it. // this might be initialized in NewFBInitialization by EW 2025.1.7
	double BackgroundAcceleration[Dim]; // this might be initialized in NewFBInitialization by EW 2025.1.7

	// For SDAR
	bool isActive;
	bool isUpdateToDate;
	double radius; // used in SEVN too
	double dm; // Stellar mass which will be distributed to nearby gas cells // used in SEVN too
	double time_check; // time to check next interrupt
	long long int binary_state; // contain two parts, low bits (first BINARY_STATE_ID_SHIFT bits) is binary interrupt state and high bits are pair ID
	double a_spin[3]; // dimensionless spin parameter a
	Group* GroupInfo;
	bool isCMptcl; // do we need this? (Query) we can simply use if GroupInfo == nullptr right?
	int CMPtclIndex; // added for write_out_group function by EW 2025.1.6
	int Members[10]; // ParticleIndex of group members; only used for cm ptcls by EW 2025.1.30
	int NumberOfMember; // Number of group members; only used for cm ptcls by EW 2025.1.30

#ifdef SEVN
	StarSEVN* StellarEvolution;
#ifdef SEVN_BINARY
	Binstar* BinaryEvolution;
#endif
	double FormationTime; // Myr // for restart
	double WorldTime; // Myr // FormationTime + EvolutionTime
#endif

	Particle() {
		Position[0] = Position[1] = Position[2] = 0.0;
		Velocity[0] = Velocity[1] = Velocity[2] = 0.0;
		PID             = -1;
		Mass            = 0;
		RadiusOfNeighbor= -1;
		NumberOfNeighbor= 0;
		NewNumberOfNeighbor= 0;
		ParticleType    = NO_FEEDBACK_STAR;
		CurrentTimeIrr  = 0.; // consistent with actual current time
		CurrentTimeReg  = 0.;
		CurrentBlockIrr = 0; // consistent with actual current time
		CurrentBlockReg = 0;
		NextBlockIrr    = 0;
		TimeStepIrr     = 0;
		TimeStepReg     = 0;
		TimeLevelIrr    = 0;
		TimeLevelReg    = 0;
		TimeBlockIrr    = 0;
		TimeBlockReg    = 0;
		for (int i=0; i<Dim; i++) {
			Velocity[i]     = 0.;
			Position[i]     = 0.;
			NewPosition[i] = 0.;
			NewVelocity[i] = 0.;
			BackgroundAcceleration[i] = 0.;
			for (int j=0; j<HERMITE_ORDER; j++) {
				a_tot[i][j] = 0.;
				a_reg[i][j] = 0.;
				a_irr[i][j] = 0.;
			}
			a_spin[i] = 0.;
		}
		isActive = true;
		ParticleIndex = -1;
		radius = 0.;
		dm = 0.0;
		time_check = NUMERIC_FLOAT_MAX;
		setBinaryInterruptState(BinaryInterruptState::none);
		GroupInfo = nullptr;
		isCMptcl = false;
		isUpdateToDate = true;
		CMPtclIndex = -1;
		NumberOfMember = 0;
#ifdef SEVN
		StellarEvolution = nullptr;
#ifdef SEVN_BINARY
		BinaryEvolution = nullptr;
#endif
		FormationTime = 0.0; // Myr
		WorldTime = 0.0; // Myr
#endif
	}

	/*
	Particle(double x, double y, double z, double vx, double vy, double vz, double m, double q)
		: Mass(m) {
			Position[0] = x; Position[1] = y; Position[2] = z;
			Velocity[0] = vx; Velocity[1] = vy; Velocity[2] = vz;
		}
		*/

	void initialize(double *data, int PID) {
		this->PID          = PID;
		this->Position[0]  = data[0];
		this->Position[1]  = data[1];
		this->Position[2]  = data[2];
		this->Velocity[0]  = data[3];
		this->Velocity[1]  = data[4];
		this->Velocity[2]  = data[5];
		this->Mass         = data[6];
		this->ParticleType = NO_FEEDBACK_STAR;
		//this->NextParticleInEnzo = NextParticleInEnzo;
		this->CurrentTimeReg		= 0;
		this->CurrentTimeIrr		= 0;
		this->RadiusOfNeighbor		= InitialNeighborRadius*InitialNeighborRadius;
		//this->RadiusOfNeighbor	= 1;
		//this->RadiusOfNeighbor	= 0;
		this->NumberOfNeighbor		= 0;
		this->NewNumberOfNeighbor	= 0;

		this->isActive				= true;
		this->ParticleIndex			= PID;
		this->dm = 0.0;
		this->time_check = NUMERIC_FLOAT_MAX;
		this->setBinaryInterruptState(BinaryInterruptState::none);
		this->GroupInfo = nullptr;
		this->isCMptcl = false;
		this->a_spin[0] = 0.;
		this->a_spin[1] = 0.;
		this->a_spin[2] = 0.;
		this->CMPtclIndex = -1;
		this->isUpdateToDate = true;
		this->NumberOfMember = 0;

#ifndef SEVN
		this->ParticleType = NO_FEEDBACK_STAR;
		this->radius = 2.25461e-8/position_unit*pow(this->Mass*1e9, 1./3); // stellar radius in code unit
		/*
		if (this->Mass*1e9 > 8) {
			this->ParticleType = Blackhole+SingleStar;
			this->radius = 2*this->Mass*1e9/mass_unit/pow(299752.458/(velocity_unit/yr*pc/1e5), 2); // Schwartzshild radius in code unit
			// initialBHspin(this);
		}
		else {
			this->ParticleType = NormalStar+SingleStar;
			this->radius = 2.25461e-8/position_unit*pow(this->Mass*1e9, 1./3); // stellar radius in code unit
		}
		*/
#endif
	}

	// Clear function // Simplified by EW 2025.1.7
    void clear() {
        PID = -1;
		ParticleIndex = -1;

        NumberOfNeighbor = 0;
        NewNumberOfNeighbor = 0;

		isUpdateToDate = true;
        isActive = false;
		GroupInfo = nullptr;
		isCMptcl = false;
		CMPtclIndex = -1;
		NumberOfMember = 0;
		setBinaryInterruptState(BinaryInterruptState::none);
		ParticleType = NO_FEEDBACK_STAR;

		a_spin[0] = 0.;
		a_spin[1] = 0.;
		a_spin[2] = 0.;
    }

	void normalizeParticle() {
		// pc to computing unit, km/s to computing unit
		this->Mass *= 1e9;
		this->Mass /= mass_unit;
		for (int dim=0; dim<Dim; dim++) {
			this->Position[dim] *= 1000; // kpc to pc
			this->Position[dim] /= position_unit;
			this->Velocity[dim] *= 1e5*yr/pc; // km/s to pc/yr
			this->Velocity[dim] /= velocity_unit;
		}
	}



	// void updateParticle(); 
	void updateParticle() { // inline function by EW 2025.3.3 to reduce function call time

		if (this->isCMptcl) {
			Particle* members;
			for (int i = 0; i < this->NumberOfMember; i++) {
				members = &particles[this->Members[i]];
				for (int dim = 0; dim < Dim; dim++) {
					members->Position[dim] += this->NewPosition[dim] - this->Position[dim];
					members->Velocity[dim] += this->NewVelocity[dim] - this->Velocity[dim];
				}
			}
		}
	
		for (int dim=0; dim<Dim; dim++) {
			this->Position[dim] = this->NewPosition[dim];
			this->Velocity[dim] = this->NewVelocity[dim];
		}
	}

	/*
	void getAcceleration(const double pos[], const double vel[]) {
		double dx[Dim], dr2;
		for (int dim=0; dim<Dim; dim++) {
			dx[dim] = pos[dim] - this->Position[dim];
			dr2 += dx[dim]*dx[dim];
		}

		for (int dim=0; dim<Dim; dim++) {
			a_reg[dim] += dx[dim]/dr2/sqrt(dr2);
		}
	}
	*/

	// inline function by EW 2025.7.15 to reduce function call time
	template <typename T>
	// void predictParticleSecondOrder(double dt, T pos[], T vel[]);
    void predictParticleSecondOrder(double dt, T pos[], T vel[]) {
        // Doubling check
        // temporary variables for calculation
    
        // only predict the positions if necessary
        // how about using polynomial correction here?
        
        dt = dt*EnzoTimeStep;
    
        if (dt == 0) {
            for (int dim=0; dim<Dim; dim++) {
                pos[dim] = static_cast<T>(Position[dim]);
                vel[dim] = static_cast<T>(Velocity[dim]);
            }
        }
        else {
            for (int dim=0; dim<Dim; dim++) {
                pos[dim] = static_cast<T>(((a_tot[dim][1]*dt/3 + a_tot[dim][0])*dt/2 + Velocity[dim])*dt + Position[dim]);
                vel[dim] = static_cast<T>((a_tot[dim][1]*dt/2 + a_tot[dim][0])*dt   + Velocity[dim]);
            }
        }
        return;
    }

	void predictParticleSecondOrder(double dt, std::vector<Jparticle>& jparticles, std::vector<Iparticle>& iparticles, std::vector<int>& localRegularList) {

		dt = dt*EnzoTimeStep;

		Iparticle iptcl;
		Jparticle jptcl;

		if (dt == 0) {
			jptcl.posx		= static_cast<CUDA_REAL>(Position[0]);
			jptcl.posy		= static_cast<CUDA_REAL>(Position[1]);
			jptcl.posz		= static_cast<CUDA_REAL>(Position[2]);
			jptcl.velx		= static_cast<CUDA_REAL>(Velocity[0]);
			jptcl.vely		= static_cast<CUDA_REAL>(Velocity[1]);
			jptcl.velz		= static_cast<CUDA_REAL>(Velocity[2]);
			jptcl.mass		= static_cast<CUDA_REAL>(Mass);
			jptcl.index		= ParticleIndex;
			jparticles.push_back(jptcl);

		} else {
			jptcl.posx		= static_cast<CUDA_REAL>(((a_tot[0][1]*dt/3 + a_tot[0][0])*dt/2 + Velocity[0])*dt + Position[0]);
			jptcl.posy		= static_cast<CUDA_REAL>(((a_tot[1][1]*dt/3 + a_tot[1][0])*dt/2 + Velocity[1])*dt + Position[1]);
			jptcl.posz		= static_cast<CUDA_REAL>(((a_tot[2][1]*dt/3 + a_tot[2][0])*dt/2 + Velocity[2])*dt + Position[2]);
			jptcl.velx		= static_cast<CUDA_REAL>((a_tot[0][1]*dt/2 + a_tot[0][0])*dt   + Velocity[0]);
			jptcl.vely		= static_cast<CUDA_REAL>((a_tot[1][1]*dt/2 + a_tot[1][0])*dt   + Velocity[1]);
			jptcl.velz		= static_cast<CUDA_REAL>((a_tot[2][1]*dt/2 + a_tot[2][0])*dt   + Velocity[2]);
			jptcl.mass		= static_cast<CUDA_REAL>(Mass);
			jptcl.index		= ParticleIndex;
			jparticles.push_back(jptcl);
		}
		if (CurrentBlockReg + TimeBlockReg == global_variable->NextRegTimeBlock) {
			iptcl.posx		= jptcl.posx;
			iptcl.posy		= jptcl.posy;
			iptcl.posz		= jptcl.posz;
			iptcl.velx		= jptcl.velx;
			iptcl.vely		= jptcl.vely;
			iptcl.velz		= jptcl.velz;
			iptcl.r2		= static_cast<CUDA_REAL>(RadiusOfNeighbor);
			iptcl.dtr		= static_cast<CUDA_REAL>(TimeBlockReg*time_step*EnzoTimeStep);
			iparticles.push_back(iptcl);
			localRegularList.push_back(ParticleIndex);
		}
	}

	void correctParticleFourthOrder(double dt, double pos[], double vel[], double a[3][4]);


	/*
	void update_timestep() {
		double acc = mag(acceleration);
		double vel = mag(Velocity);

		time_step = eta*sqrt(std::abs(vel/acc));
	}
	*/

	void updateRadius();

	//void initializeNeighbor();
	//void initializeAcceleration();
	void initializeTimeStep();

	void computeAccelerationIrr();
	void computeAccelerationReg();


	void calculateTimeStepIrr();
	void calculateTimeStepIrr2();
	void calculateTimeStepReg();

	void updateRegularParticleCuda();

	// SDAR
	void checkNewGroup();
	void checkNewGroup2();
	void checkNewGroup3();
	void checkNewGroup4();

	//! save pair id in binary_state with shift bit size of BINARY_STATE_ID_SHIFT
	void setBinaryPairID(const int _id) {
		binary_state = (binary_state&BINARY_INTERRUPT_STATE_MASKER) | (_id<<BINARY_STATE_ID_SHIFT);
	}

	//! save binary interrupt state in the first  BINARY_STATE_ID_SHIFT bit in binary_state
	void setBinaryInterruptState(const BinaryInterruptState _state) {
		binary_state = ((binary_state>>BINARY_STATE_ID_SHIFT)<<BINARY_STATE_ID_SHIFT) | int(_state);
	}

	//! get binary interrupt state from binary_state
	BinaryInterruptState getBinaryInterruptState() const {
		return static_cast<BinaryInterruptState>(binary_state&BINARY_INTERRUPT_STATE_MASKER);
	}

	//! get pair ID from binary_state 
	int getBinaryPairID() const {
		return (binary_state>>BINARY_STATE_ID_SHIFT);
	}

	double* getPos() {
		return Position;
	}

	double* getVel() {
		return Velocity;
	}

	static void printColumnTitle(std::ostream & _fout, const int _width=20) {
		_fout<<std::setw(_width)<<"mass"
			<<std::setw(_width)<<"pos.x"
			<<std::setw(_width)<<"pos.y"
			<<std::setw(_width)<<"pos.z"
			<<std::setw(_width)<<"vel.x"
			<<std::setw(_width)<<"vel.y"
			<<std::setw(_width)<<"vel.z"
			<<std::setw(_width)<<"radius"
			<<std::setw(_width)<<"id";
	}

	//! print data of class members using column style (required)
	/*! print data of class members in one line for column style. Notice no newline is printed at the end
	@param[out] _fout: std::ostream output object
	@param[in] _width: print width (defaulted 20)
	*/

	void printColumn(std::ostream & _fout, const int _width=20){
		_fout<<std::setw(_width)<<Mass
			<<std::setw(_width)<<Position[0]
			<<std::setw(_width)<<Position[1]
			<<std::setw(_width)<<Position[2]
			<<std::setw(_width)<<Velocity[0]
			<<std::setw(_width)<<Velocity[1]
			<<std::setw(_width)<<Velocity[2]
			<<std::setw(_width)<<radius
			<<std::setw(_width)<<PID;
	}

	// made by EW 2025.1.6
	/*
	void copyNewNeighbor(Particle* ptcl) {
		this->NewNumberOfNeighbor = ptcl->NewNumberOfNeighbor;
		std::memcpy(this->NewNeighbors, ptcl->NewNeighbors, sizeof(int)*ptcl->NewNumberOfNeighbor);
	}
	*/
	
#define NO_PRINT_FULL_ACC
	void printParticleInfo(FILE* file) {
		fprintf(file, "PID: %d, ParticleIndex: %d, ParticleType: %d\n", PID, ParticleIndex, ParticleType);
		fprintf(file, "Position (pc):%e, %e, %e\n", Position[0]*position_unit, Position[1]*position_unit, Position[2]*position_unit);
		fprintf(file, "Velocity (km/s): %e, %e, %e\n", Velocity[0]*velocity_unit/yr*pc/1e5, Velocity[1]*velocity_unit/yr*pc/1e5, Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(file, "Mass (Msun): %e\n", Mass*mass_unit);
		fprintf(file, "NumNeighbor: %d, ACRadius (pc): %e\n", NumberOfNeighbor, sqrt(RadiusOfNeighbor)*position_unit);
		fprintf(file, "atot (0): %e, %e, %e\n", a_tot[0][0], a_tot[1][0], a_tot[2][0]);
#ifdef PRINT_FULL_ACC
		fprintf(file, "atot (1): %e, %e, %e\n", a_tot[0][1], a_tot[1][1], a_tot[2][1]);
		fprintf(file, "atot (2): %e, %e, %e\n", a_tot[0][2], a_tot[1][2], a_tot[2][2]);
		fprintf(file, "atot (3): %e, %e, %e\n", a_tot[0][3], a_tot[1][3], a_tot[2][3]);
#endif
		fprintf(file, "areg (0): %e, %e, %e\n", a_reg[0][0], a_reg[1][0], a_reg[2][0]);
#ifdef PRINT_FULL_ACC
		fprintf(file, "areg (1): %e, %e, %e\n", a_reg[0][1], a_reg[1][1], a_reg[2][1]);
		fprintf(file, "areg (2): %e, %e, %e\n", a_reg[0][2], a_reg[1][2], a_reg[2][2]);
		fprintf(file, "areg (3): %e, %e, %e\n", a_reg[0][3], a_reg[1][3], a_reg[2][3]);
#endif
		fprintf(file, "airr (0): %e, %e, %e\n", a_irr[0][0], a_irr[1][0], a_irr[2][0]);
#ifdef PRINT_FULL_ACC
		fprintf(file, "airr (1): %e, %e, %e\n", a_irr[0][1], a_irr[1][1], a_irr[2][1]);
		fprintf(file, "airr (2): %e, %e, %e\n", a_irr[0][2], a_irr[1][2], a_irr[2][2]);
		fprintf(file, "airr (3): %e, %e, %e\n", a_irr[0][3], a_irr[1][3], a_irr[2][3]);
#endif
		fprintf(file, "TimeStepIrr (Myr): %e, TimeStepReg (Myr): %e\n", TimeStepIrr*EnzoTimeStep*1e4, TimeStepReg*EnzoTimeStep*1e4);
	}
};




#endif
