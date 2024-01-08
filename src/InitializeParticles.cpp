#include <vector>
#include <iostream>
#include "Particle/Particle.h"
#include <cmath>




void findNeighbor(std::vector<Particle*> &particle);
void calculateForce1(std::vector<Particle*> &particle);
void calculateForce2(std::vector<Particle*> &particle);



void initializeParticle(std::vector<Particle*> &particle) {

  	std::cout << "Initialization starts." << std::endl;
    findNeighbor(particle);

    calculateForce1(particle);
    calculateForce2(particle);
	  std::cout << "Timestep initializing..." << std::endl;
		for (Particle* elem:particle)
			elem->initializeTimeStep();

  	std::cout << "Initialization finished." << std::endl;
}

void findNeighbor(std::vector<Particle*> &particle) {

    double rs0 = 0.1 ; // I need to change this to 4 parsecs
    double rij;


    int ACnum;
    int j1;


  	std::cout << "finding neightbors" << std::endl;
    // initialize the guess for radius of neighbor sphere


    for (int i1=0; i1 < NNB; i1++) {
        particle[i1]->RadiusOfAC = rs0;
    }

    // search for neighbors for particle i2

    for (int i2=0; i2 < NNB; i2++) {
        ACnum = 0;
        j1 = 0;

    // increase neighbor radius until we have at least one neighbor

        while (j1 < NNB) {

            rij = particle[i2]->getDistanceTo(particle[j1]);

            if (rij < particle[i2]->RadiusOfAC && i2 != j1) {
                ACnum++;
                particle[i2]->ACList.push_back(particle[j1]);
            }
						j1++;
            if ( ACnum == 0 && j1 >= NNB ){
                j1 = 0;
                particle[i2]->RadiusOfAC *= 1.59;
            }
        }
        particle[i2]->NumberOfAC = ACnum;
    }
		/***
		if (debug) {
			std::cout << particle[2]->NumberOfAC << std::endl;
			for (Particle* element : particle[2]->ACList)
				std::cout << element->getPID() << " ";
		}
		for (int i=0; i<NNB; i++) {
			std::cout << particle[i]->NumberOfAC << ' ';
		}
		***/
		std::cout << "\n";
  	std::cout << "finding neightbors finished." << std::endl;
}


void calculateForce1(std::vector<Particle*> &particle) {

    // temporary variables for calculation
    double A[6];
    double A7,A8,A9;
    double F1[3];
    double F1DOT[3];

    int jAC;

    // softening length, hard coded
    double EPS2 = 0.0001; //need to adjust this parameter before running


  	std::cout << "entering calculateForce1" << std::endl;
    for (int i=0; i < NNB; i++) {

        jAC = 0;

        // calculate the irregular and regular force components

        for (int j = 0; j<NNB; j++) {

            if (i != j){

                for (int k1=0; k1<3 ; k1++) {
                    A[k1] = particle[j]->Position[k1] - particle[i]->Position[k1];
                    A[k1+3] = particle[j]->Velocity[k1] - particle[i]->Velocity[k1];
                }

                A7 = 1.0/(A[0]*A[0] + A[1]*A[1] + A[2]*A[2] + EPS2);
                A8 = (particle[j]->Mass)*A7*sqrt(A7);
                A9 = 3.0*(A[0]*A[3] + A[1]*A[4] + A[2]*A[5])*A7;

                for (int k2=0; k2<3 ; k2++){
                    F1[k2] = A[k2]*A8;
                    F1DOT[k2] = (A[k2+3]-A[k2]*A9)*A8;
                }
            }

            if (particle[j] != particle[i]->ACList[jAC]) {

                for (int k3=0; k3<3 ; k3++){
                    particle[i]->FReg[k3] += F1[k3];
                    particle[i]->dFReg[k3][1] += F1DOT[k3];
                }

            } else {

                for (int k4=0; k4<3 ; k4++){
                    particle[i]->FIrr[k4] += F1[k4];
                    particle[i]->dFIrr[k4][1] += F1DOT[k4];
                }

                if (jAC < (particle[i]->NumberOfAC-1)){
                    jAC++;
                }
            }
        }

        for (int dim=0; dim<3; dim++){
            particle[i]->Force[dim] = particle[i]->FIrr[dim] + particle[i]->FReg[dim];
            particle[i]->ForceDot[dim] = particle[i]->dFIrr[dim][1] + particle[i]->dFReg[dim][1];
            particle[i]->dFIrr[dim][0] = particle[i]->FIrr[dim];
            particle[i]->dFReg[dim][0] = particle[i]->FReg[dim];
            particle[i]->FIrrDot[dim] = particle[i]->dFIrr[dim][1];
            particle[i]->FRegDot[dim] = particle[i]->dFReg[dim][1];
        }

    }
}



void calculateForce2(std::vector<Particle*> &particle) {

    // temporary variables for calculation
    double A[12];
    double A13,A14,A15,A16,A17,A18,A19,A20,A21,A22;
    double F1DOTK;
    double F2DOT[3];
    double F3DOT[3];

    double Rij2;
    double DT,DT1,DTR,DT1R;

    int jAC;

    // softening length, hard coded
    double EPS2 = 4.0; //need to adjust this parameter before running


  	std::cout << "entering calculateForce2" << std::endl;

    for (int i=0; i<NNB; i++){

        jAC = 0;

        for (int j=0; j<NNB; j++){

            if (i != j){

                for (int k1=0; k1<3 ; k1++){
                    A[k1] = particle[j]->Position[k1] - particle[i]->Position[k1];
                }

                Rij2 = A[0]*A[0] + A[1]*A[1] + A[2]*A[2];

                if ((Rij2>(3*particle[i]->RadiusOfAC)) && (particle[j] != particle[i]->ACList[jAC])) {

                 // empty block

                } else {

                    if (global_time <= 0.0){

                        for (int k2=0; k2<3 ; k2++){
                            A[k2+6] =  particle[j]->Force[k2];
                            A[k2+9] =  particle[j]->ForceDot[k2];
                        }
                    } else {

                        DT = global_time - particle[j]->CurrentTimeIrr;
                        DT1 = 0.5*DT;

                        DTR = global_time - particle[j]->CurrentTimeReg;
                        DT1R = 0.5*DTR;

                        for (int k3=0; k3<3 ; k3++){
                            A[k3+6] =  (particle[j]->dFIrr[k3][2]*DT1 + particle[j]->dFIrr[k3][1])*DT + particle[j]->FIrr[k3] \
                                        + (particle[j]->dFReg[k3][2]*DT1R + particle[j]->dFReg[k3][1])*DTR + particle[j]->FReg[k3];

                            A[k3+9] = particle[j]->dFIrr[k3][2]*DT + particle[j]->dFIrr[k3][1] \
                                      + particle[j]->dFReg[k3][2]*DTR + particle[j]->dFReg[k3][1];
                        }
                    }

                    for (int k4=0; k4<3 ; k4++){
                        A[k4+3] =  particle[j]->Velocity[k4] - particle[i]->Velocity[k4];
                        A[k4+6] -=  particle[i]->Force[k4];
                        A[k4+9] -=  particle[i]->ForceDot[k4];
                    }

                    A13 = 1.0/(Rij2 + EPS2);
                    A14 = particle[j]->Mass*A13*sqrt(A13);
                    A15 = (A[0]*A[3] + A[1]*A[4] + A[2]*A[5])*A13;
                    A16 = A15*A15;
                    A17 = 3.0*A15;
                    A18 = 6.0*A15;
                    A19 = 9.0*A15;
                    A20 = (A[3]*A[3] + A[4]*A[4] + A[5]*A[5] + A[0]*A[6] + \
                           A[1]*A[7] + A[2]*A[8] )*A13 + A16;
                    A21 = 9.0*A20;
                    A20 = 3.0*A20;
                    A22 = (9.0*(A[3]*A[6] + A[4]*A[7] + A[5]*A[8] ) + \
                           3.0*(A[0]*A[9] + A[1]*A[10] + A[2]*A[11]))*A13 + \
                                                    A17*(A20 - 4.0*A16);

                    for (int k5=0; k5<3 ; k5++){
                        F1DOTK = A[k5+3] - A17*A[k5];
                        F2DOT[k5] = (A[k5+6] - A18*F1DOTK - A20*A[k5])*A14;
                        F3DOT[k5] = (A[k5+9] - A21*F1DOTK - A22*A[k5])*A14 - A19*F2DOT[k5];
                    }

                    if (particle[j] != particle[i]->ACList[jAC]){

                        for (int k6=0; k6<3 ; k6++){
                            particle[i]->dFReg[k6][2] += F2DOT[k6];
                            particle[i]->dFReg[k6][3] += F3DOT[k6];
                        }

                    } else {

                        for (int k7=0; k7<3 ; k7++){
                            particle[i]->dFIrr[k7][2] += F2DOT[k7];
                            particle[i]->dFIrr[k7][3] += F3DOT[k7];
                        }

                        if (jAC < (particle[i]->NumberOfAC-1)){
                            jAC++;
                        }
                    }
                }
            }
        }
    }
}



