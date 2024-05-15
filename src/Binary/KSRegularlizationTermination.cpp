#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"


void generate_Matrix(double a[3], double (&A)[3][4]);
void ReInitializeKSParticle(Particle* KSParticle, std::vector<Particle*> &particle);

void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle, double current_time){

    double R[Dim], Rdot[Dim];
    double Rinv;
    double ratioM;

    double L[3][4];
    int ptclCMIndex;
    int ptclBinIndex;

    bool findPtclCM;

    Particle* ptclI;
    Particle* ptclJ;

    Binary* ptclBin;

    fprintf(stdout,"--------------------------------------\n");
    fprintf(stdout,"In KSRegularlizationTermination.cpp...\n\n");

    ptclI = ptclCM->BinaryParticleI;
    ptclJ = ptclCM->BinaryParticleJ;
    ptclBin = ptclCM->BinaryInfo;

    fprintf(stdout,"Converting the KS coordinates to physical coordinates of ptclI and ptclJ\n");

    // update the values of positions of ptclI and ptcl J

    R[0]   = ptclBin->u[0]*ptclBin->u[0] - ptclBin->u[1]*ptclBin->u[1] - ptclBin->u[2]*ptclBin->u[2] + ptclBin->u[3]*ptclBin->u[3];
    R[1]   = 2*(ptclBin->u[0]*ptclBin->u[1] - ptclBin->u[2]*ptclBin->u[3]);
    R[2]   = 2*(ptclBin->u[0]*ptclBin->u[2] + ptclBin->u[1]*ptclBin->u[3]);
    ratioM = ptclJ->Mass/ptclCM->Mass;


    for (int dim=0; dim<Dim; dim++) {
        ptclI->Position[dim] = ptclCM->Position[dim] + ratioM*R[dim];
        ptclJ->Position[dim] = ptclCM->Position[dim] - R[dim];
    }


    // do the same thing for velocity components


    generate_Matrix(ptclBin->u,L);

    Rinv = 1/(ptclBin->u[0]*ptclBin->u[0] + ptclBin->u[1]*ptclBin->u[1] + ptclBin->u[2]*ptclBin->u[2] + ptclBin->u[3]*ptclBin->u[3]) ;


    for (int dim=0; dim<Dim; dim++) {

        Rdot[dim] = 0.0;

        for (int dimu=0; dimu<4; dimu++) {
            Rdot[dim] += 2*L[dim][dimu]*ptclBin->udot[dim]*Rinv;
        }
    }


    for (int dim=0; dim<Dim; dim++) {
        ptclI->Velocity[dim] = ptclCM->Velocity[dim] + ratioM*Rdot[dim];
        ptclJ->Velocity[dim] = ptclI->Velocity[dim] - Rdot[dim];
    }


    fprintf(stdout,"END CONVERTING THE COORDINATES\n \n");

    // delete the original components from the list

    fprintf(stdout,"deleting CM particle from the particle list\n");

    ptclCMIndex = -1;

    for (Particle* ptcl : particle) {

        ptclCMIndex += 1;

        if (ptcl == ptclCM) {
            break;
        }
    }

    particle.erase(particle.begin() + ptclCMIndex);


    // add the original particles

    fprintf(stdout,"add the binary components to particle list\n");

    particle.push_back(ptclI);
    particle.push_back(ptclJ);

    ptclI->CurrentTimeIrr = current_time; //ptclCM->CurrentTimeIrr;
    ptclI->CurrentTimeReg = current_time; //ptclCM->CurrentTimeReg;

    ptclJ->CurrentTimeIrr = current_time; //ptclCM->CurrentTimeIrr;
    ptclJ->CurrentTimeReg = current_time; //ptclCM->CurrentTimeReg;


    fprintf(stdout,"initialize particle I \n");
    ReInitializeKSParticle(ptclI, particle);
    fprintf(stdout,"initialize particle J \n");
    ReInitializeKSParticle(ptclJ, particle);


    // we also need to revert the neighbor list of Particles
    // assuming that all the neighbors are bidirectional
    // may need to update later if the radius for neighbor differs depending on the particle

    fprintf(stdout,"replacing CM particle in neighbor list to component particles \n");

    for (Particle* ptcl1: particle) {

        //auto it = std::find(ptcl->ACList.begin(), ptcl->ACList.end(), ptclCM);
        
        //if (it != ptclJ->ACList.end()) {
        //    ptcl->ACList.erase(it);
        //    ptcl->ACList.push_back(ptclI);
        //    ptcl->ACList.push_back(ptclJ);
        //}
	
	ptclCMIndex = -1;
	findPtclCM = false;

	for (Particle* ptcl2 : ptcl1->ACList) {

       	    ptclCMIndex += 1;

            if (ptcl2 == ptclCM) {
		findPtclCM = true;
                break;
            }
        }

	if (findPtclCM) {
	    particle.erase(particle.begin() + ptclCMIndex);
	    ptcl1->ACList.push_back(ptclI);
	    ptcl1->ACList.push_back(ptclJ);
	}

    }

    // we also need to delete it from the binary list

    fprintf(stdout,"deleting binary information from the BinaryList \n");

    ptclBinIndex = -1;

    for (Binary* ptcl : BinaryList) {

        ptclBinIndex += 1;

        if (ptcl == ptclBin) {
            break;
        }
    }

    BinaryList.erase(BinaryList.begin() + ptclBinIndex);

    fprintf(stdout,"end of KS Regularlization Termination \n ");

    fprintf(binout, "\nPosition: ptclI - x:%e, y:%e, z:%e, \n", ptclI->Position[0], ptclI->Position[1], ptclI->Position[2]);
    fprintf(binout, "Velocity: ptclI - vx:%e, vy:%e, vz:%e, \n", ptclI->Velocity[0], ptclI->Velocity[1], ptclI->Velocity[2]);
    fflush(binout);

    fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclI->a_tot[0][0], ptclI->a_tot[1][0], ptclI->a_tot[2][0]);
    fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclI->a_tot[0][1], ptclI->a_tot[1][1], ptclI->a_tot[2][1]);
    fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclI->a_tot[0][2], ptclI->a_tot[1][2], ptclI->a_tot[2][2]);
    fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclI->a_tot[0][3], ptclI->a_tot[1][3], ptclI->a_tot[2][3]);
    fprintf(binout, "Time Steps - irregular:%e, regular:%e \n", ptclI->TimeStepIrr, ptclI->TimeStepReg);
    fflush(binout);


    fprintf(binout, "\nPosition: ptclJ - x:%e, y:%e, z:%e, \n", ptclJ->Position[0], ptclJ->Position[1], ptclJ->Position[2]);
    fprintf(binout, "Velocity: ptclJ - vx:%e, vy:%e, vz:%e, \n", ptclJ->Velocity[0], ptclJ->Velocity[1], ptclJ->Velocity[2]);
    fflush(binout);

    fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclI->a_tot[0][0], ptclI->a_tot[1][0], ptclJ->a_tot[2][0]);
    fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclJ->a_tot[0][1], ptclI->a_tot[1][1], ptclI->a_tot[2][1]);
    fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclJ->a_tot[0][2], ptclJ->a_tot[1][2], ptclJ->a_tot[2][2]);
    fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclJ->a_tot[0][3], ptclJ->a_tot[1][3], ptclJ->a_tot[2][3]);
    fprintf(binout, "Time Steps - irregular:%e, regular:%e \n", ptclJ->TimeStepIrr, ptclJ->TimeStepReg);
    fflush(binout);


    fprintf(binout, "\nPosition: ptclCM - x:%e, y:%e, z:%e, \n", ptclCM->Position[0], ptclCM->Position[1], ptclCM->Position[2]);
    fprintf(binout, "Velocity: ptclCM - vx:%e, vy:%e, vz:%e, \n", ptclCM->Velocity[0], ptclCM->Velocity[1], ptclCM->Velocity[2]);
    fflush(binout);

    fprintf(binout, "Total Acceleration - ax:%e, ay:%e, az:%e, \n", ptclCM->a_tot[0][0], ptclCM->a_tot[1][0], ptclCM->a_tot[2][0]);
    fprintf(binout, "Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", ptclCM->a_tot[0][1], ptclCM->a_tot[1][1], ptclCM->a_tot[2][1]);
    fprintf(binout, "Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", ptclCM->a_tot[0][2], ptclCM->a_tot[1][2], ptclCM->a_tot[2][2]);
    fprintf(binout, "Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", ptclCM->a_tot[0][3], ptclCM->a_tot[1][3], ptclCM->a_tot[2][3]);
    fprintf(binout, "Time Steps - irregular:%e, regular:%e \n", ptclCM->TimeStepIrr, ptclCM->TimeStepReg);
    fflush(binout);

    fprintf(stdout,"--------------------------------------\n");
    fflush(stdout); 

}
