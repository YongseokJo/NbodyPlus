#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"


void generate_Matrix(double a[3], double (&A)[3][4]);
void InitializeParticle(Particle* newParticle, std::vector<Particle*> &particle);

void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle){

    double R[Dim], Rdot[Dim];
    double Rinv;
    double ratioM;

    double L[3][4];
    int ptclCMIndex;
    int ptclBinIndex;

    Particle* ptclI;
    Particle* ptclJ;

    Binary* ptclBin;



    ptclI = ptclCM->BinaryParticleI;
    ptclJ = ptclCM->BinaryParticleJ;
    ptclBin = ptclCM->BinaryInfo;

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


    // delete the original components from the list

    ptclCMIndex = -1;

    for (Particle* ptcl : particle) {

        ptclCMIndex += 1;

        if (ptcl == ptclCM) {
            break;
        }
    }

    particle.erase(particle.begin() + ptclCMIndex);


    // add the original particles

    particle.push_back(ptclI);
    particle.push_back(ptclJ);

    InitializeParticle(ptclI, particle);
    InitializeParticle(ptclJ, particle);


    // we also need to revert the neighbor list of Particles
    // assuming that all the neighbors are bidirectional
    // may need to update later if the radius for neighbor differs depending on the particle

    // first for particle I

    for (Particle* ptcl: ptclCM->ACList) {

        auto it = std::find(ptcl->ACList.begin(), ptcl->ACList.end(), ptclCM);
        
        if (it != ptclJ->ACList.end()) {
            ptcl->ACList.erase(it);
            ptcl->ACList.push_back(ptclI);
            ptcl->ACList.push_back(ptclJ);
        }
    }

    // we also need to delete it from the binary list

    ptclBinIndex = -1;

    for (Binary* ptcl : BinaryList) {

        ptclBinIndex += 1;

        if (ptcl == ptclBin) {
            break;
        }
    }

    BinaryList.erase(BinaryList.begin() + ptclBinIndex);


}
