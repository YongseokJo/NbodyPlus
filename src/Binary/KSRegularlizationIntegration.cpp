#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "../global.h"
#include "../defs.h"


void generate_Matrix(double a[3], double (&A)[3][4]);
void direct_sum(double *x, double *v, double r2, double vx,
		        double mass, double (&a)[3], double (&adot)[3]);


// refer to ksint.f

void Binary::KSIntegration(Particle* ptclCM){

    double next_time;


    next_time = CurrentTime + TimeStep;

    // first predict the future position of the binary particle. 

    predictBinary(ptclCM, next_time);

    // if there are zero neighbors for binary particles, calculate by unperturbed equations

    IntegrateBinary(ptclCM, next_time);

}



void Binary::predictBinary(Particle* ptclCM, double next_time) {

    double dt, dt2, dt3;
    double dtau, dtau2, dtau3, dtau4;
    double rinv, rinv2, rinv3, rinv4, rinv5;

    double R[Dim], Rdot[Dim];
    double Rinv;
    double ratioM;

    double L[3][4];

    Particle* ptclI;
    Particle* ptclJ;

    ptclI = ptclCM->BinaryParticleI;
    ptclJ = ptclCM->BinaryParticleJ;


    // time interval for prediction

    dt = next_time - CurrentTime;
    dt2 = dt*dt;
    dt3 = dt2*dt;

    // inverse of r

    rinv = 1/r;
    rinv2 = rinv/r;
    rinv3 = rinv2/r;
    rinv4 = rinv3/r;
    rinv5 = rinv4/r;


    // the binary particle prediction consists of two parts
    // first, we need to predict the future position of the center of mass

    for (int dim=0; dim<Dim; dim++) {
		PredPosition[dim] = Position[dim] + Velocity[dim]*dt + a_tot[dim][0]*dt2/2 + a_tot[dim][1]*dt3/6;
		PredVelocity[dim] = Velocity[dim] + a_tot[dim][0]*dt + a_tot[dim][1]*dt2/2;
	}


    // then we need to convey this information to the single particles as well
    // originally, the prediction order differs depending on conditions
    // but for now I would only use the 3rd order prediction

    // dtau = dtau0 - 

    dtau = dt*rinv - t2dot*dt2*rinv3/2 + t2dot*t2dot*dt3*rinv5/2 - t3dot*dt3*rinv4/6;

    // if (abs(dtau)>dTau) {
    //     dtau = 0.8*dTau;
    // }

    dtau2 = dtau*dtau;
    dtau3 = dtau2*dtau;
    dtau4 = dtau3*dtau;


    for (int dim=0; dim<4; dim++) {
        u_pred[dim]    = u[dim]    + udot[dim]*dtau  + u2dot[dim]*dtau2/2 + u3dot[dim]*dtau3/6 + u4dot[dim]*dtau4/24;
        udot_pred[dim] = udot[dim] + u2dot[dim]*dtau + u3dot[dim]*dtau2/2 + u4dot[dim]*dtau3/6;
    }

    
    // convert things back to global coordinates
    // r is the distance vector between the two binaries, and R is the distance between them
    // using the predicted R and predicted cm values, we obtain the position of each binary component at a future time


    R[0]   = u_pred[0]*u_pred[0] - u_pred[1]*u_pred[1] - u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3];
    R[1]   = 2*(u_pred[0]*u_pred[1] - u_pred[2]*u_pred[3]);
    R[2]   = 2*(u_pred[0]*u_pred[2] + u_pred[1]*u_pred[3]);
    ratioM = ptclJ->Mass/Mass;


    for (int dim=0; dim<Dim; dim++) {
        ptclI->PredPosition[dim] = PredPosition[dim] + ratioM*R[dim];
        ptclJ->PredPosition[dim] = ptclI->PredPosition[dim] - R[dim];
    }


    // do the same thing for velocity components


    generate_Matrix(u_pred,L);

    Rinv = 1/(u_pred[0]*u_pred[0] + u_pred[1]*u_pred[1] + u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3]) ;


    for (int dim=0; dim<Dim; dim++) {

        Rdot[dim] = 0.0;

        for (int dimu=0; dimu<4; dimu++) {
            Rdot[dim] += 2*L[dim][dimu]*udot_pred[dim]*Rinv;
        }
    }


    for (int dim=0; dim<Dim; dim++) {
        ptclI->PredVelocity[dim] = PredVelocity[dim] + ratioM*Rdot[dim];
        ptclJ->PredVelocity[dim] = ptclI->PredVelocity[dim] - Rdot[dim];
    }

}








void Binary::IntegrateBinary(Particle* ptclCM, double next_time) {

    // variables for position, velocity prediction

    double dt, dt2, dt3;
    double rinv, rinv2, rinv3, rinv4, rinv5;

    double R[Dim], Rdot[Dim];
    double Rinv;
    double ratioM;

    double L[3][4];
    double dtau, dtau2, dtau3, dtau4, dtau5, dtau6;


    // variables for perturbing acceleration calculation

    double P[Dim], Pdot[Dim];  // perturbing acceleration and jerk
    double dxi[Dim],dvi[Dim],dxj[Dim],dvj[Dim];
    double dr2i, dr2j, dxdvi, dxdvj;


    // variables for calculating correction factors

    double dh0, dh;
    double h3dot_hermite, h4dot_hermite; 

    double Q_pred[4], Qdot_pred[4], Q2dot_pred[4];
    double L[3][4], Ldot[3][4], L2dot[3][4]; 

    double r_pred, rdot_pred, hdot_pred, h2dot_pred;
    double u2dot_pred[4], u3dot_pred[4];
    double u4dot_hermite[4], u5dot_hermite[4];

    double z;


    Particle* ptclI;
    Particle* ptclJ;

    ptclI = ptclCM->BinaryParticleI;
    ptclJ = ptclCM->BinaryParticleJ;


    dtau = dTau;
    dtau2 = dtau*dtau;
    dtau3 = dtau2*dtau;
    dtau4 = dtau3*dtau;
    dtau5 = dtau4*dtau;
    dtau6 = dtau5*dtau;


    // refer to kspred.f

    // first predict the positions of neighbors

    // initialize the perturbing acceleration and jerk

    for (int dim=0; dim<Dim; dim++) {
        P[dim] = 0.0;
        Pdot[dim] = 0.0;
    }


    if (ptclCM->NumberOfAC >1) {

        for (Particle* ptcl: ptclCM->ACList) {
            ptcl->predictParticleSecondOrder(next_time);
        }

        // predict the positions of binary pair particles to highest order possible
        // using stumpff coefficients


        for (int dimu=0; dimu<4; dimu++) {
            u_pred[dimu]    = u[dimu] + udot[dimu]*dtau + u2dot[dimu]*dtau2/2 + u3dot[dimu]*dtau3/6 \
                    + cn[4]*u4dot[dimu]*dtau4/24 + cn[5]*u5dot[dimu]*dtau5/120;
            udot_pred[dimu] = udot[dimu] + u2dot[dimu]*dtau + u3dot[dimu]*dtau2/2 \
                    + cn[4]*u4dot[dimu]*dtau3/6+ cn[5]*u5dot[dimu]*dtau4/24;
        }


        R[0]   = u_pred[0]*u_pred[0] - u_pred[1]*u_pred[1] - u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3];
        R[1]   = 2*(u_pred[0]*u_pred[1] - u_pred[2]*u_pred[3]);
        R[2]   = 2*(u_pred[0]*u_pred[2] + u_pred[1]*u_pred[3]);
        ratioM = ptclJ->Mass/Mass;


        for (int dim=0; dim<Dim; dim++) {
            ptclI->PredPosition[dim] = PredPosition[dim] + ratioM*R[dim];
            ptclJ->PredPosition[dim] = ptclI->PredPosition[dim] - R[dim];
        }

        generate_Matrix(u_pred,L);

        r_pred = (u_pred[0]*u_pred[0] + u_pred[1]*u_pred[1] + u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3]) ;


        for (int dim=0; dim<Dim; dim++) {

            Rdot[dim] = 0.0;

            for (int dimu=0; dimu<4; dimu++) {
                Rdot[dim] += 2*L[dim][dimu]*udot_pred[dim]/r_pred;
            }
        }

        for (int dim=0; dim<Dim; dim++) {
            ptclI->PredVelocity[dim] = PredVelocity[dim] + ratioM*Rdot[dim];
            ptclJ->PredVelocity[dim] = ptclI->PredVelocity[dim] - Rdot[dim];
        }


        // calculate perturbation from single particles

        for (Particle* ptcl: ACList) {

            dr2i = 0;
            dr2j = 0;

            dxdvi = 0;
            dxdvj = 0;

            for (int dim=0; dim<Dim; dim++) {

                dxi[dim] = ptcl->PredPosition[dim] - ptclI->PredPosition[dim];
                dvi[dim] = ptcl->PredVelocity[dim] - ptclJ->PredVelocity[dim];
                dr2i += dxi[dim]*dxi[dim];
                dxdvi += dxi[dim]*dvi[dim];

                dxj[dim] = ptcl->PredPosition[dim] - ptclI->PredPosition[dim];
                dvj[dim] = ptcl->PredVelocity[dim] - ptclJ->PredVelocity[dim];
                dr2j += dxj[dim]*dxj[dim];
                dxdvj += dxj[dim]*dvj[dim];

            }

            direct_sum(dxi ,dvi, dr2i, dxdvi, ptcl->Mass, P, Pdot);
            direct_sum(dxj ,dvj, dr2j, dxdvj, -ptcl->Mass, P, Pdot);

        }

        
        for (int dim=0; dim<Dim; dim++) {
            Pdot[dim] *= r_pred;
        }

    }

    // calculating perturbation from particle systems
    // add later on

    gamma = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2])*r_pred*r_pred/Mass;



    // refer to kscorr.f

    dh = hdot*dtau + h2dot*dtau2/2 + h3dot*dtau3/6 + h4dot*dtau4/24;
    dh0 = hdot*dtau + h2dot*dtau2/2;


    // caculation of lower order derivatives

    // initialize the variables

    for(int dimu=0; dimu<4; dimu++) {
        Q_pred[dimu] = 0.0;
        Qdot_pred[dimu] = 0.0;
        u2dot_pred[dimu] = 0.0;
        u3dot_pred[dimu] = 0.0;
        u4dot_hermite[dimu] = 0.0;
        u5dot_hermite[dimu] = 0.0;
    }

    rdot_pred = 0.0; 
    hdot_pred = 0.0;

    generate_Matrix(u_pred,L);
    generate_Matrix(udot_pred,Ldot);

     
    
    for (int dimu=0; dimu<4; dimu++) {

        Q_pred[dimu]     = L[0][dimu]*P[0] + L[1][dimu]*P[1] + L[2][dimu]*P[2];
        Qdot_pred[dimu]  = Ldot[0][dimu]*P[0] + Ldot[1][dimu]*P[1] + Ldot[2][dimu]*P[2] \
                         + L[0][dimu]*Pdot[0] + L[1][dimu]*Pdot[1] + L[2][dimu]*Pdot[2];

        u2dot_pred[dimu] = 0.5*(h+dh)*u_pred[dimu] + 0.5*r_pred*Q_pred[dimu];

        rdot_pred       += 2*u_pred[dimu]*udot_pred[dimu];
        hdot_pred       += 2*udot_pred[dimu]*Q_pred[dimu];
        h2dot_pred      += 2*(u2dot_pred[dimu]*Q_pred[dimu] + udot_pred[dimu]*Qdot_pred[dimu]);
        
    }

    for (int dimu = 0; dimu<4; dimu++) {
        u3dot_pred[dimu] = 0.5*((h+dh)*udot_pred[dimu] + hdot_pred*u_pred[dimu]) \
                         + 0.5*(rdot_pred*Q_pred[dimu] + r_pred*Qdot_pred[dimu]);
    }


    // apply the hermite corrector method to calculate 4th and 5th derivatives of u
    // and also the 3rd and 4th derivatives of h

    for (int dimu = 0; dimu<4; dimu++) {
        u4dot_hermite[dimu] = - 6*(u2dot[dimu] - u2dot_pred[dimu])/dtau2 \
                              - 2*(2*u3dot[dimu] + u3dot_pred[dimu])/dtau;
        u5dot_hermite[dimu] =   12*(u2dot[dimu] - u2dot_pred[dimu])/dtau3 \
                              + 6*(u3dot[dimu] + u3dot_pred[dimu])/dtau2;
    }

    h3dot_hermite = -6*(hdot-hdot_pred)/dtau2 -2*(2*h2dot + h2dot_pred)/dtau;
    h4dot_hermite = 12*(hdot-hdot_pred)/dtau3 + 6*(h2dot + h2dot_pred)/dtau2;


    // based on the higher-order derivatives, correct the values of u, udot and r

    for (int dimu=0; dimu<4; dimu++) {
        u_pred[dimu]    = u[dimu] + udot[dimu]*dtau + u2dot[dimu]*dtau2/2 + u3dot[dimu]*dtau3/6 \
                        + cn[4]*u4dot_hermite[dimu]*dtau4/24 + cn[5]*u5dot_hermite[dimu]*dtau5/120;
        udot_pred[dimu] = udot[dimu] + u2dot[dimu]*dtau + u3dot[dimu]*dtau2/2 \
                        + cn[4]*u4dot_hermite[dimu]*dtau3/6+ cn[5]*u5dot_hermite[dimu]*dtau4/24;           
    }

    r_pred = (u_pred[0]*u_pred[0] + u_pred[1]*u_pred[1] + u_pred[2]*u_pred[2] + u_pred[3]*u_pred[3]);


    // intialize and save the updated values of all relevant derivatives

    r = r_pred;
    rdot = 0.0;
    r2dot = 0.0;
    r3dot = 0.0;
    r4dot = 0.0;
    r5dot = 0.0;

    h = h + dh0 + h3dot_hermite*dtau3/6 + h4dot_hermite*dtau4/24;
    hdot = hdot_pred;
    h2dot = h2dot_pred;
    h3dot = h3dot_hermite + h4dot_hermite*dtau;
    h4dot = h4dot_hermite;

    for (int dimu=0; dimu<4 ;dimu++) {

        Q[dimu] = Q_pred[dimu];
        Qdot[dimu] = Qdot_pred[dimu];

        u[dimu] = u_pred[dimu];
        udot[dimu] = udot_pred[dimu];
        u2dot[dimu] = 0.5*h*u_pred[dimu] + 0.5*r_pred*Q_pred[dimu];
        u3dot[dimu] = 0.5*(h*udot_pred[dimu] + hdot_pred*u_pred[dimu]) \
                    + 0.5*(rdot_pred*Q_pred[dimu] + r_pred*Qdot_pred[dimu]);
        u4dot[dimu] = u4dot_hermite[dimu] + u5dot_hermite[dimu]*dTau;
        u5dot[dimu] = u5dot_hermite[dimu];
    }


    for (int dimu=0; dimu<4; dimu++) {
        rdot      += 2*u[dimu]*udot[dimu];
        r2dot     += 2*(udot[dimu]*udot[dimu] + u[dimu]*u2dot[dimu]);
        r3dot     += 2*(3*u2dot[dimu]*udot[dimu] + u3dot[dimu]*u[dimu]);
        r4dot     += 2*(u4dot[dimu]*u[dimu] + 4*u3dot[dimu]*udot[dimu] + 3*u2dot[dimu]*u2dot[dimu]);
        r5dot     += 2*(u5dot[dimu]*u[dimu] + 5*u4dot[dimu]*udot[dimu] + 10*u3dot[dimu]*u2dot[dimu]);
    }

    tdot = r;
    t2dot = rdot;
    t3dot = r2dot;
    t4dot = r3dot;
    t5dot = r4dot;
    t6dot = r5dot;


    // update the time

    CurrentTime = next_time;


    // generate new stumpff coefficients
    // refer to stumpff.for

    z = -0.5*h*dtau;
    getStumpffCoefficients(z);

    // TimeStep = tdot*dtau + t2dot*dtau2/2 + t3dot*dtau3/6 + t4dot*dtau4/24 \
    //          + t5dot*cn_4z[5]*dtau5/120 + t6dot*cn_4z[6]*dtau6/720;

    // generate new list of perturbers

    CurrentTime = next_time;
    CurrentTau += dtau;
    
}