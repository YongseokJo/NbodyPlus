#include "Particle.h"
#include <vector>
#include <iostream>
#include <cmath>


// remeber that ACList would be a vector list

void Particle::calculateIrrForce(double next_time) {

    // temporary variables for calculation

    double dt;
    double dtsq,dt6,dt2,dtsq12,dt13;

    double dx[3];
    double dv[3];
    double firr[3],firrdot[3],fdum[3];

    double rij2;
    double dr2i,dr3i,drdv,drdp;

    double df,fid,sum,at3,bt2;

    // predict poistion and velocity of ith particle
    dt = next_time - CurrentTimeIrr;
    predictPosAndVel(next_time);

    // initialize irregular force terms for ith particle
    for (int dim=0; dim<Dim; dim++){
        FIrr[dim] = 0.0;
    }

    // add the contribution of jth particle to irregular force
    for (Particle* element:ACList) {

        rij2 = 0;
        element->predictPosAndVel(next_time);
        
        for (int dim=0; dim<Dim; dim++) {

            dx[dim] = element->PredPosition[dim] - PredPosition[dim];
            dv[dim] = element->PredVelocity[dim] - PredVelocity[dim];

            rij2 += dx[dim]*dx[dim];
            drdv += dx[dim]*dv[dim];
        }

        dr2i = 1.0/rij2;
        dr3i = element->Mass*dr2i*sqrt(dr2i);
        drdp = 3.0*drdv*dr2i;

        for (int dim=0; dim<Dim; dim++) {
            firr[dim] += dx[dim]*dr3i;  // to SY is it correct? there was 3 in the indicies.
            firrdot[dim] += (dv[dim] - dx[dim]*drdp)*dr3i;
        }
    }

    // correct the force using the 4th order hermite method

    dtsq = dt*dt;
    dt6 = 6.0/(dt*dtsq); // to SY is it dtsq? it was dts1.
    dt2 = 2.0/dtsq;
    dtsq12 = dtsq/12.0;
    dt13 = dt/3.0;

    for (int dim=0; dim<Dim; dim++) {

        df  = FIrr[dim] - firr[dim];
        fid = FIrrDot[dim];
        sum = FIrrDot[dim] + firrdot[dim];
        at3 = 2.0*df + dt*sum;
        bt2 = -3.0*df - dt*(sum + fid);

				//to SY why do you need it? you don't use it.
        //xn[dim]    = PredPosition[dim] + (0.6*at3 + bt2)*dtsq12;
        //xndot[dim] = PredVelocity[dim] + (0.75*at3 + bt2)*dt13; 

        FIrr[dim] = firr[dim];
        FIrrDot[dim] = firrdot[dim];

        fdum[dim] = FIrr[dim] + FReg[dim];

        dFIrr[dim][2] = (3.0*at3 + bt2)*dt2;
        dFIrr[dim][3] = at3*dt6;
    }

}


void Particle::predictPosAndVel(double next_time){

    // temporary variables for calculation

    double dt,dt1,dt2;

    dt = next_time - CurrentTimeIrr;
    dt1 = 1.5*dt;
    dt2 = 2.0*dt;

    for (int dim=0; dim<Dim; dim++) {
        PredPosition[dim] = ((ForceDot[dim]*dt + Force[dim])*dt + Velocity[dim])*dt + Position[dim];
        PredVelocity[dim] = (ForceDot[dim]*dt1 + Force[dim])*dt2 + Velocity[dim];
    }

}
