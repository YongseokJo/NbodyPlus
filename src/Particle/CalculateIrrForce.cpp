#include "Particle.h"
#include <vector>
#include <iostream>
#include <cmath>

// remeber that ACList would be a vector list

void Particle::calculateIrrForce() {

	if (this->NumberOfAC == 0) {
		CurrentTimeIrr += TimeStepIrr;
		return;
	}
	// temporary variables for calculation

	double dt,dtr;
	double dtsq,dt6,dt2,dtsq12,dt13;
	double irr_next_time;

	double dx[3];
	double dv[3];
	double firr[3],firrdot[3],fdum[3];

	double rij2;
	double dr2i,dr3i,drdv,drdp;

	double df,fid,sum,at3,bt2;

	// predict poistion and velocity of ith particle

	dt = TimeStepIrr; // interval of time step
	irr_next_time = CurrentTimeIrr + TimeStepIrr; // the time to be advanced to


	// only predict when regular block doesn't come
	// b/c full prediction would be performed at regular block


	predictParticleSecondOrder(irr_next_time);


	// initialize irregular force terms for ith particle
	for (int dim=0; dim<Dim; dim++){
		firr[dim] = 0.0;
		firrdot[dim] = 0.0;
	}

	// add the contribution of jth particle to irregular force
	for (Particle* element:ACList) {

		rij2 = 0;
		drdv = 0;

		// only predict when regular block doesn't come
		// b/c full prediction would be performed at regular block

		element->predictPosAndVel(irr_next_time);


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
	fprintf(stdout, "PID=%d\n", this->getPID());
	fprintf(stdout, "%lf, %lf, %lf, %lf, %lf, %lf\n",
			FIrr[0], FIrr[1], FIrr[2], FIrrDot[0], FIrrDot[1], FIrrDot[2]);
	fprintf(stdout, "%lf, %lf, %lf, %lf, %lf, %lf\n",
			firr[0], firr[1], firr[2], firrdot[0], firrdot[1], firrdot[2]);
	fprintf(stdout, "%lf, %lf, %lf, %lf, %lf, %lf\n\n",
			PredPosition[0], PredPosition[1], PredPosition[2], PredVelocity[0], PredVelocity[1], PredVelocity[2]);
	fprintf(stdout, "%lf, %lf, %lf, %lf, %lf, %lf\n\n",
			Position[0], Position[1], Position[2], Velocity[0], Velocity[1], Velocity[2]);
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

		FIrr[dim] = firr[dim];
		FIrrDot[dim] = firrdot[dim];

		fdum[dim] = FIrr[dim] + FReg[dim];

		// save the newly calculated values to variables in particle class

		Position[dim] = PredPosition[dim] + (0.6*at3 + bt2)*dtsq12;
		Velocity[dim] = PredVelocity[dim] + (0.75*at3 + bt2)*dt13;

		dFIrr[dim][0] = FIrr[dim];
		dFIrr[dim][1] = FIrrDot[dim];
		dFIrr[dim][2] = (3.0*at3 + bt2)*dt2;
		dFIrr[dim][3] = at3*dt6;

		fdum[dim] = FIrr[dim] + FReg[dim];
	}

	// save the current and next time steps
	CurrentTimeIrr = irr_next_time;

	// calculate the new time step
	this->calculateTimeStepIrr(fdum, dFIrr);

	// if regular force is not updated, update the total force
	// with the regular force that the particle holds now
	if (CurrentTimeIrr < (CurrentTimeReg + TimeStepReg)){
		dtr = CurrentTimeIrr - CurrentTimeReg;
		for (int dim=0; dim<Dim; dim++){
			Force[dim] = 0.5*(FRegDot[dim]*dtr + FReg[dim] + FIrr[dim]);
			ForceDot[dim] = (FRegDot[dim]+FIrrDot[dim]) / 6.0 ;
		}
	}

}


