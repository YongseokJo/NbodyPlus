#include "Particle.h"
#include "../defs.h"
#include <vector>
#include <iostream>
#include <cmath>


void predictAll(std::vector<Particle*> &particle) {
	for (int i=0; i<NNB; i++){
		particle[i]->predictPosAndVel(global_time);
	}
}


void Particle::calculateRegForce(std::vector<Particle*> &particle) {

	// temporary variables for calculation, regint.for
	// variables for obtaining predicted values of new forces

	double dtr;
	double dx[3];
	double dv[3];

	double xnew[3],xnewdot[3];
	double freg[3],fregdot[3];
	double firr[3],firrdot[3];
	double dfirr[3],dfirrdot[3];
	double pot = 0.0; // potential 

	double rij2;
	double dr2i,dr3i,drdv,drdp;
	double dtsq,dt6,dt2,dtsq12,dt13;
	double dfr,fdr0,frd,sum,at3,bt2;

	double ri2;
	double rscale = 0.1;
	double rh2;
	double vrfac;


	// temporary variables for handeling neighbors

	int j = 0;
	int newACnum = 0;  // temporary variable for saving new neighbor number
	std::vector<int> newACList; // temporary list for saving new neighbor list

	bool isGain,isLoss;

	int maxNeighborNum = 250;
	double rs2 = 0.1;  // scale radius...? need to check


	// temporary variables for calculation in regcor_gpu.f

	double fij;
	double rdot;

	predictAll(particle);

	// codes from regint.f in NBODY6++GPU
	// initialize regular force terms for ith particle

	// note that all positions are already predicted, and we are calculating the predicted components

	dtr = TimeStepReg;

	for (int dim=0; dim<Dim; dim++){
		firr[dim] = 0.0;
		firrdot[dim] = 0.0;
		freg[dim] = 0.0;
		fregdot[dim] = 0.0;
	}

	// scan over all single particles to find the regular force components

	std::cout <<  "Looping single particles for force...\n" << std::flush;
	while (j<NNB) {
		if (particle[j] != this){
			j++;
			continue;
		}
		rij2 = 0.0;
		drdv = 0.0;

		for (int dim=0; dim<Dim; dim++) {

			dx[dim] = particle[j]->PredPosition[dim] - PredPosition[dim];
			dv[dim] = particle[j]->PredVelocity[dim] - PredVelocity[dim];

			rij2 += dx[dim]*dx[dim];
			drdv += dx[dim]*dv[dim];
		}

		rij2 += EPS2;  // add softening length
		dr2i = 1.0/rij2;
		dr3i = particle[j]->Mass*dr2i*sqrt(dr2i);
		drdp = 3.0*drdv*dr2i;

		if (rij2>RadiusOfAC*RadiusOfAC){

			for (int dim=0; dim<Dim; dim++){
				freg[dim] += dx[dim]*dr3i;
				fregdot[dim] += (dv[dim] - dx[dim]*drdp)*dr3i;
			}
		} else if (NNB > 100) {
			newACnum += 1;
			newACList.push_back(j) ;

			for (int dim=0; dim<Dim; dim++){
				firr[dim] += dx[dim]*dr3i;
				fregdot[dim] += (dv[dim] - dx[dim]*drdp)*dr3i;
			}
		}

		pot += dr3i*rij2;
		j++;

		// Maybe need to add new criterias present in Nbody2 version
		// which is too big neighbor radius and distant particles with large steps
		// I would include routines for zero neighbors later on...

		//if ( (newACnum>maxNeighborNum || newACnum == 0) && j>=NNB) {
		if ( (newACnum > maxNeighborNum) && j>=NNB) {

			for (int dim=0; dim<Dim; dim++){
				firr[dim] = 0.0;
				freg[dim] = 0.0;
				fregdot[dim] = 0.0;
			}

			newACList.clear();
			newACnum = 0;
			pot = 0;
			j = 0;

			if (newACnum>maxNeighborNum){
				RadiusOfAC *= 0.949;
			} 
		}
	}


	// i need to see the changed members between old and new neighbor lists

	// find the gained members
	std::cout <<  "finding new neighbors...\n" << std::flush;

	for (j=0; j<newACnum; j++){

		isGain = true;

		for (Particle* element : ACList){
			if (element->PID == particle[j]->PID) {
				isGain = false;
			}
		}

		if (isGain) {

			rij2 = 0.0;
			drdv = 0.0;

			for (int dim=0; dim<Dim; dim++) {

				dx[dim] = particle[j]->PredPosition[dim] - PredPosition[dim];
				dv[dim] = particle[j]->PredVelocity[dim] - PredVelocity[dim];

				rij2 += dx[dim]*dx[dim];
				drdv += dx[dim]*dv[dim];
			}

			rij2 += EPS2;  // add softening length
			dr2i = 1.0/rij2;
			dr3i = particle[j]->Mass*dr2i*sqrt(dr2i);
			drdp = 3.0*drdv*dr2i;

			for (int dim=0; dim<Dim; dim++){
				dfirr[dim] += dx[dim]*dr3i;
				dfirrdot[dim] += (dv[dim] - dx[dim]*drdp)*dr3i;
			}

		}
	}

	// find the lost members

	for (Particle* element : ACList){

		isLoss = true;

		for (j=0; j<newACnum; j++){
			if (particle[j]->PID == element->PID) {
				isLoss = false;
			}
		}

		if (isLoss) {

			rij2 = 0.0;
			drdv = 0.0;

			for (int dim=0; dim<Dim; dim++) {

				dx[dim] = particle[j]->PredPosition[dim] - PredPosition[dim];
				dv[dim] = particle[j]->PredVelocity[dim] - PredVelocity[dim];

				rij2 += dx[dim]*dx[dim];
				drdv += dx[dim]*dv[dim];
			}

			rij2 += EPS2;  // add softening length
			dr2i = 1.0/rij2;
			dr3i = particle[j]->Mass*dr2i*sqrt(dr2i);
			drdp = 3.0*drdv*dr2i;

			for (int dim=0; dim<Dim; dim++){
				dfirr[dim] -= dx[dim]*dr3i;
				dfirrdot[dim] -= (dv[dim] - dx[dim]*drdp)*dr3i;
			}

		}
	}

	// calculated the final corrected forces

	if (NumberOfAC != 0) {
		dtsq = dtr*dtr;
		dt6 = 6.0/(dtr*dtsq); // to SY is it dtsq? it was dts1.
		dt2 = 2.0/dtsq;
		dtsq12 = dtsq/12.0;
		dt13 = dtr/3.0;

		for (int dim=0; dim<Dim; dim++) {

			dfr = FReg[dim] - freg[dim] - dfirr[dim];
			fdr0 = fregdot[dim] + dfirrdot[dim];
			frd = FRegDot[dim];
			sum = frd + fdr0;
			at3 = 2.0*dfr + dtr*sum;
			bt2 = -3.0*dfr - dtr*(sum+frd);

			Position[dim] = Position[dim] + (0.6*at3 + bt2)*dtsq12;
			Velocity[dim] = Velocity[dim] + (0.75*at3 + bt2)*dt13;

			FIrr[dim] += dfirr[dim];
			FReg[dim] = freg[dim];
			FIrrDot[dim] += dfirrdot[dim];
			FRegDot[dim] += fregdot[dim];

			dFIrr[dim][0] = FIrr[dim];
			dFReg[dim][0] = FReg[dim]; 
			dFReg[dim][1] = FRegDot[dim];
			dFReg[dim][2] = (3.0*at3 + bt2)*dt2;
			dFReg[dim][3] = at3*dt6;

			Force[dim] = 0.5*(FIrr[dim] + FReg[dim]);
			ForceDot[dim] = (FIrrDot[dim] + FRegDot[dim]);

		}
	}
	else {
		dtsq   = TimeStepReg*TimeStepReg;
		dt6    = 6.0/(TimeStepReg*dtsq); // to SY is it dtsq? it was dts1.
		dt2    = 2.0/dtsq;
		dtsq12 = dtsq/12.0;
		dt13   = TimeStepReg/3.0;

		for (int dim=0; dim<Dim; dim++) {
			dfr  = FReg[dim] - freg[dim];
			fdr0 = fregdot[dim];
			frd  = FRegDot[dim];
			sum  = frd + fdr0;
			at3  = 2.0*dfr + dtr*sum;
			bt2  = -3.0*dfr - dtr*(sum+frd);

			Position[dim] = PredPosition[dim] + (0.6*at3 + bt2)*dtsq12;
			Velocity[dim] = PredVelocity[dim] + (0.75*at3 + bt2)*dt13;

			FReg[dim] = freg[dim];
			FRegDot[dim] += fregdot[dim];

			dFReg[dim][0] = FReg[dim];
			dFReg[dim][1] = FRegDot[dim];
			dFReg[dim][2] = (3.0*at3 + bt2)*dt2;
			dFReg[dim][3] = at3*dt6;

			Force[dim] = 0.5*FReg[dim];
			ForceDot[dim] = FRegDot[dim];
		}
	}


	// update the new neighbor list

	NumberOfAC = newACnum;
	ACList.clear();

	for (j=0; j<newACnum; j++){
		ACList.push_back(particle[j]);
	}

	// update the regular time step

	CurrentTimeReg += TimeStepReg;
	this->calculateTimeStepReg(FReg,dFReg);
	this->isRegular = 0;
}



