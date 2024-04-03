#include <vector>
#include <iostream>
#include <cmath>
#include "../global.h"
#include "../defs.h"



void direct_sum(double *x, double *v, double r2, double vx,
		        double mass, double (&a)[3], double (&adot)[3]);



//       /////////////////////        //
//       /////////////////////        //
//       CALCULATEACCELERATION        //
//       /////////////////////        //
//       /////////////////////        //


void CalculateKSAcceleration(Particle* ptclI, Particle* ptclJ, Binary* ptclBin, std::vector<Particle*> &particle, double current_time) {

    int j=0;
    Particle* ptcl1;
	double x[Dim], v[Dim], a[Dim], adot[Dim];
	double r2 = 0;
	double vx = 0;
    double v2 = 0;
    double m_r3, vx_r2, v2x2_r4,v2_r2__ax_r2__v2x2_r4, a2dot, a3dot;
    double A,B;

	for (int dim=0; dim<Dim; dim++) {
		x[dim]    = 0.;
		v[dim]    = 0.;
		a[dim]    = 0.;
		adot[dim] = 0.;
	}

	//std::cout << "nbody+: Entering CalculateInitialAcceleration  ..." << std::endl;


    for (int i=0; i<2; i++) {

        if (i=0) {
            ptcl1 = ptclI;

        } else {
            ptcl1 = ptclJ;
        }

        // updated the predicted positions and velocities just in case
            
        ptcl1->predictParticleSecondOrder(current_time);

        for (Particle *ptcl2: particle) {

            r2 = 0;
            vx = 0;
            v2 = 0;


            // if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
            ptcl2->predictParticleSecondOrder(current_time);

            for (int dim=0; dim<Dim; dim++) {
                x[dim] = ptcl1->PredPosition[dim] - ptcl2->PredPosition[dim];
                v[dim] = ptcl1->PredVelocity[dim] - ptcl2->PredVelocity[dim];
                r2    += x[dim]*x[dim];
                vx    += v[dim]*x[dim];
                v2    += v[dim]*v[dim];
            }

            //r2  += EPS2;
            m_r3 = ptcl2->Mass/r2/sqrt(r2); 

            for (int dim=0; dim<Dim; dim++) {
                // Calculate 0th and 1st derivatives of acceleration
                if ((ptcl1->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
                    ptcl1->a_reg[dim][0] += m_r3*x[dim];
                    ptcl1->a_reg[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
                }
                else {
                    ptcl1->a_irr[dim][0] += m_r3*x[dim];
                    ptcl1->a_irr[dim][1] += m_r3*(v[dim] - 3*x[dim]*vx/r2);
                    j++;
                }
            }
    
        } // end of loop for ptcl2 (full particle)


        // update total acceleration as well

        for (int dim=0; dim<Dim; dim++)	 {
            for (int order=0; order<HERMITE_ORDER; order++) {
                ptcl1->a_tot[dim][order] = ptcl1->a_reg[dim][order] + ptcl1->a_irr[dim][order]; 
            }
	    }

    }// end of loop for pair particle, ptclI and ptclJ


    // copy the calculated values to CM particle
    // and initialize the 3rd and 4th order derivatives just in case

    for (int dim=0; dim<Dim; dim++)	 {
		for (int order=0; order<HERMITE_ORDER; order++) {
            if (order<2) {
			    ptclBin->a_reg[dim][order] = (ptclI->a_reg[dim][order]*ptclI->Mass + ptclJ->a_reg[dim][order]*ptclJ->Mass)/(ptclBin->Mass); 
                ptclBin->a_irr[dim][order] = (ptclI->a_irr[dim][order]*ptclI->Mass + ptclJ->a_irr[dim][order]*ptclJ->Mass)/(ptclBin->Mass); 
                ptclBin->a_tot[dim][order] = (ptclI->a_tot[dim][order]*ptclI->Mass + ptclJ->a_tot[dim][order]*ptclJ->Mass)/(ptclBin->Mass); 
            } else {
                ptclBin->a_reg[dim][order] = 0.0; 
                ptclBin->a_irr[dim][order] = 0.0;
                ptclBin->a_tot[dim][order] = 0.0;
            }
		}
	}


    // updated the predicted positions and velocities just in case
            
    //ptcl1 = ptclCM;
    //ptcl1->predictParticleSecondOrder(current_time);


    for (Particle *ptcl2: particle) {

            r2 = 0;
            vx = 0;
            v2 = 0;


            // if current time = the time we need, then PredPosition and PredVelocity is same as Position and Velocity
            ptcl2->predictParticleSecondOrder(current_time);

            for (int dim=0; dim<Dim; dim++) {
                x[dim] = ptclBin->PredPosition[dim] - ptcl2->PredPosition[dim];
                v[dim] = ptclBin->PredVelocity[dim] - ptcl2->PredVelocity[dim];
                r2    += x[dim]*x[dim];
                vx    += v[dim]*x[dim];
                v2    += v[dim]*v[dim];
            }

            //r2  += EPS2;
            m_r3 = ptcl2->Mass/r2/sqrt(r2);

            vx_r2   = vx/r2;
            v2x2_r4 = vx_r2*vx_r2;
            v2_r2__ax_r2__v2x2_r4 = (v2+x[0]*a[0]+x[1]*a[1]+x[2]*a[2])/r2+v2x2_r4;
            A = (9*(v[0]*a[0]+v[1]*a[1]+v[2]*a[2]) + 3*(x[0]*adot[0]+x[1]*adot[1]+x[2]*adot[2]))/r2\
                    +3*vx_r2*(3*v2_r2__ax_r2__v2x2_r4 - 4*v2x2_r4);

            for (int dim=0; dim<Dim; dim++) {
                B     = v[dim] - 3*x[dim]*vx_r2;
                a2dot = (v[dim] - 6*B*vx_r2                 - 3*v2_r2__ax_r2__v2x2_r4*x[dim])*m_r3;
                a3dot = (a[dim] - 9*B*v2_r2__ax_r2__v2x2_r4 - A*x[dim]                      )*m_r3\
                                - 9*vx_r2*a2dot;
                if ((ptclBin->NumberOfAC==0) || (ptcl2 != ptcl1->ACList[j])) {
                    ptclBin->a_reg[dim][2] += a2dot;
                    ptclBin->a_reg[dim][3] += a3dot;
                }
                else {
                    ptclBin->a_irr[dim][2] += a2dot;
                    ptclBin->a_irr[dim][3] += a3dot;
                    j++;
                }
            } // endfor dim

    } // end of loop for ptcl2 (full particle)

    for (int dim=0; dim<Dim; dim++)	 {
		for (int order=2; order<HERMITE_ORDER; order++) {
            ptclBin->a_tot[dim][order] = (ptclBin->a_reg[dim][order] + ptclBin->a_reg[dim][order]*ptclJ->Mass); 
		}
	}


}



//        ///////////////////        //
//        ///////////////////        //
//           ISKSCANDIDATE           //
//        ///////////////////        //
//        ///////////////////        //



// made 2024.02.19 by Seoyoung Kim

// from particles with shortest time steps...
// see if they meet the KS regularlization conditions
// reference - search.f

void Particle::isKSCandidate(double next_time) {

    // temporary calculation variables

    double x[Dim],v[Dim];
    double r2;

    double rmin;
    Particle* minPtcl;

    // t
    int numberOfCloseParticle;
    std::vector<Particle*> closeParticleList;



    // initialize variables and count numbers just in case

    r2 = 0.;
    rmin = 1e8;
    numberOfCloseParticle = 0;
    closeParticleList.clear();


    // predict the particle position to obtain more information
    // particle regularlized if the conditions are satisfied at a future time


    this->predictParticleSecondOrder(next_time);

    // need to consider case when c.m particles are the cause of the small steps
    // need to be added later - CHECK

    for (Particle* ptcl: ACList) {

        // if particle time step is too large, skip
        // if the neighbor step is larger then 8 times of the candidate particle, then skip
        
        if ((ptcl->TimeStepIrr) > (8*TimeStepIrr)) {
            continue;
        }


        // find out what the paired particle is

        ptcl->predictParticleSecondOrder(next_time);

        for (int dim=0; dim<Dim; dim++) {
				
            // calculate position and velocity differences
		    x[dim] = ptcl->PredPosition[dim] - this->PredPosition[dim];
			v[dim] = ptcl->PredVelocity[dim] - this->PredVelocity[dim];

			// calculate the square of radius and inner product of r and v for each case
			r2 += x[dim]*x[dim];
        }

        // find out the close particles

        if (r2<KSdistance) {

            numberOfCloseParticle += 1;
            closeParticleList.push_back(ptcl);

            if (r2<rmin) {
                rmin = r2;
                minPtcl = ptcl;
            }
        }

        // save the KS pair information

        if (numberOfCloseParticle>0) {
            isKS = true;
            KSPairParticle = minPtcl;
        }

    }    
}




//        ///////////////////        //
//        ///////////////////        //
//           STARTNEWKSREG           //
//        ///////////////////        //
//        ///////////////////        //



// start new KS regularlization



//        ///////////////////        //
//        ///////////////////        //
//        NEWKSINITIALIZATION        //
//        ///////////////////        //
//        ///////////////////        //



// initialize conditions for new KS regularlization
// reference: ksinit.F

void NewKSInitialization(Particle* ptclI, std::vector<Particle*> &particle, double current_time) {

    // basic variables for calculation

    Particle *ptclJ;
    Binary *ptclBin;
    std::vector<Particle*> KSNeighborList;


    // BASIC DEFINITIONS

    // define the pair particle for particle I

    ptclJ = ptclI->KSPairParticle;

    ptclI->predictParticleSecondOrder(current_time);
    ptclJ->predictParticleSecondOrder(current_time);


    // need to put option if there aren't any close neighbors


    // define the new center of mass particle

    ptclBin = new Binary;

    // calculate the values of the center of mass particle
    // and save it to the new center of mass particle

    ptclBin->Mass = ptclI->Mass + ptclJ->Mass;

    for (int dim=0; dim<Dim; dim++) {
        ptclBin->Position[dim] = (ptclI->PredPosition[dim]*ptclI->Mass + ptclJ->PredPosition[dim]*ptclJ->Mass)/ptclBin->Mass;
        ptclBin->Velocity[dim] = (ptclI->PredVelocity[dim]*ptclI->Mass + ptclJ->PredVelocity[dim]*ptclJ->Mass)/ptclBin->Mass;
        ptclBin->PredPosition[dim] = ptclBin->Position[dim];
        ptclBin->PredVelocity[dim] = ptclBin->Velocity[dim];
    }

    ptclBin->CurrentTime = current_time;
    ptclBin->PredTime = current_time;

    ptclBin->BinaryParticleI = ptclI;
    ptclBin->BinaryParticleI = ptclJ;


    // copy the neighbor list for c.m particle

    ptclBin->RadiusOfAC = ptclI->RadiusOfAC;

    for (Particle* ptcl: ptclI->ACList) {

        if (ptcl == ptclJ) {
            continue;
        }

        ptclBin->ACList.push_back(ptcl);
        ptclBin->NumberOfAC += 1;

    }


    // predict coordinates of neighbor particles

    for (Particle* ptcl: ptclI->ACList) {

        ptcl->predictParticleSecondOrder(current_time);

        // need to make routine for resolving KS components if a KS pair is in the neighbor

    }


    // calculate the 0th, 1st, 2nd, 3rd derivative of accleration accurately for the binary pair particle and the cm particle

    CalculateKSAcceleration(ptclI,ptclJ,ptclBin,particle,current_time);


    // calculate the initial values of relevant variables

    ptclBin->InitializeBinary(current_time);

}

