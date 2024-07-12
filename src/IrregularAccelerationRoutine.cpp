#include <iostream>
#include "global.h"
#include <cassert>

#define no_chain_debug

std::vector<Particle*> ComputationList{};
int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
bool CreateComputationChain(std::vector<Particle*> &particle);
bool UpdateComputationChain(Particle* ptcl);
bool CreateComputationList(Particle* ptcl);
bool AddNewBinariesToList(std::vector<Particle*> &particle);
void BinaryAccelerationRoutine(double next_time, std::vector<Particle*> &particle);
void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle, double current_time, ULL current_block);



Particle *FirstComputation;
// 1. sort particles according to next irregular time
// 2. the sorted list should be update regularly
bool IrregularAccelerationRoutine(std::vector<Particle*> &particle)
{

	//fprintf(binout, "Starting irregular force\n");
#ifdef time_trace
	_time.irr_chain.markStart();
#endif
	//std::cout << "Creating a computation chain ...\n";

	bool first = true;
	Particle *ParticleForComputation;
	//AddNewBinariesToList(particle, particle);
	//std::cout << "Create Chain\n" << std::flush;
	if (CreateComputationChain(particle) == false) {
		std::cout << "No irregular particle to update ...\n";
		return true;
	}

	/*
	fprintf(stderr,"Entire IrrChain=");
	Particle* next=FirstComputation;
	while (next!= nullptr) {
		fprintf(stderr,"%.3e (%d), ", next->CurrentTimeIrr+next->TimeStepIrr, next->PID);
		next = next->NextParticleForComputation;
	}
	fprintf(stderr,"\n");
	*/


	//std::cout << "Calculating irregular force ...\n" << std::endl;
#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif
	// Caculating irregular acceleration


	while (CreateComputationList(FirstComputation) && ComputationList.size() != 0) {
			
		//std::cout << "List size=" << ComputationList.size() << std::endl;

#define binary
#ifdef binary
		//if (AddNewBinariesToList(ComputationList, particle) && ComputationList.size() == 0) {
		if (AddNewBinariesToList(particle) && ComputationList.size() == 0) {
			//std::cout << "No irregular particle to update afte binary formation." << std::endl;
			fprintf(stdout, "No irregular particle to update afte binary formation.\n");
			break;
		}
		//fprintf(binout, "in the irregular force\n");

		/*
			fprintf(stderr, "in irr, particle: ");
			for (Particle* ptcl: particle) {
				fprintf(stderr,"%d, ",ptcl->PID);	
			}
			fprintf(stderr,"\n");	
			*/


		//fflush(stdout);
		for (Particle* ptcl:particle) {
				if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
						|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
					fprintf(stdout, "before, myself = %d\n", ptcl->PID);
					fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
					fflush(stdout);
					//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}
		for (Particle* ptcl:particle) {
			if (ptcl->CurrentTimeIrr > 1) 
				fprintf(stderr, "outside before bin, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
		}

		if ((BinaryList.size()>0)) { //&(binary_time_prev != binary_time)) {
			fprintf(binout, "-------------------------------------\n");
			fprintf(binout, "irr_time = %e \n",
					ComputationList[0]->CurrentTimeIrr*1e10/1e6*EnzoTimeStep);
			fprintf(binout, "binary_time = %e \n",
					binary_time*1e10/1e6*EnzoTimeStep);


			fprintf(stdout, "Integrating Binaries ...\n");
			fprintf(binout, "Evolve.cpp: integrating binaries\n");
			fprintf(binout, "# of binaries = %d \n",int(BinaryList.size()));
			fflush(stdout);
			fflush(stderr);
			fflush(binout);
			BinaryAccelerationRoutine(
					ComputationList[0]->CurrentTimeIrr+ComputationList[0]->TimeStepIrr,
				 	particle);
			fprintf(stdout, "Finishing Binaries ...\n");
			fflush(stdout);
			fflush(binout);

		}
		for (Particle* ptcl:particle) {
			if (ptcl->CurrentTimeIrr > 1)  {
				fprintf(stderr, "outside after bin, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
				fflush(stderr);
			}
		}

		for (Particle* ptcl:particle) {
				if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
						|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
					fprintf(stdout, "after, myself = %d\n", ptcl->PID);
					fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
					fflush(stdout);
					//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}
#endif

		//std::cout << "Start IRR\n" << std::flush;
		//while (ParticleForComputation != nullptr) {

#ifdef chain_debug
		std::cerr << "BU, ComputationList=" << std::flush;
		for (Particle* ptcl:ComputationList) {
			fprintf(stderr,"%.3e (%d), ", ptcl->CurrentTimeIrr+ptcl->TimeStepIrr, ptcl->PID);
		}
		std::cerr << std::endl;

		fprintf(stderr,"BU,IrrChain=");
		Particle* next=FirstComputation;
		while (next!= nullptr) {
			fprintf(stderr,"%.3e (%d), ", next->CurrentTimeIrr+next->TimeStepIrr, next->PID);
			next = next->NextParticleForComputation;
		}
		fprintf(stderr,"\n");

#endif

		for (Particle* ptcl:ComputationList) {

			//if(ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6 < 1e-7) {
			if (ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg < ptcl->CurrentBlockIrr)
				fprintf(stdout,"--------------------error--------------------------------------------------------------------\n");

			fprintf(stdout, "PID=%d, NextRegTimeBlock= %llu, NextIrrBlock = %llu (%.2e Myr)\n",
					ptcl->getPID(), NextRegTimeBlock,
					ptcl->CurrentBlockIrr+ptcl->TimeBlockIrr,
					(ptcl->CurrentTimeIrr + ptcl->TimeStepIrr)*EnzoTimeStep*1e10/1e6);
			fprintf(stdout, "CurrentTimeIrr = %.2e Myr (%llu), CurrentTimeReg = %.2e Myr (%llu), TimeStepIrr = %.2e Myr (%llu), TimeStepReg= %.2e Myr (%llu)\n",
					ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6, ptcl->CurrentBlockIrr,
					ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6, ptcl->CurrentBlockReg,
					ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6, ptcl->TimeBlockIrr, 
					ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6, ptcl->TimeBlockReg
					); 

			if (ptcl->CurrentBlockReg > ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg < ptcl->CurrentBlockIrr)
				fprintf(stdout,"----------------------------------------------------------------------------------------\n");
			//}
			/*
			for (Particle* ptcl:particle) {
				if (ptcl->CurrentTimeIrr > 1) {
					fprintf(stderr, "inside, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
					fflush(stderr);
				}
			}
			*/
			/*
				 fprintf(stdout, "a_tot = (%.2e,%.2e,%.2e), a_irr = (%.2e,%.2e,%.2e), a1_irr = (%.2e,%.2e,%.2e), a2_irr = (%.2e,%.2e,%.2e), a3_irr = (%.2e,%.2e,%.2e), n_n=%d\n",
				 ptcl->a_tot[0][0],
					ptcl->a_tot[1][0],
					ptcl->a_tot[2][0],
					ptcl->a_irr[0][0],
					ptcl->a_irr[1][0],
					ptcl->a_irr[2][0],
					ptcl->a_irr[0][1],
					ptcl->a_irr[1][1],
					ptcl->a_irr[2][1],
					ptcl->a_irr[0][2],
					ptcl->a_irr[1][2],
					ptcl->a_irr[2][2],
					ptcl->a_irr[0][3],
					ptcl->a_irr[1][3],
					ptcl->a_irr[2][3],
					ptcl->NumberOfAC);
			fprintf(stdout, "x = (%.2e,%.2e,%.2e), v = (%.2e,%.2e,%.2e)\n",
					ptcl->Position[0],ptcl->Position[1],ptcl->Position[2],
					ptcl->Velocity[0],ptcl->Velocity[1],ptcl->Velocity[2]
					);
			fflush(stdout);
			*/


			/*
			fprintf(stdout, "in irr, PID=%d : ", ptcl->PID);
			for (Particle* nn:ptcl->ACList) {
				fprintf(stdout,"%d, ",nn->PID);	
			}
			fprintf(stdout,"\n");	
			*/

#ifdef time_trace
			_time.irr_force.markStart();
#endif
			ptcl->calculateIrrForce(); // this includes particle position

#ifdef time_trace
			_time.irr_force.markEnd();
			_time.irr_force.getDuration();
#endif


			//ParticleForComputation = ParticleForComputation->NextParticleForComputation;// this includes time evolution.
			//ParticleForComputation = SortComputationChain(ParticleForComputation);
		}




#ifdef time_trace
		_time.irr_sort.markStart();
#endif
		//std::cout << "Update and Sort\n" << std::flush;
		//FirstComputation =  ComputationList.back()->NextParticleForComputation;
		/*
		if (FirstComputation != nullptr) {
			fprintf(stdout, "FirstComputation: PID=%d, NextBlockIrr=%llu\n", 
					FirstComputation->PID, FirstComputation->CurrentBlockIrr + FirstComputation->TimeBlockIrr);
			fprintf(stderr, "FirstComputation: PID=%d, NextTimeIrr=%.3e\n", 
					FirstComputation->PID, (FirstComputation->CurrentBlockIrr + FirstComputation->TimeBlockIrr)*time_step);
		}
		else
			fprintf(stdout, "FirstComputation is nullptr\n");
			*/
		// see if binary should be terminated or not 
		
		/*
		fprintf(stderr,"AC,BU,IrrChain=");
		next=FirstComputation;
		while (next!= nullptr) {
			fprintf(stderr,"%.3e (%d), ", next->CurrentTimeIrr+next->TimeStepIrr, next->PID);
			next = next->NextParticleForComputation;
		}
		fprintf(stderr,"\n");
		*/

		/*
		for (Binary* ptclBin: BinaryList) {
			if ((ptclBin->r>(ptclBin->r0*2.0)) || (ptclBin->TimeStep > 2.0*KSTime)) {
				std::cout << "Terminating Binary ..." << std::endl;
				fprintf(binout, "Terminating Binary at time : %e \n", binary_time);
				KSTermination(ptclBin->ptclCM, particle, binary_time, binary_block);
			}
		}
		*/


		//fflush(stdout);
		for (Particle* ptcl:particle) {
			if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
					|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
				fprintf(stdout, "before term, myself = %d\n", ptcl->PID);
				fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
				fflush(stdout);
				//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}

		// update particles and chain
		// The next particle of the particle calculated lastly should be the start of the next iteration.
#ifdef binary
		binary_time_prev = binary_time;
		binary_block = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
		binary_time  = binary_block*time_step;
		global_time_irr = ComputationList[0]->CurrentBlockIrr+ComputationList[0]->TimeBlockIrr;
		//std::cout << "ComputationList of " << ComputationList.size() << " : " ;
		bool bin_termination=false;	
#endif
		for (Particle* ptcl:ComputationList) {

			if (ptcl->CurrentTimeIrr > 1) 
				fprintf(stderr, "outside, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
			//std::cout << ptcl->PID << " " ;
			ptcl->updateParticle();
			ptcl->CurrentBlockIrr += ptcl->TimeBlockIrr;
			ptcl->CurrentTimeIrr   = ptcl->CurrentBlockIrr*time_step;
			//std::cout << "before TimeStepCal\n" << std::flush;
			ptcl->calculateTimeStepIrr(ptcl->a_tot, ptcl->a_irr);
			//std::cout << "after TimeStepCal\n" << std::flush;
#ifdef binary
			if (ptcl->isCMptcl) {
				if (ptcl->BinaryInfo->r > ptcl->BinaryInfo->r0*2.0
						|| ptcl->BinaryInfo->TimeStep*EnzoTimeStep*1e4 > 2.0*KSTime) {
					fprintf(binout, "Terminating Binary at time : %e \n", binary_time);
					fprintf(stdout, "Terminating Binary at time : %e \n", binary_time);
					KSTermination(ptcl, particle, binary_time, binary_block);
					bin_termination=true;
					continue;
				}
			}
#endif
			UpdateComputationChain(ptcl);
			if (ptcl->CurrentTimeIrr > 1) {
				fprintf(stderr, "outside after, PID=%d, CurrentTimeIrr=%e\n", ptcl->PID, ptcl->CurrentTimeIrr);
				fflush(stderr);
			}
		}

		for (Particle* ptcl:particle) {
			if (ptcl->Position[0] !=  ptcl->Position[0] || ptcl->Velocity[0] !=  ptcl->Velocity[0] 
					|| ptcl->a_tot[0][0] !=  ptcl->a_tot[0][0]) {
				fprintf(stdout, "after term, myself = %d\n", ptcl->PID);
				fprintf(stdout, "x[0] = %e\n", ptcl->Position[0]);
				fflush(stdout);
				//assert(ptcl->Position[0] ==  ptcl->Position[0]);
			}
		}
#ifdef binary
		if (bin_termination)  {
			//std::cerr << "After termination, NextRegTimeBlock=" << NextRegTimeBlock << std::endl; 
			CreateComputationChain(particle); // because NextRegTime might have been changed.
			/* 
				 fprintf(stderr, "in irr, particle: ");
				 for (Particle* ptcl: particle) {
				 fprintf(stderr,"%d, ",ptcl->PID);	
				 }
				 fprintf(stderr,"\n");	
				 */
		}
#endif



		//std::cout << "\ndone!" << std::endl ;

#ifdef chain_debug
		std::cerr << "AU, ComputationList=" << std::flush;
		for (Particle* ptcl:ComputationList) {
			fprintf(stderr,"%.3e (%d), ", ptcl->CurrentTimeIrr+ptcl->TimeStepIrr, ptcl->PID);
		}
		fprintf(stderr,"\n");

		fprintf(stderr,"AC,AU,IrrChain=");
		next=FirstComputation;
		while (next!= nullptr) {
			fprintf(stderr,"%.3e (%d), ", next->CurrentTimeIrr+next->TimeStepIrr, next->PID);
			next = next->NextParticleForComputation;
		}
		fprintf(stderr,"\n");
#endif


#define no_IRR_TEST
#ifdef IRR_TEST
		// create output at appropriate time intervals just for irr
		if (outputTime <= ComputationList[0]->CurrentTimeIrr ) {
			writeParticle(particle, ComputationList[0]->CurrentTimeIrr, outNum++);
			outputTime += outputTimeStep;
		}
#endif

#ifdef time_trace
		_time.irr_sort.markEnd();
		_time.irr_sort.getDuration();
#endif
		//fflush(stdout);
	}
	//std::cout << "Finishing irregular force ..." << std::endl;
	//kstd::cerr << "Finishing irregular force ..." << std::endl;
	//fflush(stderr);


	/*
		 for (Particle *ptcl : ComputationChain)
		 {
		 fprintf(stdout, "PID=%d, NextRegTime= %e, NextIrrTime = %e\n",
		 ptcl->getPID(), NextRegTime, ptcl->CurrentTimeIrr + ptcl->TimeStepIrr);
		 fprintf(stdout, "CurrentTimeIrr = %e, TimeStepIrr = %e, CurrentTimeReg=%e, TimeStepReg=%e\n",
		 ptcl->CurrentTimeIrr, ptcl->TimeStepIrr, ptcl->CurrentTimeReg, ptcl->TimeStepReg);

		 ptcl->calculateIrrForce(); // this includes particle position and time evolution.
		 SortComputationChain(ComputationChain);
		 }
		 */

	//std::cout << "Finishing irregular force ...\n" << std::endl;
	return true;
}



