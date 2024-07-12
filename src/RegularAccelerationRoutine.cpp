#include <iostream>
#include <cmath>
#include "global.h"

void UpdateNextRegTime(std::vector<Particle*> &particle);
void CalculateRegAccelerationOnGPU(std::vector<Particle*> RegularList, std::vector<Particle*> &particle);

bool RegularAccelerationRoutine(std::vector<Particle*> &particle)
{
	fprintf(stdout, "Calculating regular force ... particle size=%lu, regular size=%lu\n", particle.size(), RegularList.size());
	fflush(stdout);
	// Calulating regular acceleration of the particles
	if (RegularList.size() > 0) {
		for (Particle* ptcl: RegularList) {
			//if(ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6 < 1e-6) {
			if (ptcl->CurrentBlockReg >= ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg != ptcl->CurrentBlockIrr) {
				fprintf(stdout,"-------error--------------------------------------------------------------------------------- \n");
			}
			fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextRegTime= %.3e Myr(%llu),\n"\
					"dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n",
					ptcl->PID,
					ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
					ptcl->CurrentBlockIrr,
					ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
					ptcl->CurrentBlockReg,
					NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
					NextRegTimeBlock,
					ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
					ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
					ptcl->TimeBlockIrr,
					ptcl->TimeLevelIrr,
					ptcl->TimeBlockReg,
					ptcl->TimeLevelReg
					);
			if (ptcl->CurrentBlockReg >= ptcl->CurrentBlockIrr || ptcl->CurrentBlockReg+ptcl->TimeBlockReg != ptcl->CurrentBlockIrr) {
				fprintf(stdout,"--------------------------------------------------------------------------------------------- \n");
			}
			//}
			/*
			fprintf(stdout, "In reg, PID=%d : ", ptcl->PID);
			for (Particle* nn:ptcl->ACList) {
				fprintf(stdout,"%d, ",nn->PID);	
			}
			fprintf(stdout,"\n");	
			*/
			fprintf(stdout, " a_tot = (%.2e,%.2e,%.2e), a_reg = (%.2e,%.2e,%.2e), a_irr = (%.2e,%.2e,%.2e), n_n=%d, R=%.3e\n\
					a1_reg = (%.2e,%.2e,%.2e), a2_reg = (%.2e,%.2e,%.2e), a3_reg = (%.2e,%.2e,%.2e)\n\
					a1_irr = (%.2e,%.2e,%.2e), a2_irr = (%.2e,%.2e,%.2e), a3_irr = (%.2e,%.2e,%.2e)\n", 
					ptcl->a_tot[0][0],
					ptcl->a_tot[1][0],
					ptcl->a_tot[2][0],
					ptcl->a_reg[0][0],
					ptcl->a_reg[1][0],
					ptcl->a_reg[2][0],
					ptcl->a_irr[0][0],
					ptcl->a_irr[1][0],
					ptcl->a_irr[2][0],
					ptcl->NumberOfAC,
					ptcl->RadiusOfAC,
					ptcl->a_reg[0][1],
					ptcl->a_reg[1][1],
					ptcl->a_reg[2][1],
					ptcl->a_reg[0][2],
					ptcl->a_reg[1][2],
					ptcl->a_reg[2][2],
					ptcl->a_reg[0][3],
					ptcl->a_reg[1][3],
					ptcl->a_reg[2][3],
					ptcl->a_irr[0][1],
					ptcl->a_irr[1][1],
					ptcl->a_irr[2][1],
					ptcl->a_irr[0][2],
					ptcl->a_irr[1][2],
					ptcl->a_irr[2][2],
					ptcl->a_irr[0][3],
					ptcl->a_irr[1][3],
					ptcl->a_irr[2][3]
						);
		}
		fflush(stdout);
		CalculateRegAccelerationOnGPU(RegularList, particle);
		fflush(stdout);
	}

	/*
		 fprintf(stdout, "in reg, particle: ");
		 for (Particle* ptcl: particle) {
		 fprintf(stdout,"%d, ",ptcl->PID);	
		 }
		 fprintf(stdout,"\n");	
		 */


	global_time = NextRegTimeBlock*time_step;
	// update the next regular time step
	fprintf(stdout, "NextRegTimeBlock=%llu\n", NextRegTimeBlock);
	fprintf(binout, "Finishing regular force\n");
	std::cerr << "Finishing regular force ...\n" << std::flush;
	//fprintf(stderr, "Finishing regular force.. NextRegTime=%.3e\n", NextRegTimeBlock*time_step);
	fprintf(stdout, "Finishing regular force.. NextRegTime=%.3e\n", NextRegTimeBlock*time_step);
	return true;
}



void UpdateNextRegTime(std::vector<Particle*> &particle) {

	ULL time_tmp=0, time=block_max;

	RegularList.clear();
	for (Particle *ptcl:particle)
	{
		// Next regular time step
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;

		// Find the minum regular time step
		if (time_tmp <= time) {
			//fprintf(stderr, "PID=%d, time_tme=%llu\n", ptcl->PID, time_tmp);
			if ( time_tmp < time) {
				RegularList.clear();
				time = time_tmp;
			}
			RegularList.push_back(ptcl);
		}
	}
	NextRegTimeBlock = time;
}


