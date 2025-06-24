#ifndef AR_INTERACTION_H
#define AR_INTERACTION_H
#ifdef FEWBODY
//#pragma once

extern double EnzoTimeStep;
extern FILE* workerout;
extern Particle *particles;

#include "ar_perturber.hpp"
#include <cassert>

#define ASSERT(x) assert(x)

//! AR interaction clas
class Interaction{
public:

    Float gravitational_constant;
    // TwoBodyTide tide;

    Interaction(): gravitational_constant(Float(1.0)) {}


    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(gravitational_constant>0.0);
        return true;
    }        

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"G      : "<<gravitational_constant<<std::endl;
    }    

    //! (Necessary) calculate inner member acceleration, potential and inverse time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _f1: force for particle 1 to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard are overwritten, not accummulating old values)
      @param[out] _f2: force for particle 2
      @param[out] _epot: total inner potential energy
      @param[in] _p1: particle 1
      @param[in] _p2: particle 2
      \return the inverse time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInvTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const Particle& _p1, const Particle& _p2) {
        // acceleration
        const Float mass1 = _p1.Mass;
        const Float* pos1 = _p1.Position;

        const Float mass2 = _p2.Mass;
        const Float* pos2 = _p2.Position;

        Float gm1 = gravitational_constant*mass1;
        Float gm2 = gravitational_constant*mass2;
        Float gm1m2 = gm1*mass2;
        
        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _f1.acc_in;
        Float* acc2 = _f2.acc_in;


        Float gmor3_1 = gm2*inv_r3;
        Float gmor3_2 = gm1*inv_r3;

        Float gm1or =  gm1*inv_r;
        Float gm2or =  gm2*inv_r;
        Float gm1m2or = gm1m2*inv_r;


        acc1[0] = gmor3_1 * dr[0];
        acc1[1] = gmor3_1 * dr[1];
        acc1[2] = gmor3_1 * dr[2];

        _f1.pot_in = -gm2or;

        acc2[0] = - gmor3_2 * dr[0];
        acc2[1] = - gmor3_2 * dr[1];
        acc2[2] = - gmor3_2 * dr[2];

        _f2.pot_in = -gm1or;


#ifdef AR_TTL 

        Float gm1m2or3 = gm1m2*inv_r3;
        Float* gtgrad1 = _f1.gtgrad;
        Float* gtgrad2 = _f2.gtgrad;
        gtgrad1[0] = gm1m2or3 * dr[0];
        gtgrad1[1] = gm1m2or3 * dr[1];
        gtgrad1[2] = gm1m2or3 * dr[2];

        gtgrad2[0] = - gtgrad1[0];
        gtgrad2[1] = - gtgrad1[1];
        gtgrad2[2] = - gtgrad1[2];
#endif

        // potential energy
        _epot = - gm1m2or;

        // transformation factor for kick
        Float gt_kick_inv = gm1m2or;
        // fprintf(stderr, "gt_kick_inv: %e\n", gt_kick_inv); // Eunwoo debug

        return gt_kick_inv;
    }

    //! calculate inner member acceleration, potential and inverse time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the inverse time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInv(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        _epot = Float(0.0);
        Float gt_kick_inv = Float(0.0);

        for (int i=0; i<_n_particle; i++) {
            const Float massi = _particles[i].Mass;
            const Float* posi = _particles[i].Position;
            Float* acci = _force[i].acc_in;
            acci[0] = acci[1] = acci[2] = Float(0.0);

#ifdef AR_TTL 
            Float* gtgradi = _force[i].gtgrad;
            gtgradi[0] = gtgradi[1] = gtgradi[2] = Float(0.0);
#endif

            Float poti = Float(0.0);
            Float gtki = Float(0.0);

            for (int j=0; j<_n_particle; j++) {
                if (i==j) continue;
                const Float massj = _particles[j].Mass;
                const Float* posj = _particles[j].Position; 
                Float dr[3] = {posj[0] -posi[0],
                               posj[1] -posi[1],
                               posj[2] -posi[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;

                Float gmor3 = gravitational_constant*massj*inv_r3;
                Float gmor = gravitational_constant*massj*inv_r;

                acci[0] += gmor3 * dr[0];
                acci[1] += gmor3 * dr[1];
                acci[2] += gmor3 * dr[2];

#ifdef AR_TTL                     
                Float gmimjor3 = massi*gmor3;
                gtgradi[0] += gmimjor3 * dr[0];
                gtgradi[1] += gmimjor3 * dr[1];
                gtgradi[2] += gmimjor3 * dr[2];
#endif

                poti -= gmor;
                gtki += gmor;
                    
            }
            _epot += poti * massi;
            gt_kick_inv += gtki * massi;
        }
        _epot   *= 0.5;
        gt_kick_inv *= 0.5;
        // fprintf(stderr, "gt_kick_inv: %e\n", gt_kick_inv); // Eunwoo debug

        return gt_kick_inv;
    }

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*!@param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
    */
    void calcAccPert(AR::Force* _force, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        static const Float inv3 = 1.0 / 3.0;

        // perturber force
        // const int n_pert = _perturber.neighbor_address.getSize();
        const int n_pert = _particle_cm.NumberOfNeighbor;
        // const int n_pert_single = _perturber.n_neighbor_single;
        // const int n_pert_group = _perturber.n_neighbor_group;

        if (n_pert>0) {

            Float time = _time;

            // auto* pert_adr = _perturber.neighbor_address.getDataAddress();
            auto pert_adr = _particle_cm.Neighbors;

            Float xp[n_pert][3], xcm[3], m[n_pert];
            // ChangeOver* changeover[n_pert_single];
            // H4::NBAdr<Particle>::Group* ptclgroup[n_pert_group];

            // int n_single_count=0;
            // int n_group_count=0;
            for (int j=0; j<n_pert; j++) {
                // H4::NBAdr<Particle>::Single* pertj;
                Particle* pertj;
                pertj = &particles[pert_adr[j]];
                // int k; // index of predicted data
                // if (pert_adr[j].type==H4::NBType::group) {
                //     pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                //     k = n_group_count + n_pert_single;
                //     ptclgroup[n_group_count] = (H4::NBAdr<Particle>::Group*)pert_adr[j].adr;
                //     n_group_count++;
                // }
                // else {
                //     pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                //     k = n_single_count;
                //     changeover[n_single_count] = &pertj->changeover;
                //     n_single_count++;
                // }

                Float dt = time - pertj->CurrentTimeIrr*EnzoTimeStep;
                // ASSERT(dt>=0.0); // Eunwoo debug // Is this right?
                //ASSERT(dt>=-1e-7);
                xp[j][0] = pertj->Position[0] + dt*(pertj->Velocity[0] + 0.5*dt*(pertj->a_irr[0][0] + inv3*dt*pertj->a_irr[0][1]));
                xp[j][1] = pertj->Position[1] + dt*(pertj->Velocity[1] + 0.5*dt*(pertj->a_irr[1][0] + inv3*dt*pertj->a_irr[1][1]));
                xp[j][2] = pertj->Position[2] + dt*(pertj->Velocity[2] + 0.5*dt*(pertj->a_irr[2][0] + inv3*dt*pertj->a_irr[2][1]));


                m[j] = pertj->Mass;
            }
            // ASSERT(n_single_count == n_pert_single);
            // ASSERT(n_group_count == n_pert_group);

            Float dt = time - _particle_cm.CurrentTimeIrr*EnzoTimeStep;
            // ASSERT(dt>=0.0); // Eunwoo debug // Is this right?

            xcm[0] = _particle_cm.Position[0] + dt*(_particle_cm.Velocity[0] + 0.5*dt*(_particle_cm.a_irr[0][0] + inv3*dt*_particle_cm.a_irr[0][1]));
            xcm[1] = _particle_cm.Position[1] + dt*(_particle_cm.Velocity[1] + 0.5*dt*(_particle_cm.a_irr[1][0] + inv3*dt*_particle_cm.a_irr[1][1]));
            xcm[2] = _particle_cm.Position[2] + dt*(_particle_cm.Velocity[2] + 0.5*dt*(_particle_cm.a_irr[2][0] + inv3*dt*_particle_cm.a_irr[2][1]));


            Float acc_pert_cm[3]={0.0, 0.0, 0.0};
            Float mcm = 0.0;
            // if (_perturber.need_resolve_flag) {
            // calculate component perturbation
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                Float& pot_pert = _force[i].pot_pert;
                const auto& pi = _particles[i];
                // auto& chi = pi.changeover;
                acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
                pot_pert = 0.0;

                Float xi[3];
                xi[0] = pi.Position[0] + xcm[0];
                xi[1] = pi.Position[1] + xcm[1];
                xi[2] = pi.Position[2] + xcm[2];

                // single perturber
                for (int j=0; j<n_pert; j++) {
                    Float dr[3] = {xp[j][0] - xi[0],
                                   xp[j][1] - xi[1],
                                   xp[j][2] - xi[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    Float r  = sqrt(r2);
                    Float r3 = r*r2;
                    Float gm = gravitational_constant*m[j];
                    Float gmor3 = gm/r3;

                    acc_pert[0] += gmor3 * dr[0];
                    acc_pert[1] += gmor3 * dr[1];
                    acc_pert[2] += gmor3 * dr[2];

                }

                acc_pert_cm[0] += pi.Mass *acc_pert[0];
                acc_pert_cm[1] += pi.Mass *acc_pert[1];
                acc_pert_cm[2] += pi.Mass *acc_pert[2];

                mcm += pi.Mass;

            }
//#ifdef AR_DEBUG
//            ASSERT(abs(mcm-_particle_cm.mass)<1e-10);
//#endif
                
            // get cm perturbation (exclude soft pert)
            acc_pert_cm[0] /= mcm;
            acc_pert_cm[1] /= mcm;
            acc_pert_cm[2] /= mcm;

            // remove cm. perturbation
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                Float& pot_pert = _force[i].pot_pert;
                const auto& pi = _particles[i];
                acc_pert[0] -= acc_pert_cm[0]; 
                acc_pert[1] -= acc_pert_cm[1];        
                acc_pert[2] -= acc_pert_cm[2]; 
                
                pot_pert -= acc_pert[0]*pi.Position[0] + acc_pert[1]*pi.Position[1] + acc_pert[2]*pi.Position[2];

            }

        }
        else {
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
                _force[i].pot_pert = 0.0;
            }  
        }
    }
    
    //! (Necessary) calculate acceleration from internal members and perturbers
    /*! The Force class acc_pert should be updated
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential 
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
      \return perturbation energy to calculate slowdown factor
    */
    Float calcAccPotAndGTKickInv(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        // inner force
        Float gt_kick_inv;
        if (_n_particle==2) gt_kick_inv = calcInnerAccPotAndGTKickInvTwo(_force[0], _force[1], _epot, _particles[0], _particles[1]);
        else gt_kick_inv = calcInnerAccPotAndGTKickInv(_force, _epot, _particles, _n_particle);

        calcAccPert(_force, _particles, _n_particle, _particle_cm, _perturber, _time);

        return gt_kick_inv;
    }


    //! calculate perturbation from c.m. acceleration
    Float calcPertFromForce(const Float* _force, const Float _mp, const Float _mpert) {
        Float force2 = _force[0]*_force[0]+_force[1]*_force[1]+_force[2]*_force[2];
#ifdef AR_SLOWDOWN_PERT_R4
        return force2/(gravitational_constant*_mp*_mpert);
#else
        Float force = sqrt(force2)/gravitational_constant;
        return sqrt(force/(_mp*_mpert))*force;
#endif
    }

    //! calculate perturbation from binary tree
    static Float calcPertFromBinary(const COMM::Binary& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        Float apo2 = apo*apo;
#ifdef AR_SLOWDOWN_PERT_R4
        return (_bin.m1*_bin.m2)/(apo2*apo2);
#else
        return (_bin.m1*_bin.m2)/(apo2*apo);
#endif
    }

    //! calculate perturbation from distance to perturber and masses of particle and perturber 
    static Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
        Float r2 = _r*_r;
#ifdef AR_SLOWDOWN_PERT_R4
        return _mp*_mpert/(r2*r2);
#else
        return (_mp*_mpert)/(r2*_r);
#endif
    }

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)

    //! calculate slowdown timescale
    void calcSlowDownTimeScale(Float& _t_min_sq, const Float dv[3], const Float dr[3], const Float& r, const Float& gm) {

        Float r2 = r*r;
        Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
        Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];

        Float semi = 1.0/(2.0/r - v2/gm);
        //hyperbolic, directly use velocity v
        if (semi<0) 
            _t_min_sq = std::min(_t_min_sq, r2/v2);
        else {
            Float ra_fact = (1 - r/semi); 
            Float e2 = drdv*drdv/(gm*semi) + ra_fact*ra_fact; // ecc^2
            Float r_vrmax = semi*(1-e2);
            if (r<r_vrmax) {
                // avoid decrese of vr once the orbit pass, calculate vr max at cos(E)=e (r==semi*(1-e^2))
                // vr_max = sqrt(er*(drdv^2*er + r*vcr2^2))/(G(m1+m2)r)
                //        = e*sqrt[G(m1+m2)/(a*(1-e^2)]
                Float vrmax_sq = e2*gm/r_vrmax;
                //Float rv2 = r*v2;
                //Float er = 2*gm - rv2;
                //Float vcr2 = gm - rv2;
                //Float vrmax_sq = er*(drdv*drdv*er + r*vcr2*vcr2)/(gm*gm*r2);
                _t_min_sq = std::min(_t_min_sq, semi*semi/vrmax_sq);
            }
            else {
                // r/vr
                Float rovr = r2/abs(drdv);
                _t_min_sq = std::min(_t_min_sq, rovr*rovr);
            }
        }
    }

    //! calculate slowdown perturbation and timescale from particle j to particle i
    /*! 
      @param[out] _pert_out: perturbation from particle j
      @param[out] _t_min_sq: timescale limit from particle j
      @param[in] _pi: particle i (cm of binary)
      @param[in] _pj: particle j 
     */
    void calcSlowDownPertOne(Float& _pert_out, Float& _t_min_sq, const Particle& pi, const Particle& pj) {
        Float dr[3] = {pj.Position[0] - pi.Position[0],
                       pj.Position[1] - pi.Position[1],
                       pj.Position[2] - pi.Position[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float r = sqrt(r2);
        _pert_out += calcPertFromMR(r, pi.Mass, pj.Mass);
            
#ifdef AR_SLOWDOWN_TIMESCALE
        Float dv[3] = {pj.Velocity[0] - pi.Velocity[0],
                       pj.Velocity[1] - pi.Velocity[1],
                       pj.Velocity[2] - pi.Velocity[2]};

        // identify whether hyperbolic or closed orbit
        Float gm = gravitational_constant*(pi.Mass+pj.Mass);

        calcSlowDownTimeScale(_t_min_sq, dv, dr, r, gm);
        // force dependent method
        // min sqrt(r^3/(G m))
        //Float gmor3 = (mp+mcm)*r*r2/(sdt->G*mp*mcm);
        //sdt->trf2_min =  std::min(sdt->trf2_min, gmor3);
#endif
    }

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[out] _pert_out: perturbation 
      @param[out] _t_min_sq: timescale limit 
      @param[in] _time: physical time for prediction
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(Float& _pert_out, Float& _t_min_sq, const Float& _time, const Particle& _particle_cm, const Perturber& _perturber) {
        static const Float inv3 = 1.0 / 3.0;

        // const int n_pert = _perturber.neighbor_address.getSize();
        const int n_pert = _particle_cm.NumberOfNeighbor;

        if (n_pert>0) {

            auto pert_adr = _particle_cm.Neighbors;

            Float xp[3], xcm[3];
            Float dt = _time - _particle_cm.CurrentTimeIrr*EnzoTimeStep;
            // ASSERT(dt>=0.0); // Eunwoo debug // Is this necessary?
            xcm[0] = _particle_cm.Position[0] + dt*(_particle_cm.Velocity[0] + 0.5*dt*(_particle_cm.a_irr[0][0] + inv3*dt*_particle_cm.a_irr[0][1]));
            xcm[1] = _particle_cm.Position[1] + dt*(_particle_cm.Velocity[1] + 0.5*dt*(_particle_cm.a_irr[1][0] + inv3*dt*_particle_cm.a_irr[1][1]));
            xcm[2] = _particle_cm.Position[2] + dt*(_particle_cm.Velocity[2] + 0.5*dt*(_particle_cm.a_irr[2][0] + inv3*dt*_particle_cm.a_irr[2][1]));

            Float mcm = _particle_cm.Mass;
            // auto& chi = _particle_cm.changeover;

#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method 
            Float vp[3], vcm[3];

            vcm[0] = _particle_cm.Velocity[0] + dt*(_particle_cm.a_irr[0][0] + 0.5*dt*_particle_cm.a_irr[0][1]);
            vcm[1] = _particle_cm.Velocity[1] + dt*(_particle_cm.a_irr[1][0] + 0.5*dt*_particle_cm.a_irr[1][1]);
            vcm[2] = _particle_cm.Velocity[2] + dt*(_particle_cm.a_irr[2][0] + 0.5*dt*_particle_cm.a_irr[2][1]);
#endif

            for (int j=0; j<n_pert; j++) {
                Particle* pertj;
                pertj = &particles[pert_adr[j]];

                Float dt = _time - pertj->CurrentTimeIrr*EnzoTimeStep;
                // ASSERT(dt>=0.0); // Eunwoo debug // Is this necessary?
                xp[0] = pertj->Position[0] + dt*(pertj->Velocity[0] + 0.5*dt*(pertj->a_irr[0][0] + inv3*dt*pertj->a_irr[0][1]));
                xp[1] = pertj->Position[1] + dt*(pertj->Velocity[1] + 0.5*dt*(pertj->a_irr[1][0] + inv3*dt*pertj->a_irr[1][1]));
                xp[2] = pertj->Position[2] + dt*(pertj->Velocity[2] + 0.5*dt*(pertj->a_irr[2][0] + inv3*dt*pertj->a_irr[2][1]));

                Float mj = pertj->Mass;

                // auto& chj = pertj->changeover;

                Float dr[3] = {xp[0] - xcm[0],
                               xp[1] - xcm[1],
                               xp[2] - xcm[2]};

                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float r = sqrt(r2);
                // Float k  = ChangeOver::calcAcc0WTwo(chi, chj, r);
                // _pert_out += calcPertFromMR(r, mcm, k*mj);
                _pert_out += calcPertFromMR(r, mcm, mj);

#ifdef AR_SLOWDOWN_TIMESCALE
                // velocity dependent method 
                vp[0] = pertj->Velocity[0] + dt*(pertj->a_irr[0][0] + 0.5*dt*pertj->a_irr[0][1]);
                vp[1] = pertj->Velocity[1] + dt*(pertj->a_irr[1][0] + 0.5*dt*pertj->a_irr[1][1]);
                vp[2] = pertj->Velocity[2] + dt*(pertj->a_irr[2][0] + 0.5*dt*pertj->a_irr[2][1]);

                Float dv[3] = {vp[0] - vcm[0],
                               vp[1] - vcm[1],
                               vp[2] - vcm[2]};

                // identify whether hyperbolic or closed orbit
                Float gm = gravitational_constant*(mcm+mj);

                calcSlowDownTimeScale(_t_min_sq, dv, dr, r, gm);
#endif
            }
        }
        else {
            _pert_out = 0.0;
            _t_min_sq = 0.0;
        }

        // add soft perturbation
        _pert_out += _perturber.soft_pert_min;

    }
#endif

    //! (Necessary) modify one particle function
    /*!
      @param[in] _p: particle
      @param[in] _time_now: current time (not physical time, NB unit)
      @param[in] _time_end: required evolved time (not physical time, do not directly use, NB unit)
      \return 0: no modification; 1: modify mass; 2: modify mass and velocity; 3: mass become zero
     */
    // template <class Tparticle>
    int modifyOneParticle(Particle& _p, const Float& _time_now, const Float& _time_end) {
        return 0;
    }

  
// This function is used if manager->interrupt_detection_option > 0
//! (Necessary) modify the orbits and interrupt check 
    /*! check the inner left binary whether their separation is smaller than particle radius sum and become close, if true, set one component stauts to merger with cm mass and the other unused with zero mass. Return the binary tree address 
      @param[in] _bin_interrupt: interrupt binary information: adr: binary tree address; time_now: current physical time; time_end: integration finishing time; status: interrupt status: change, merge,none
      @param[in] _bin: binarytree to check iteratively
     */
    void modifyAndInterruptIter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin) {

        if (_bin_interrupt.status==AR::InterruptStatus::none) {
            auto* p1 = _bin.getLeftMember();
            auto* p2 = _bin.getRightMember();

            if(_bin.getMemberN()==2) {

                Float radius = mergerRadius(p1, p2); // p1->radius + p2->radius;

                if (p1->getBinaryInterruptState()== BinaryInterruptState::collisioncandidate && 
                    p2->getBinaryInterruptState()== BinaryInterruptState::collisioncandidate &&
                    (p1->time_check<_bin_interrupt.time_end || p2->time_check<_bin_interrupt.time_end) &&
                    (p1->getBinaryPairID()==p2->ParticleIndex||p2->getBinaryPairID()==p1->ParticleIndex)) {

                        p1->setBinaryInterruptState(BinaryInterruptState::collision);
                        p2->setBinaryInterruptState(BinaryInterruptState::collision);
                        _bin_interrupt.status = AR::InterruptStatus::merge;
                        _bin_interrupt.adr = &_bin;

                        // merge(_bin_interrupt, _bin);
                    }
                else {
                    Float semi, ecc, dr, drdv;
                    _bin.particleToSemiEcc(semi, ecc, dr, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
                    Float peri = semi*(1 - ecc);

                    // if (peri<radius && p1->getBinaryPairID()!=p2->PID&&p2->getBinaryPairID()!=p1->PID) { // original
                    if (peri<radius) { // Eunwoo modified
                        Float ecc_anomaly  = _bin.calcEccAnomaly(dr);
                        Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, ecc);
                        Float mean_motion  = sqrt(gravitational_constant*_bin.Mass/(fabs(_bin.semi*_bin.semi*_bin.semi))); 
                        Float t_peri = mean_anomaly/mean_motion;
                        if (drdv<0 && t_peri<_bin_interrupt.time_end-_bin_interrupt.time_now) {
                            fprintf(workerout, "Merger1. peri: %e pc, radius: %e pc\n", peri*position_unit, radius*position_unit);
                            fprintf(workerout, "PID: %d and %d\n", p1->PID, p2->PID);
                            fflush(workerout);

                            p1->setBinaryInterruptState(BinaryInterruptState::collision);
                            p2->setBinaryInterruptState(BinaryInterruptState::collision);
                            p1->setBinaryPairID(p2->ParticleIndex);
                            p2->setBinaryPairID(p1->ParticleIndex);
                            _bin_interrupt.status = AR::InterruptStatus::merge;
                            _bin_interrupt.adr = &_bin;

                            // merge(_bin_interrupt, _bin);
                        }
                            
                        else if (semi>0||(semi<0&&drdv<0)) {
                            p1->setBinaryPairID(p2->ParticleIndex);
                            p2->setBinaryPairID(p1->ParticleIndex);
                            p1->setBinaryInterruptState(BinaryInterruptState::collisioncandidate);
                            p2->setBinaryInterruptState(BinaryInterruptState::collisioncandidate);
                            p1->time_check = std::min(p1->time_check, _bin_interrupt.time_now + (drdv<0 ? t_peri : (_bin.period - t_peri)));
                            p2->time_check = std::min(p1->time_check, p2->time_check);
                            fprintf(workerout, "Merger2. peri: %e pc, radius: %e pc\n", peri*position_unit, radius*position_unit);
                            fprintf(workerout, "PID: %d and %d might merge soon!\n", p1->PID, p2->PID);
                            fflush(workerout);
                        }
                    }
                }
            }
            else {
                for (int k=0; k<2; k++) 
                    if (_bin.isMemberTree(k)) modifyAndInterruptIter(_bin_interrupt, *_bin.getMemberAsTree(k));
            }
        }
    }

// This function is used if manager->interrupt_detection_option > 0
//! Modify the orbits and interrupt check for Kepler solver 
    /*! check the inner left binary whether their separation is smaller than particle radius sum and become close, if true, set one component stauts to merger with cm mass and the other unused with zero mass. Return the binary tree address 
      @param[in] _bin_interrupt: interrupt binary information: adr: binary tree address; time_now: current physical time; time_end: integration finishing time; status: interrupt status: change, merge,none
      @param[in] _bin: binarytree to check iteratively
      @param[in] _dt: integration time
     */
    void modifyAndInterruptKepler(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, double _dt) {

        if (_bin_interrupt.status==AR::InterruptStatus::none) {
            auto* p1 = _bin.getLeftMember();
            auto* p2 = _bin.getRightMember();

            Float radius = mergerRadius(p1, p2);

            Float semi, ecc, dr, drdv;
            _bin.particleToSemiEcc(semi, ecc, dr, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
            Float peri = semi*(1 - ecc);
            Float ecc_anomaly  = _bin.calcEccAnomaly(dr); // 0 ~ pi
            Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, ecc); // 0 ~ pi
            if (drdv < 0)
                mean_anomaly = 2*M_PI - mean_anomaly;
            mean_anomaly = fmod(mean_anomaly + 2*M_PI, 2*M_PI);
            Float mean_motion  = sqrt(gravitational_constant*_bin.Mass/(fabs(_bin.semi*_bin.semi*_bin.semi))); 
            Float t_peri = abs(mean_anomaly/mean_motion); // always smaller than half the period
            Float period = 2*M_PI/mean_motion;

            if (peri < radius && _dt > t_peri) {
                fprintf(workerout, "Merger in KeplerSolver. peri: %e pc, radius: %e pc\n", peri*position_unit, radius*position_unit);
                fprintf(workerout, "PID: %d and %d\n", p1->PID, p2->PID);
                fflush(workerout);

                p1->setBinaryInterruptState(BinaryInterruptState::collision);
                p2->setBinaryInterruptState(BinaryInterruptState::collision);
                p1->setBinaryPairID(p2->ParticleIndex);
                p2->setBinaryPairID(p1->ParticleIndex);
                _bin_interrupt.status = AR::InterruptStatus::merge;
                _bin_interrupt.adr = &_bin;

                // merge(_bin_interrupt, _bin);
            }
        }
    }    

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    }

    Float mergerRadius(Particle *p1, Particle *p2) {

        Float radius = 0.0;

        if (p1->ParticleType == (Blackhole+SingleStar) && p2->ParticleType == (Blackhole+SingleStar)) {
            radius = (p1->radius > p2->radius) ? 3*p1->radius : 3*p2->radius; // r_ISCO == 3 * Schwartzschild raiuds
        }
        else if (p1->ParticleType == (Blackhole+SingleStar) && p2->ParticleType == (NormalStar+SingleStar)) {
            radius = 1.3*pow((p1->Mass + p2->Mass)/p2->Mass, 1./3)*p2->radius; // TDE radius
        }
        else if (p1->ParticleType == (NormalStar+SingleStar) && p2->ParticleType == (Blackhole+SingleStar)) {
            radius = 1.3*pow((p1->Mass + p2->Mass)/p1->Mass, 1./3)*p1->radius; // TDE radius
        }
        else if (p1->ParticleType == (NormalStar+SingleStar) && p2->ParticleType == (NormalStar+SingleStar)) {
            radius = p1->radius + p2->radius; // Sum of two stellar radius
        }

        return radius;
    }
    
};
#endif
#endif