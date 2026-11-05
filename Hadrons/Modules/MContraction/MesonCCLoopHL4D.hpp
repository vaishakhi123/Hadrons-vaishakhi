/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/MesonCCLoopHL.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>
Author: Vera Guelpers <vmg1n14@soton.ac.uk>
Author: Tom Blum (conserved currents)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MContraction_MesonLoopCCHL4D_hpp_
#define Hadrons_MContraction_MesonLoopCCHL4D_hpp_

#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/Solver.hpp>
#include <Grid/lattice/Lattice_reduction.h>
#include <ctime>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             TMesonLoopCCHL                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


class MesonLoopCCHL4DPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonLoopCCHL4DPar,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, eigenPack,
                                    std::string, solver,
                                    std::string, action,
                                    double, mass,
                                    int, tinc,
                                    int, block,
                                    int, tblock,
                                    int, hits);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonLoopCCHL4D: public Module<MesonLoopCCHL4DPar>
{
public:
    typedef typename FImpl1::FermionField FermionField;
    typedef A2AVectorsSchurStaggered<FImpl1> A2A;
    typedef FermionOperator<FImpl1>          FMat;
    FERM_TYPE_ALIASES(FImpl1,1);
    FERM_TYPE_ALIASES(FImpl2,2);
    SOLVER_TYPE_ALIASES(FImpl1,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<Complex>, corr,
                                        std::vector<std::vector<Complex>>, corr_all);
    };
public:
    // constructor
    TStagMesonLoopCCHL4D(const std::string name);
    // destructor
    virtual ~TStagMesonLoopCCHL4D(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    protected:
    // execution
    virtual void setup(void);
    // execution
    virtual void execute(void);
    inline bool exists (const std::string& name) {
      struct stat buffer;
      return (stat (name.c_str(), &buffer) == 0);
    }

private:
    //FMat         *action_{nullptr};
    //Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(StagMesonLoopCCHL4D, ARG(TStagMesonLoopCCHL4D<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonLoopCCHL implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonLoopCCHL4D<FImpl1, FImpl2>::TStagMesonLoopCCHL4D(const std::string name)
: Module<MesonLoopCCHL4DPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL4D<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().action};
    std::string sub_string;
    
    input.push_back(par().eigenPack);
    //input.push_back(par().solver + "_subtract");
    input.push_back(par().solver);
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL4D<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};

    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL4D<FImpl1, FImpl2>::setup(void)
{
    
    auto        &action     = envGet(FMat, par().action);
    //auto        &solver     = envGet(Solver, par().solver + "_subtract");
    auto        &solver     = envGet(Solver, par().solver);
    envTmp(A2A, "a2a", 1, action, solver);
    
    envTmpLat(FermionField, "source");
    envTmpLat(FermionField, "sink");
    envTmpLat(FermionField, "tmp");
    envTmpLat(FermionField, "tmp2");
    envTmpLat(FermionField, "tmp3");
    envTmpLat(FermionField, "sol");
    envTmpLat(FermionField, "solshift");
    envTmpLat(FermionField, "sourceshift");
    envTmpLat(FermionField, "w");
    envTmpLat(LatticeComplex, "eta");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL4D<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing High-Low Conserved Current Stag meson contractions " << std::endl;

    std::vector<ComplexD>  corr;
    std::vector<Result> result(3);
    int nt = env().getDim(Tp);
    int ns = env().getDim(Xp);
    // init
    for(int mu=0;mu<3;mu++){
        result[mu].corr.resize(nt);
        result[mu].corr_all.resize(nt, std::vector<ComplexD>(nt, ComplexD(0.,0.)));
        for(int t=0;t<nt;t++){
            result[mu].corr[t]=(ComplexD)(0.,0.);
            for (int dt = 0; dt < nt; dt++) {
                    result[mu].corr_all[t][dt] = (ComplexD)(0.,0.);
            }
            
        }
    }

    auto &U       = envGet(LatticeGaugeField, par().gauge);
    auto        &action      = envGet(FMat, par().action);
    //auto        &solver    = envGet(Solver, par().solver + "_subtract");
    auto        &solver    = envGet(Solver, par().solver);
    auto &epack   = envGet(BaseFermionEigenPack<FImpl1>, par().eigenPack);
    
    double mass = par().mass;
    int block = par().block;
    int hits = par().hits;

    int Nl_ = epack.evec.size();
    std::vector<double> mlsq(Nl_);
    
    FermionField tmp_e(env().getRbGrid());
    FermionField tmp_o(env().getRbGrid());

    envGetTmp(A2A, a2a);
    
    // Do spatial gammas only
    Lattice<iScalar<vInteger>> x(U.Grid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger>> y(U.Grid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger>> z(U.Grid()); LatticeCoordinate(z,2);
    Lattice<iScalar<vInteger>> t(U.Grid()); LatticeCoordinate(t,3);
    Lattice<iScalar<vInteger>> lin_z(U.Grid()); lin_z=x+y;
    Lattice<iScalar<vInteger>> lin_t(U.Grid()); lin_t=x+y+z;
    LatticeComplex phases(U.Grid());
    std::vector<LatticeColourMatrix> Umu(3,U.Grid());
    
    // source, solution
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sink);
    envGetTmp(FermionField, tmp);
    envGetTmp(FermionField, tmp2);
    envGetTmp(FermionField, tmp3);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, solshift);
    envGetTmp(FermionField, sourceshift);
    envGetTmp(FermionField, w);


    
    
    for(int i=0;i<Nl_;i++){
        mlsq[i]=(epack.eval[i]-mass*mass) * mass;
    }
    DeflatedGuesser<FermionField> LLsub(epack.evec, mlsq);
    
    std::string outFileName;
    std::vector<std::vector<std::vector<ComplexD>>> all_results(3, 
        std::vector<std::vector<ComplexD>>(hits,
            std::vector<ComplexD>(nt, ComplexD(0., 0.))));
    
    for(int mu=0;mu<3;mu++){

        //staggered phases go into links
        Umu[mu] = PeekIndex<LorentzIndex>(U,mu);
        phases=1.0;
        if(mu==0){
        }else if(mu==1){
            phases = where( mod(x    ,2)==(Integer)0, phases,-phases);
        }else if(mu==2){
            phases = where( mod(lin_z,2)==(Integer)0, phases,-phases);
        }else if(mu==3){
            phases = where( mod(lin_t,2)==(Integer)0, phases,-phases);
        }else assert(0);
        Umu[mu] *= phases;
    }

    
    std::vector<Complex> eta(nt*3*hits*Nl_*2);
    int tblock = par().tblock;
    std::vector<FermionField> source_list(tblock*3*hits, env().getGrid()); 
    std::vector<FermionField> sink_list(tblock*3*hits,env().getGrid());


    // Precompute and store random numbers for eta for all time slices, mu, hits and eigenvectors
    for(int ts = 0; ts < nt; ts++){
        for(int mu=0; mu< 3; mu++){
            for(int hit = 0; hit < hits; hit++){
                for(unsigned int il = 0; il < Nl_; il += block){
                    for(int iv = il; iv < il + block; iv++){
                        for(int pm = 0; pm < 2; pm++){
                            int idx = ts * hits * 3 * Nl_ * 2 + 
                                      mu * hits * Nl_ * 2 +
                                      hit * Nl_ * 2 +
                                      iv * 2 + pm;
                            bernoulli(rngSerial(), eta[idx]);
                            Complex shift(1., 1.);
                            eta[idx] = (2.*eta[idx] - shift)*(1./::sqrt(2.));
                        }
                    }
                }
            }
        }
    }
    
    FermionField sub(env().getGrid());

    // loop over time slice
    for(int its=0; its<nt;its+=tblock){
        for (int ts=its; ts< its+tblock; ts+=par().tinc){
            // loop over directions
            for(int mu=0;mu<3;mu++){
                
                // loop over hits
                for(int hit = 0; hit < hits; hit++)
                {
                        int s_idx = (ts-its)*3*hits + mu*hits + hit;
                        source_list[s_idx] = 0.0;
                        sink_list[s_idx] = 0.0;
                }
            }
        }
        for(int iv=0;iv<Nl_;iv++){
            std::complex<double> eval(mass,sqrt(epack.eval[iv]-mass*mass));
            for(int pm=0;pm<2;pm++){
                
                a2a.makeLowModeW(w, epack.evec[iv], eval, pm);
                 
                if(pm){
                    eval = conjugate(eval);
                }
                std::complex<double> iota_angle(0.0, std::arg(eval));
                // coefficients as Grid's ComplexD (std::complex on CPU, thrust::complex on GPU)
                std::complex<double> csrc = std::exp(-iota_angle)/std::sqrt(std::abs(eval));
                ComplexD             coefSrc(csrc.real(), csrc.imag());
                RealD                coefSnk = 1./std::sqrt(std::abs(eval));
                for (int ts=its; ts< its+tblock; ts+=par().tinc){
    
                    // lopp over directions
                    for(int mu=0;mu<3;mu++){
                    
                        for(int hit = 0; hit < hits; hit++)
                        {
                            int s_idx = (ts-its)*3*hits + mu*hits + hit;
                            int idx = ts * hits * 3 * Nl_ * 2 + 
                                      mu * hits * Nl_ * 2 +
                                      hit * Nl_ * 2 +
                                      iv * 2 + pm;
                
                            source_list[s_idx] += (eta[idx]*coefSrc)*w;
                            sink_list[s_idx] += (eta[idx]*coefSnk)*w;
                        }
                    }
                }
            }
        }
        // loop over time slices
        for (int ts=its; ts< its+tblock; ts+=par().tinc){
            
            LOG(Message) << "StagMesonLoopCCHLHL src_ts " << ts << std::endl;
            // loop over directions
            for(int mu=0;mu<3;mu++){
            

                LOG(Message) << "StagMesonLoopCCHLHL src_mu " << mu << std::endl;
                 
                for(int hit = 0; hit < hits; hit++){
                    int s_idx = (ts-its)*3*hits + mu*hits + hit;
                    
                    tmp = where(t == ts, source_list[s_idx], source_list[s_idx]*0.);
                    tmp2 = adj(Umu[mu]) * tmp;
                    
                    // shift source x-mu to x
                    tmp = Cshift(tmp2, mu, -1);
                    
                    tmp3 = tmp;

                    tmp = where(t == ts, source_list[s_idx], source_list[s_idx]*0.);
                    // shift source
                    tmp2 = Cshift(tmp, mu, 1);
                    tmp = Umu[mu] * tmp2;

                    tmp3 += tmp;

                    solver(sol, tmp3);

                    // subtract the low modes
                    sub = Zero();
                    pickCheckerboard(Even,tmp_e,tmp3);
                    action.Meooe(tmp_e,tmp_o);
                    LLsub(tmp_o,tmp_e);
                    action.Meooe(tmp_e,tmp_o);// tmp_o is now even
                    setCheckerboard(sub,tmp_o);
                    sol += sub;
                
             
                    // take inner-product with eigenbra on all time slices
                    solshift = Cshift(sol,mu,1);
                    solshift = Umu[mu]*solshift;
                    sliceInnerProductVector(corr,sink_list[s_idx],solshift,3); //first term + second term
                    
                    for(int tsnk=0; tsnk<nt; tsnk++){
                        result[mu].corr[(tsnk-ts+nt)%nt] += (corr[tsnk]);
                        result[mu].corr_all[ts][(tsnk - ts + nt) % nt] += (corr[tsnk]);
                    }
            
                    // take inner-product with eigenmode on all time slices
                        
                    sourceshift = Cshift(sink_list[s_idx],mu,1);
                    sourceshift = Umu[mu]*sourceshift;
                    sliceInnerProductVector(corr,sourceshift,sol,3); //fourth term
                    for(int tsnk=0; tsnk<nt; tsnk++){
                        result[mu].corr[(tsnk-ts+nt)%nt] += (corr[tsnk]); 
                        result[mu].corr_all[ts][(tsnk - ts + nt) % nt] += (corr[tsnk]);
                    }
                }
            }
        }
    }
    for (int i = 0; i < 3; ++i){
        if(U.Grid()->IsBoss()){
            makeFileDir(par().output);
            outFileName = par().output+"HLcc_2pt_mu"+std::to_string(i);
            for(int t=0; t<nt; t++)
                result[i].corr[t] = result[i].corr[t]/ComplexD(hits, 0.0);
            for (int ts = 0; ts < nt; ts++) {
                for (int dt = 0; dt < nt; dt++) {
                    result[i].corr_all[ts][dt] /= ComplexD(hits, 0.0);
                }
            }
            
            saveResult(outFileName, "HLCC", result[i]);
        }
    }
    
}
    
END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
