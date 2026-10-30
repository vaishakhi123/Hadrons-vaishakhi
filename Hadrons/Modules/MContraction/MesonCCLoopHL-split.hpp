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

#ifndef Hadrons_MContraction_MesonLoopCCHL_hpp_
#define Hadrons_MContraction_MesonLoopCCHL_hpp_

#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/Solver.hpp>
#include <Grid/lattice/Lattice_reduction.h>
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>
#include <utility>


BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             TMesonLoopCCHL                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


class MesonLoopCCHLPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonLoopCCHLPar,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, eigenPack,
                                    std::string, action,
                                    double, mass,
                                    std::string, solver,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    double     , c1,
                                    double     , tad,
                                    bool, solInitGuess , // true for making guess solution
                                    bool, subGuess , //true for subtract
                                    int, tinc,
                                    int, block,
                                    int, hits,
                                    std::string, mpi_split);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonLoopCCHL: public Module<MesonLoopCCHLPar>
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
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TStagMesonLoopCCHL(const std::string name);
    // destructor
    virtual ~TStagMesonLoopCCHL(void) {};
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

MODULE_REGISTER_TMP(StagMesonLoopCCHL, ARG(TStagMesonLoopCCHL<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonLoopCCHL implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonLoopCCHL<FImpl1, FImpl2>::TStagMesonLoopCCHL(const std::string name)
: Module<MesonLoopCCHLPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().action};
    std::string sub_string;
    
    input.push_back(par().eigenPack);
    //input.push_back(par().solver + "_subtract");
    input.push_back(par().solver);
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};

    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL<FImpl1, FImpl2>::setup(void)
{
    //Grid_init(&argc,&argv);
  
    auto        &action     = envGet(FMat, par().action);
    //auto        &solver     = envGet(Solver, par().solver + "_subtract");
    auto        &solver     = envGet(Solver, par().solver);
    envTmp(A2A, "a2a", 1, action, solver);
    
    envTmpLat(FermionField, "source");
    
    envTmpLat(FermionField, "sink");
    envTmpLat(FermionField, "tmp1");
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
void TStagMesonLoopCCHL<FImpl1, FImpl2>::execute(void)
{

    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    //Coordinate mpi_split (mpi_layout.size(),1);
    

    LOG(Message) << "Computing High-Low Conserved Current Stag-meson contractions " << std::endl;

    std::vector<ComplexD>  corr;
    std::vector<Result> result(3);
    int nt = env().getDim(Tp);
    int ns = env().getDim(Xp);
    // init
    for(int mu=0;mu<3;mu++){
        result[mu].corr.resize(nt);
        for(int t=0;t<nt;t++){
            result[mu].corr[t]=(ComplexD)(0.,0.);
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
    
    double c1 = par().c1;
    double tad = par().tad;
    
    
    LOG(Message) << "test" << std::endl;
    std::istringstream iss(par().mpi_split);
    
    // Vector to store the integers
    std::vector<int> mpi_split;
    
    // Temporary variable to hold each integer
    int temp;
    
    // Read integers from the string stream
    while (iss >> temp) {
        // Add the integer to the vector
        mpi_split.push_back(temp);
    }
    
    std::vector<double> mlsq(epack.eval.size());
    for(int i=0;i<epack.eval.size();i++){
        mlsq[i]=(epack.eval[i]-mass*mass) * mass;
    }
    LOG(Message) << "split_layout" << mpi_split << std::endl;
    LOG(Message) << "mpi_layout" << mpi_layout << std::endl;
    int nrhs = 1;
    for(int i=0;i<mpi_layout.size();i++){ 
        nrhs *= (mpi_layout[i]/mpi_split[i]);
        LOG(Message) << "nrhs" << nrhs << std::endl;
    }
    typedef DeflatedGuesser<FermionField>  Guesser;
    // for solution
    Guesser LMA(epack.evec, epack.eval);

    DeflatedGuesser<FermionField> LLsub(epack.evec, mlsq);
    
    typedef DeflatedGuesser<FermionField>   FineGuesser;
    //auto guesserPt = std::make_shared<DeflatedGuesser>(epack.evec, epack.eval);
    //std::shared_ptr<LinearFunction<FermionField>> guesserPt(new  FineGuesser(epack.evec, epack.eval));
    
    FermionField tmp_e(env().getRbGrid());
    FermionField tmp_o(env().getRbGrid());
    FermionField tmpRB(env().getRbGrid());
    FermionField tmpRB2(env().getRbGrid());

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
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp1);
    envGetTmp(FermionField, tmp2);
    envGetTmp(FermionField, tmp3);
    envGetTmp(FermionField, solshift);
    envGetTmp(FermionField, sourceshift);
    envGetTmp(FermionField, w);

    std::vector<FermionField>    TMP(nrhs,env().getGrid());
    
    std::vector<FermionField>    Sink(nrhs,env().getGrid());
    std::vector<FermionField>    SOL(nrhs,env().getGrid());
    for(int s=0;s<nrhs;s++){ 
        SOL[s]=Zero();
    }
    std::string outFileName;
    std::vector<std::vector<std::vector<ComplexD>>> all_results(3, 
        std::vector<std::vector<ComplexD>>(hits,
            std::vector<ComplexD>(nt, ComplexD(0., 0.))));

    
    
    GridCartesian         * UGrid   = env().getGrid();
    GridRedBlackCartesian         * RBGrid   = env().getRbGrid();

    GridCartesian * SGrid = new GridCartesian(GridDefaultLatt(),GridDefaultSimd(Nd,vComplex::Nsimd()),mpi_split,* UGrid ); 
    GridRedBlackCartesian * SrbGrid  = SpaceTimeGrid::makeFourDimRedBlackGrid(SGrid); 
    //new GridCartesian(GridDefaultLatt(),GridDefaultSimd(Nd,vComplex::Nsimd()),mpi_split,* RbGrid ); 



    
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

    int Nl_ = epack.evec.size();
    std::vector<Complex> eta(nt*3*hits*Nl_*2);
            

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
    
    int idx = 0.;
    FermionField sub(env().getGrid());
    std::vector<std::pair<int, int>> ts_mu_pairs;

    // lopp over time slice
    for(int ts=0; ts<nt;ts+=par().tinc){ //ss=split size
        
        
        LOG(Message) << "StagMesonLoopCCHLHL src_t " << ts << std::endl;
        //std::complex<double> eta = precomputedRandomNumbers[ts / par().tinc];
        
        // loop over directions
        for(int mu=0;mu<3;mu++){

            LOG(Message) << "StagMesonLoopCCHLHL src_mu " << mu << std::endl;
            ts_mu_pairs.push_back({ts, mu});

            // lopp over hits
            LOG(Message) << "Total " << hits << "hits" <<std::endl;
            for(int hit = 0; hit < hits; hit++)
            {
                // loop over evecs
                for (unsigned int il = 0; il < Nl_; il+=block)
                {
                    source = 0.0;
                    sink = 0.0;

                    //loop over blocks
                    for(int iv=il;iv<il+block;iv++){
            
                        std::complex<double> eval(mass,sqrt(epack.eval[iv]-mass*mass));
                        for(int pm=0;pm<2;pm++){
                            LOG(Message) << "Eigenvector " << 2*iv+pm << std::endl;
                            a2a.makeLowModeW(w, epack.evec[iv], eval, pm);
                            if(pm){
                                eval = conjugate(eval);
                            }
                            std::complex<double> iota_angle(0.0, std::arg(eval));
                            // source coefficient as Grid's ComplexD (std::complex on CPU, thrust::complex on GPU)
                            std::complex<double> csrc = std::exp(-iota_angle)/std::sqrt(std::abs(eval));
                            ComplexD             coefSrc(csrc.real(), csrc.imag());
                            int idx = ts * hits * 3 * Nl_ * 2 + 
                                      mu * hits * Nl_ * 2 +
                                      hit * Nl_ * 2 +
                                      iv * 2 + pm;
                
                            source += ((eta[idx])*coefSrc)*w;
                            sink += ((eta[idx])*(1./std::sqrt(std::abs(eval))))*w;
                        }
                    } 
                    
                    tmp1 = where(t == ts, source, source*0.);
                    tmp2 = adj(Umu[mu]) * tmp1;
                    // shift source at x-mu to x
                    tmp1 = Cshift(tmp2, mu, -1);
                    tmp3 = tmp1;

                    tmp1 = where(t == ts, source, source*0.);
                    // shift source
                    tmp2 = Cshift(tmp1, mu, 1);
                    tmp1 = Umu[mu] * tmp2;

                    tmp3 += tmp1;

                    sol = tmp3;
                    
                    pickCheckerboard(Even,tmp_e,sol);
                    pickCheckerboard(Odd ,tmp_o,sol);
                    
                    //once the source on full grid is made, make the preconditioned source and solution, hit the                          guesser with the solution.
                    LOG(Message) << GridLogMessage<< "make the preconditioned source"<<std::endl;
                    /////////////////////////////////////////////////////
                    // src_o = (source_o - Moe MeeInv source_e) pc
                    /////////////////////////////////////////////////////
                    action.MooeeInv(tmp_e,tmpRB);  assert( tmpRB.Checkerboard() ==Even);
                    action.Meooe   (tmpRB,tmpRB2); assert( tmpRB2.Checkerboard() ==Odd);
                    tmpRB2=tmp_o-tmpRB2;           assert( tmpRB2.Checkerboard() ==Odd);
                    action.Mooee(tmpRB2,tmp_o);

                    tmpRB2.Grid()->show_decomposition();
                    //setCheckerboard(tmp1,tmpRB2);
                    TMP[idx] = tmp3 ;
                    
                    LMA(tmp_o, tmpRB);
                    setCheckerboard(sol,tmpRB);
                    
                    // zero out even part
                    tmpRB2.Checkerboard()=Even;
                    tmpRB2=Zero();
                    setCheckerboard(sol,tmpRB2);

                    SOL[idx] = sol;
                    
                    LOG(Message) << GridLogMessage<< "Guess sol made"<<std::endl;
                    
                    Sink[idx] = sink;
                    LOG(Message) << GridLogMessage<< "nrhs sources prepared"<<std::endl;

                    
                    LatticeGaugeField s_U(SGrid);
                    FermionField s_source(SGrid);
                    FermionField s_tmp(SGrid);
                    FermionField s_sol(SGrid);

                    if ((idx+1)%nrhs ==0){ //block=16; so nrh should be 3*16 = 48

                        s_sol= 0;
                        LOG(Message) << GridLogMessage<< "split the grid"<<std::endl;
                        Grid_split  (TMP,s_tmp);
                        Grid_split  (SOL,s_sol);
                        Grid_split  (U, s_U);
                        
                        // Fermionic matrix on split grid
                        NaiveStaggeredFermion<FImpl1> Ds(s_U,*SGrid,*SrbGrid,mass,c1,tad);

                        // CG on split grid
                        ConjugateGradient<FermionField> CG(par().residual,par().maxIteration);
                        HADRONS_DEFAULT_SCHUR_SOLVE<FermionField> schurSolver(CG,par().subGuess,par().solInitGuess); 
                        schurSolver(Ds, s_tmp, s_sol);
                        

                        LOG(Message) << GridLogMessage<< "Unsplitting the result"<<std::endl;
                        Grid_unsplit(SOL,s_sol);

                        for (int s=0;s<nrhs;s++){

                           
                            int s_mu = ts_mu_pairs[s].second;
                            int s_ts = ts_mu_pairs[s].first;
                            
                            // subtract the low modes
                            sub = Zero();
                            pickCheckerboard(Even,tmp_e,TMP[s]);
                            action.Meooe(tmp_e,tmp_o);
                            LLsub(tmp_o,tmp_e);
                            action.Meooe(tmp_e,tmp_o);// tmp_o is now even
                            setCheckerboard(sub,tmp_o);
                            SOL[s] += sub;
                    
                
                            // take inner-product with eigenbra on all time slices
                            solshift = Cshift(SOL[s],s_mu,1);
                            solshift = Umu[s_mu]*solshift;
                            sliceInnerProductVector(corr,Sink[s],solshift,3); //first+second term

                        
                            for(int tsnk=0; tsnk<nt; tsnk++){
                                result[s_mu].corr[(tsnk-s_ts+nt)%nt] += (corr[tsnk]);
                            }
                    
                    
                            sourceshift = Cshift(Sink[s],s_mu,1);
                            sourceshift = Umu[s_mu]*sourceshift;
                            sliceInnerProductVector(corr,sourceshift,SOL[s],3); //third+fourth term
                    
                            // take inner-product with eigenmode on all time slices
                            for(int tsnk=0; tsnk<nt; tsnk++){
                                result[s_mu].corr[(tsnk-s_ts+nt)%nt] += (corr[tsnk]);
                            }
                           
                        }
                    
                        ts_mu_pairs.clear();
                        idx=-1;
                        
                    }
                }
            }
            LOG(Message) << GridLogMessage<< "idx = "<<idx<<std::endl;
            idx++;
        }
    }
    for (int i = 0; i < 3; ++i){
        if(U.Grid()->IsBoss()){
            makeFileDir(par().output);
            outFileName = par().output+"HLcc_2pt_mu"+std::to_string(i);
            for(int t=0; t<nt; t++)
                result[i].corr[t] = result[i].corr[t]/ComplexD(hits, 0.0);
            saveResult(outFileName, "HLCC", result[i]);
        }
    }
    
}
    
END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

#endif
// instead of 1 CG on all nodes, we want multiples CGs across all nodes
//check Grids, e-Pack l empty means no deflation thats what we want as we want to deflate before split manually so solver should be without deflation
