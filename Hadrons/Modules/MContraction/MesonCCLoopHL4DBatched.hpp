/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/MesonCCLoopHL4DBatched.hpp

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

#ifndef Hadrons_MContraction_MesonLoopCCHL4DBatched_hpp_
#define Hadrons_MContraction_MesonLoopCCHL4DBatched_hpp_

#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/Modules/MGuesser/ExactDeflationBLAS.hpp>
#include <Hadrons/Solver.hpp>
#include <Grid/lattice/Lattice_reduction.h>
#include <ctime>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             TStagMesonLoopCCHL4DBatched                                    *
 *                                                                            *
 * Same contraction as StagMesonLoopCCHL4D, but all sources of a time block   *
 * (tblock/tinc timeslices x 3 directions x hits) go to the solver in one     *
 * vector call, so a batched solver (deflation guess, multi-RHS HDCG) handles *
 * them together.  Results are identical to StagMesonLoopCCHL4D.              *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)


class MesonLoopCCHL4DBatchedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonLoopCCHL4DBatchedPar,
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
class TStagMesonLoopCCHL4DBatched: public Module<MesonLoopCCHL4DBatchedPar>
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
    TStagMesonLoopCCHL4DBatched(const std::string name);
    // destructor
    virtual ~TStagMesonLoopCCHL4DBatched(void) {};
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

MODULE_REGISTER_TMP(StagMesonLoopCCHL4DBatched, ARG(TStagMesonLoopCCHL4DBatched<STAGIMPL, STAGIMPL>), MContraction);

/******************************************************************************
 *                           TStagMesonLoopCCHL implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonLoopCCHL4DBatched<FImpl1, FImpl2>::TStagMesonLoopCCHL4DBatched(const std::string name)
: Module<MesonLoopCCHL4DBatchedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL4DBatched<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().gauge, par().action};
    std::string sub_string;
    
    input.push_back(par().eigenPack);
    //input.push_back(par().solver + "_subtract");
    input.push_back(par().solver);
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonLoopCCHL4DBatched<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};

    return output;
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonLoopCCHL4DBatched<FImpl1, FImpl2>::setup(void)
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
void TStagMesonLoopCCHL4DBatched<FImpl1, FImpl2>::execute(void)
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
    // low-mode subtraction |l><l|/mlsq for all sources of a time block at once (Grid BLAS deflation, whole basis)
    MGuesser::BLASDeflatedGuesser<FermionField> LLsubBLAS(epack.evec, mlsq, Nl_, 0);
    
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
    std::vector<FermionField> sol_list(tblock*3*hits,env().getGrid());


    // stage timers (seconds), reported at the end of execute()
    double tEta = 0., tLowBuild = 0., tLowW = 0., tVertex = 0., tSolve = 0., tLLsub = 0., tContr = 0.;

    tEta -= usecond();
    // Precompute and store random numbers for eta for all time slices, mu, hits and eigenvectors
    for(int ts = 0; ts < nt; ts++){
        for(int mu=0; mu< 3; mu++){
            for(int hit = 0; hit < hits; hit++){
                for(unsigned int il = 0; il < Nl_; il += block){
                    for(int iv = il; iv < std::min<int>(il + block, Nl_); iv++){
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
    
    tEta += usecond();
    FermionField sub(env().getGrid());

    // Low-mode sources and sinks of a time block as one matrix product,
    //   [source_s | sink_s](x) = sum_k A(k, s) W_k(x),   k = 2 iv + pm,
    // with W_k = makeLowModeW(evec[iv], eval, pm) and A = eta x coefficient. The W_k are built lowModeChunk at a
    // time into a contiguous buffer and accumulated with one GEMM per chunk (column-major, as MultiRHSDeflation).
    typedef typename FermionField::scalar_type   Scalar;
    typedef typename FermionField::scalar_object ScalarObject;
    const int      lowModeChunk = std::min(64, 2*Nl_);
    const int      nBlk = tblock*3*hits, ncol = 2*nBlk;
    const uint64_t vw   = env().getGrid()->lSites()*(sizeof(ScalarObject)/sizeof(Scalar));
    deviceVector<Scalar>   Wbuf(vw*lowModeChunk), Abuf(lowModeChunk*ncol), Sbuf(vw*ncol);
    deviceVector<Scalar *> Wd(1), Ad(1), Sd(1);
    std::vector<Scalar>    Ahost(lowModeChunk*ncol);
    GridBLAS               lowModeBLAS;

    acceleratorPut(Wd[0], &Wbuf[0]);
    acceleratorPut(Ad[0], &Abuf[0]);
    acceleratorPut(Sd[0], &Sbuf[0]);

    // loop over time slice
    for(int its=0; its<nt;its+=tblock){
        tLowBuild -= usecond();
        for (int k0 = 0; k0 < 2*Nl_; k0 += lowModeChunk)
        {
            const int n = std::min(lowModeChunk, 2*Nl_ - k0);

            std::fill(Ahost.begin(), Ahost.end(), Scalar(0., 0.));
            for (int kk = 0; kk < n; ++kk)
            {
                const int iv = (k0 + kk)/2, pm = (k0 + kk)%2;
                // W uses the unconjugated eval for both pm; the coefficients use conj(eval) for pm = 1
                std::complex<double> evalW(mass,sqrt(epack.eval[iv]-mass*mass));
                std::complex<double> eval = pm ? std::conj(evalW) : evalW;

                tLowW -= usecond();
                a2a.makeLowModeW(w, epack.evec[iv], evalW, pm);
                tLowW += usecond();
                {
                    autoView(wv, w, AcceleratorRead);
                    acceleratorCopyDeviceToDevice(&wv[0], &Wbuf[kk*vw], vw*sizeof(Scalar));
                }
                std::complex<double> iota_angle(0.0, std::arg(eval));
                // coefficients as Grid's ComplexD (std::complex on CPU, thrust::complex on GPU)
                std::complex<double> csrc = std::exp(-iota_angle)/std::sqrt(std::abs(eval));
                ComplexD             coefSrc(csrc.real(), csrc.imag());
                RealD                coefSnk = 1./std::sqrt(std::abs(eval));
                for (int ts=its; ts< its+tblock; ts+=par().tinc){
                    for(int mu=0;mu<3;mu++){
                        for(int hit = 0; hit < hits; hit++)
                        {
                            int s_idx = (ts-its)*3*hits + mu*hits + hit;
                            int idx = ts * hits * 3 * Nl_ * 2 + 
                                      mu * hits * Nl_ * 2 +
                                      hit * Nl_ * 2 +
                                      (k0 + kk);
                
                            Ahost[kk + n*s_idx]          = eta[idx]*coefSrc;
                            Ahost[kk + n*(nBlk + s_idx)] = eta[idx]*coefSnk;
                        }
                    }
                }
            }
            acceleratorCopyToDevice(&Ahost[0], &Abuf[0], n*ncol*sizeof(Scalar));
            // S (vw x ncol) [+]= W (vw x n) A (n x ncol); the first chunk overwrites S
            lowModeBLAS.gemmBatched(GridBLAS_OP_N, GridBLAS_OP_N, vw, ncol, n,
                                    Scalar(1.0), Wd, Ad, Scalar(k0 == 0 ? 0.0 : 1.0), Sd);
            lowModeBLAS.synchronise();
        }
        for (int b = 0; b < nBlk; ++b)
        {
            autoView(sv, source_list[b], AcceleratorWrite);
            acceleratorCopyDeviceToDevice(&Sbuf[b*vw], &sv[0], vw*sizeof(Scalar));
            autoView(kv, sink_list[b], AcceleratorWrite);
            acceleratorCopyDeviceToDevice(&Sbuf[(nBlk + b)*vw], &kv[0], vw*sizeof(Scalar));
        }
        tLowBuild += usecond();
        tVertex   -= usecond();
        // conserved-current vertex on the low-mode sources, in place: the low-mode
        // source itself is not needed again, only the vertex source
        std::vector<int> used;
        for (int ts=its; ts< its+tblock; ts+=par().tinc){
            for(int mu=0;mu<3;mu++){
                for(int hit = 0; hit < hits; hit++){
                    int s_idx = (ts-its)*3*hits + mu*hits + hit;

                    tmp = where(t == ts, source_list[s_idx], source_list[s_idx]*0.);
                    tmp2 = adj(Umu[mu]) * tmp;

                    // shift source x-mu to x
                    tmp3 = Cshift(tmp2, mu, -1);

                    // shift source
                    tmp2 = Cshift(tmp, mu, 1);
                    tmp3 += Umu[mu] * tmp2;

                    source_list[s_idx] = tmp3;
                    used.push_back(s_idx);
                }
            }
        }

        tVertex += usecond();
        tSolve  -= usecond();
        // one batched solve for all sources of this time block; the solver
        // deflates (and, with _subtract, removes the guess) for the whole batch
        LOG(Message) << "StagMesonLoopCCHL4DBatched solving " << used.size()
                     << " sources for time block " << its << std::endl;
        if (used.size() == source_list.size())
        {
            solver(sol_list, source_list);
        }
        else
        {
            std::vector<FermionField> src_batch(used.size(), env().getGrid());
            std::vector<FermionField> sol_batch(used.size(), env().getGrid());
            for (unsigned int b = 0; b < used.size(); ++b) src_batch[b] = source_list[used[b]];
            solver(sol_batch, src_batch);
            for (unsigned int b = 0; b < used.size(); ++b) sol_list[used[b]] = sol_batch[b];
        }

        tSolve += usecond();

        // low-mode subtraction for all sources of the block: sol += Meooe LL Meooe (even part of the vertex source)
        tLLsub -= usecond();
        {
            std::vector<FermionField> llIn(used.size(), env().getRbGrid()), llOut(used.size(), env().getRbGrid());

            for (unsigned int b = 0; b < used.size(); ++b)
            {
                pickCheckerboard(Even,tmp_e,source_list[used[b]]);
                action.Meooe(tmp_e,llIn[b]);
            }
            LLsubBLAS(llIn, llOut);
            for (unsigned int b = 0; b < used.size(); ++b)
            {
                action.Meooe(llOut[b],tmp_o);// tmp_o is now even
                sub = Zero();
                setCheckerboard(sub,tmp_o);
                sol_list[used[b]] += sub;
            }
        }
        tLLsub += usecond();
        tContr -= usecond();
        // loop over time slices
        for (int ts=its; ts< its+tblock; ts+=par().tinc){

            LOG(Message) << "StagMesonLoopCCHLHL src_ts " << ts << std::endl;
            // loop over directions
            for(int mu=0;mu<3;mu++){


                LOG(Message) << "StagMesonLoopCCHLHL src_mu " << mu << std::endl;

                for(int hit = 0; hit < hits; hit++){
                    int s_idx = (ts-its)*3*hits + mu*hits + hit;

                    sol = sol_list[s_idx];   // low modes already subtracted above


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
        tContr += usecond();
    }
    {
        double tot = tEta + tLowBuild + tVertex + tSolve + tLLsub + tContr;
        auto pr = [tot](const std::string &n, double t)
        {
            LOG(Message) << "StagMesonLoopCCHL4DBatched profile: " << n << " " << t/1.e6 << " s ("
                         << 100.*t/tot << "%)" << std::endl;
        };
        pr("random numbers (eta)               ", tEta);
        pr("low-mode sources: makeLowModeW     ", tLowW);
        pr("low-mode sources: GEMM + copies    ", tLowBuild - tLowW);
        pr("conserved-current vertex           ", tVertex);
        pr("solve (incl. deflation guess)      ", tSolve);
        pr("low-mode subtraction (LLsub, BLAS) ", tLLsub);
        pr("contractions                       ", tContr);
        pr("total of the above                 ", tot);
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
