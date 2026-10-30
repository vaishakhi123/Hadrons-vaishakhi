/*
 * StagSparseA2AVectorsGridIo.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Vaishakhi Moningi <vaishu.moningi@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MSolver_StagSparseA2AVectorsGridIo_hpp_
#define Hadrons_MSolver_StagSparseA2AVectorsGridIo_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 * Sparsened staggered A2A vectors streaming eigenvectors from Grid LIME files *
 * Single-file mode: opens the file once; uses skipScidacFieldRecord() to     *
 * seek past evecStart records in O(1) I/O — no data reads during the skip.  *
 * Multi-file mode: opens each eigenvector file directly.                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class StagSparseA2AVectorsGridIoPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagSparseA2AVectorsGridIoPar,
                                    std::string, action,
                                    std::string, gauge,
                                    std::string, evecPath,
                                    std::string, output,
                                    int,         numEvecs,
                                    int,         evecStart,
                                    int,         inc,
                                    int,         tinc,
                                    double,      mass,
                                    bool,        multiFile,
                                    bool,        milcEvecs);
};

template <typename FImpl>
class TStagSparseA2AVectorsGridIo : public Module<StagSparseA2AVectorsGridIoPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsLowStaggered<FImpl> A2A;
    typedef typename Grid::NaiveStaggeredFermionD::FermionField SparseFermionField;

public:
    TStagSparseA2AVectorsGridIo(const std::string name);
    virtual ~TStagSparseA2AVectorsGridIo(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
private:
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(StagSparseA2AVectorsGridIo,
                    ARG(TStagSparseA2AVectorsGridIo<STAGIMPL>), MSolver);

/******************************************************************************
 *                  TStagSparseA2AVectorsGridIo implementation                *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagSparseA2AVectorsGridIo<FImpl>::TStagSparseA2AVectorsGridIo(const std::string name)
: Module<StagSparseA2AVectorsGridIoPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagSparseA2AVectorsGridIo<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    in.push_back(par().gauge);
    in.push_back(par().action);
    return in;
}

template <typename FImpl>
std::vector<std::string> TStagSparseA2AVectorsGridIo<FImpl>::getOutput(void)
{
    return {getName()+"_v", getName()+"_w0", getName()+"_w1", getName()+"_w2"};
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AVectorsGridIo<FImpl>::setup(void)
{
    auto &action = envGet(FMat, par().action);

    Nl_ = par().numEvecs;
    LOG(Message) << "StagSparseA2AVectorsGridIo setup: numEvecs=" << par().numEvecs
                 << " evecStart=" << par().evecStart
                 << " inc=" << par().inc
                 << " tinc=" << par().tinc << std::endl;

    envTmp(A2A, "a2a", 1, action);
    envTmp(FermionField, "tempEvec", 1, envGetRbGrid(FermionField));

    int bsinc  = par().inc  > 0 ? par().inc  : 1;
    int bstinc = par().tinc > 0 ? par().tinc : 1;
    GridBase *sgrid;
    if (bsinc > 1 || bstinc > 1)
    {
        std::vector<int> blocksize = {bsinc, bsinc, bsinc, bstinc};
        sgrid = envGetCoarseGrid(SparseFermionField, blocksize);
    }
    else
    {
        sgrid = envGetGrid(SparseFermionField);
    }
    envCreate(std::vector<SparseFermionField>, getName() + "_v",  1, 2*Nl_, sgrid);
    envCreate(std::vector<SparseFermionField>, getName() + "_w0", 1, 2*Nl_, sgrid);
    envCreate(std::vector<SparseFermionField>, getName() + "_w1", 1, 2*Nl_, sgrid);
    envCreate(std::vector<SparseFermionField>, getName() + "_w2", 1, 2*Nl_, sgrid);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagSparseA2AVectorsGridIo<FImpl>::execute(void)
{
    auto    &action = envGet(FMat, par().action);
    auto    &U      = envGet(LatticeGaugeField, par().gauge);
    double   mass   = par().mass;
    uint64_t nt     = env().getDim(Tp);
    uint64_t ns     = env().getDim(Xp);
    envGetTmp(A2A, a2a);

    auto &v  = envGet(std::vector<SparseFermionField>, getName() + "_v");
    auto &w0 = envGet(std::vector<SparseFermionField>, getName() + "_w0");
    auto &w1 = envGet(std::vector<SparseFermionField>, getName() + "_w1");
    auto &w2 = envGet(std::vector<SparseFermionField>, getName() + "_w2");
    const int traj = vm().getTrajectory();

    envGetTmp(FermionField, tempEvec);
    tempEvec.Checkerboard() = Odd;
    assert(tempEvec.Checkerboard() == Odd);
    LOG(Message) << "tempEvec grid (RbGrid):" << std::endl;
    tempEvec.Grid()->show_decomposition();
    LOG(Message) << "tempEvec checker_dim=" << tempEvec.Grid()->_checker_dim << std::endl;

    RealD    currentEval = 0.;
    PackRecord packRecord;
    ScidacReader binReader;

    std::string t       = "." + std::to_string(traj);
    std::string stem    = par().evecPath + t;           // directory for multiFile
    std::string binFile = par().evecPath + t + ".bin";  // single-file path

    LOG(Message) << "Computing sparse A2A vectors streaming " << 2*Nl_
                 << " low modes from " << (par().multiFile ? stem : binFile) << std::endl;

    LOG(Message) << " Full grid: " << std::endl;
    U.Grid()->show_decomposition();
    LOG(Message) << " Sparse grid: " << std::endl;
    v[0].Grid()->show_decomposition();

    // Staggered phases (spatial gamma only)
    Lattice<iScalar<vInteger>> x(U.Grid());     LatticeCoordinate(x, 0);
    Lattice<iScalar<vInteger>> y(U.Grid());     LatticeCoordinate(y, 1);
    Lattice<iScalar<vInteger>> lin_z(U.Grid()); lin_z = x + y;

    ComplexField phases(U.Grid());
    FermionField temp(U.Grid());
    FermionField temp2(U.Grid());

    int step = 2 * par().inc;
    std::uniform_int_distribution<uint32_t> uid(0, step - 1);
    std::vector<uint32_t> xshift(nt), yshift(nt), zshift(nt);
    if (par().inc != 1)
    {
        for (int tt = 0; tt < (int)nt; tt++)
        {
            xshift[tt] = uid(rngSerial()._generators[0]);
            yshift[tt] = uid(rngSerial()._generators[0]);
            zshift[tt] = uid(rngSerial()._generators[0]);
        }
    }
    else
    {
        for (int tt = 0; tt < (int)nt; tt++)
        { xshift[tt] = 0; yshift[tt] = 0; zshift[tt] = 0; }
    }
    CartesianCommunicator::BroadcastWorld(0, (void *)&xshift[0], sizeof(uint32_t)*xshift.size());
    CartesianCommunicator::BroadcastWorld(0, (void *)&yshift[0], sizeof(uint32_t)*yshift.size());
    CartesianCommunicator::BroadcastWorld(0, (void *)&zshift[0], sizeof(uint32_t)*zshift.size());
    LOG(Message) << "xshift" << xshift << std::endl;
    LOG(Message) << "yshift" << yshift << std::endl;
    LOG(Message) << "zshift" << zshift << std::endl;

    std::vector<ComplexD> evalM(2*Nl_);   // ComplexD: saveEvalBlock takes Grid ComplexD (thrust::complex on GPU)

    int locx    = U.Grid()->_ldimensions[0];
    int locy    = U.Grid()->_ldimensions[1];
    int locz    = U.Grid()->_ldimensions[2];
    int loct    = U.Grid()->_ldimensions[3];
    int lstartx = U.Grid()->_lstart[0];
    int lstarty = U.Grid()->_lstart[1];
    int lstartz = U.Grid()->_lstart[2];
    int lstartt = U.Grid()->_lstart[3];

    // Detect convention from evec 0 regardless of evecStart.
    // For evecStart>0 single-file: open, read header+evec0, close — then reopen
    // for the main loop.  For multiFile: open v0.bin directly.
    bool masslessDdagD = false;
    {
        RealD eval0 = 0.;
        if (par().evecStart == 0)
        {
            // Convention detected inside the main loop at il==0.
        }
        else if (par().multiFile)
        {
            std::string fname0 = stem + "/v0.bin";
            ScidacReader r0;
            FermionField evec0(tempEvec.Grid());
            PackRecord   pr0;
            r0.open(fname0);
            EigenPackIo::readHeader(pr0, r0);
            EigenPackIo::readElement(evec0, eval0, 0, r0);
            r0.close();
            masslessDdagD = (eval0 < mass * mass);
            LOG(Message) << "Eigenpack convention (from v0.bin): "
                         << (masslessDdagD ? "massless DdagD" : "massive (D+m)dag(D+m)")
                         << " (eval0=" << eval0 << ", m^2=" << mass*mass << ")" << std::endl;
        }
        else
        {
            // Single-file: open, read header, read evec 0, close.
            ScidacReader r0;
            FermionField evec0(tempEvec.Grid());
            PackRecord   pr0;
            r0.open(binFile);
            EigenPackIo::readHeader(pr0, r0);
            EigenPackIo::readElement(evec0, eval0, 0, r0);
            r0.close();
            masslessDdagD = (eval0 < mass * mass);
            LOG(Message) << "Eigenpack convention (from evec 0): "
                         << (masslessDdagD ? "massless DdagD" : "massive (D+m)dag(D+m)")
                         << " (eval0=" << eval0 << ", m^2=" << mass*mass << ")" << std::endl;
        }
    }

    // Open single file once; use skipScidacFieldRecord() to seek past evecStart
    // records — this reads only LIME headers (~bytes) and calls fseek() past the
    // binary data, so it is near-instantaneous even for large evecStart.
    if (!par().multiFile)
    {
        binReader.open(binFile);
        EigenPackIo::readHeader(packRecord, binReader);
        if (par().evecStart > 0)
        {
            LOG(Message) << "Fast-seeking past " << par().evecStart
                         << " eigenvectors (LIME header seek, no data read)" << std::endl;
            for (int sk = 0; sk < par().evecStart; sk++)
                binReader.skipScidacFieldRecord();
            LOG(Message) << "Seek complete." << std::endl;
        }
    }

    for (unsigned int il = 0; il < 2*Nl_; il++)
    {
        if (il % 2 == 0)
        {
            int k    = il / 2;
            int kabs = k + par().evecStart;
            startTimer("evec read");
            LOG(Message) << "Reading eigenvector " << kabs << std::endl;
            if (par().multiFile)
            {
                std::string fname = stem + "/v" + std::to_string(kabs) + ".bin";
                binReader.open(fname);
                EigenPackIo::readHeader(packRecord, binReader);
                EigenPackIo::readElement(tempEvec, currentEval, kabs, binReader);
                binReader.close();
            }
            else
            {
                EigenPackIo::readElement(tempEvec, currentEval, kabs, binReader);
            }
            stopTimer("evec read");
            if (il == 0)
            {
                LOG(Message) << "tempEvec grid dimensions: " << tempEvec.Grid()->GlobalDimensions() << std::endl;
                LOG(Message) << "norm2(tempEvec)=          " << norm2(tempEvec) << std::endl;
                LOG(Message) << "tempEvec checkerboard after readElement: " << tempEvec.Checkerboard() << std::endl;
            }
        }

        if (il == 0 && par().evecStart == 0)
        {
            masslessDdagD = (currentEval < mass * mass);
            LOG(Message) << "Eigenpack convention: "
                         << (masslessDdagD ? "massless DdagD" : "massive (D+m)dag(D+m)")
                         << " (eval0=" << currentEval << ", m^2=" << mass*mass << ")" << std::endl;
        }

        double lambda = masslessDdagD ? sqrt(currentEval)
                                      : sqrt(currentEval - mass * mass);
        std::complex<double> eval(mass, lambda);
        std::complex<double> eval_for_W = par().milcEvecs
                                          ? std::complex<double>(mass, lambda / 2.0)
                                          : eval;

        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low modes)" << std::endl;
        a2a.makeLowModeW(temp, tempEvec, eval_for_W, il % 2);
        if (il < 2) LOG(Message) << "norm2(W low mode, il=" << il << ") = " << norm2(temp) << std::endl;
        stopTimer("W low mode");

        il % 2 ? eval = conjugate(eval) : eval;
        evalM[il] = ComplexD(eval.real(), eval.imag());

        v[il]  = Zero();
        w0[il] = Zero();
        w1[il] = Zero();
        w2[il] = Zero();

        for (int mu = 0; mu < 3; mu++)
        {
            phases = 1.0;
            if (mu == 1)
            {
                phases = where(mod(x,     2) == (Integer)0, phases, -phases);
            }
            else if (mu == 2)
            {
                phases = where(mod(lin_z, 2) == (Integer)0, phases, -phases);
            }
            LatticeColourMatrix Umu(U.Grid());
            Umu  = PeekIndex<LorentzIndex>(U, mu);
            Umu *= phases;

            temp2 = Umu * Cshift(temp, mu, 1);

            thread_for(tt, loct, {
                int tglb = tt + lstartt;
                if (tt % 2 == 1) continue;

                Coordinate site(Nd);
                Coordinate sparseSite(Nd);
                ColourVector vec;

                for (int z = 0; z < (int)ns; z += step) {
                    int zg = (zshift[tglb] + z) % ns;
                    for (int zl = 0; zl < locz; zl++) {
                        int zgp = zl + lstartz;
                        if (zgp == zg || zgp == (zg+1) % ns) {
                            site[2] = zl;
                            if (par().inc == 1) {
                                sparseSite[2] = site[2];
                            } else if (zshift[tglb] != step-1) {
                                sparseSite[2] = 2*int(site[2]/step) + (site[2]+zshift[tglb])%2;
                            } else {
                                sparseSite[2] = 2*int(site[2]/step) + (site[2])%2;
                            }
                            for (int y = 0; y < (int)ns; y += step) {
                                int yg = (yshift[tglb] + y) % ns;
                                for (int yl = 0; yl < locy; yl++) {
                                    int ygp = yl + lstarty;
                                    if (ygp == yg || ygp == (yg+1) % ns) {
                                        site[1] = yl;
                                        if (par().inc == 1) {
                                            sparseSite[1] = site[1];
                                        } else if (yshift[tglb] != step-1) {
                                            sparseSite[1] = 2*int(site[1]/step) + (site[1]+yshift[tglb])%2;
                                        } else {
                                            sparseSite[1] = 2*int(site[1]/step) + (site[1])%2;
                                        }
                                        for (int xx = 0; xx < (int)ns; xx += step) {
                                            int xg = (xshift[tglb] + xx) % ns;
                                            for (int xl = 0; xl < locx; xl++) {
                                                int xgp = xl + lstartx;
                                                if (xgp == xg || xgp == (xg+1) % ns) {
                                                    site[0] = xl;
                                                    if (par().inc == 1) {
                                                        sparseSite[0] = site[0];
                                                    } else if (xshift[tglb] != step-1) {
                                                        sparseSite[0] = 2*int(site[0]/step) + (site[0]+xshift[tglb])%2;
                                                    } else {
                                                        sparseSite[0] = 2*int(site[0]/step) + (site[0])%2;
                                                    }
                                                    for (int that = 0; that < 2; that++) {
                                                        site[3]       = tt + that;
                                                        sparseSite[3] = site[3];
                                                        if (mu == 0) {
                                                            peekLocalSite(vec, temp,  site);
                                                            pokeLocalSite(vec, v[il],  sparseSite);
                                                            peekLocalSite(vec, temp2, site);
                                                            pokeLocalSite(vec, w0[il], sparseSite);
                                                        } else if (mu == 1) {
                                                            peekLocalSite(vec, temp2, site);
                                                            pokeLocalSite(vec, w1[il], sparseSite);
                                                        } else if (mu == 2) {
                                                            peekLocalSite(vec, temp2, site);
                                                            pokeLocalSite(vec, w2[il], sparseSite);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }); // thread_for
        } // mu
    } // il

    if (!par().multiFile)
        binReader.close();

    std::string dir = dirname(par().output);
    if (!par().output.empty())
    {
        int status = mkdir(dir);
        if (status)
        {
            HADRONS_ERROR(Io, "cannot create directory '" + dir
                          + "' (" + std::strerror(errno) + ")");
        }
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v",  v,  par().multiFile, traj);
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w0", w0, par().multiFile, traj);
        A2AVectorsIo::write(par().output + "_w1", w1, par().multiFile, traj);
        A2AVectorsIo::write(par().output + "_w2", w2, par().multiFile, traj);
        stopTimer("W I/O");
    }

    if (env().getGrid()->IsBoss())
    {
        std::string eval_filename;
        if (!par().output.empty())
            eval_filename = A2AVectorsIo::evalFilename(par().output, traj);
        else
            eval_filename = A2AVectorsIo::evalFilename("evals", traj);
        A2AVectorsIo::initEvalFile(eval_filename, evalM.size());
        A2AVectorsIo::saveEvalBlock(eval_filename, evalM.data(), 0, 2*Nl_);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_StagSparseA2AVectorsGridIo_hpp_
