/*
 * ExactDeflationBLAS.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MGuesser_ExactDeflationBLAS_hpp_
#define Hadrons_MGuesser_ExactDeflationBLAS_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

// needs Grid's MultiRHSDeflation (BLAS deflation)
#if __has_include(<Grid/algorithms/deflation/MultiRHSDeflation.h>)
#define HADRONS_HAVE_MRHS_DEFLATION
#endif
#ifdef HADRONS_HAVE_MRHS_DEFLATION

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                  Exact deflation guesser, BLAS version                     *
 *                                                                            *
 * Same guess as ExactDeflation, guess_j = sum_i v_i <v_i, src_j> / lambda_i, *
 * but computed for all sources of a vector solve at once with Grid's         *
 * MultiRHSDeflation::DeflateSources: two GEMMs (C = E^dag R / lambda,        *
 * G = E C) and one global sum, instead of one inner product (with its own    *
 * global sum) and one axpy per (eigenvector, source) pair.                   *
 *                                                                            *
 * MultiRHSDeflation keeps its own copy of the eigenvectors it uses. With     *
 * chunkSize = 0 the whole basis is copied once and kept; with chunkSize > 0  *
 * only chunkSize eigenvectors are held at a time and the guess is            *
 * accumulated over chunks (each call then re-imports every chunk).          *
 *                                                                            *
 * Double precision only: older Grid MultiRHSDeflation does not compile for   *
 * single-precision fields.                                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGuesser)

class ExactDeflationBLASPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ExactDeflationBLASPar,
                                    std::string, eigenPack,
                                    unsigned int, size,
                                    unsigned int, chunkSize);
};

template <typename EPack>
class TExactDeflationBLAS: public Module<ExactDeflationBLASPar>
{
public:
    typedef typename EPack::Field Field;
public:
    // constructor
    TExactDeflationBLAS(const std::string name);
    // destructor
    virtual ~TExactDeflationBLAS(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ExactDeflationBLAS, TExactDeflationBLAS<BaseFermionEigenPack<FIMPL>>, MGuesser);
MODULE_REGISTER_TMP(StagExactDeflationBLAS, TExactDeflationBLAS<BaseFermionEigenPack<STAGIMPL>>, MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template<class Field>
class BLASDeflatedGuesser: public LinearFunction<Field> {
private:
    std::vector<Field>          &evec;
    std::vector<RealD>          &eval;
    const unsigned int          epackSize, chunkSize;
    MultiRHSDeflation<Field>    defl;
    int                         imported;  // first eigenvector of the chunk in defl, -1 if none

    void importChunk(const unsigned int ev0, const unsigned int n)
    {
        if (imported != (int)ev0 || (unsigned int)defl.nev != n)
        {
            defl.ImportEigenBasis(evec, eval, ev0, n);
            imported = ev0;
        }
    }
public:
    using LinearFunction<Field>::operator();

    BLASDeflatedGuesser(std::vector<Field>& evec_, std::vector<RealD>& eval_,
                        const unsigned int epackSize_, const unsigned int chunkSize_)
    : evec(evec_), eval(eval_), epackSize(epackSize_)
    , chunkSize((chunkSize_ == 0 || chunkSize_ > epackSize_) ? epackSize_ : chunkSize_)
    , imported(-1)
    {
        assert(evec.size()==eval.size());
        assert(epackSize <= evec.size());
    }

    virtual void operator()(const Field &src,Field &guess) {
        std::vector<Field> srcVec   = {src};
        std::vector<Field> guessVec = {guess};

        (*this)(srcVec,guessVec);

        guess = guessVec[0];
    }

    virtual void operator() (const std::vector<Field> &src, std::vector<Field> &guess)
    {
        assert(src.size() == guess.size());

        unsigned int       sourceSize = src.size();
        // DeflateSources takes non-const sources
        std::vector<Field> srcCopy(src);
        std::vector<Field> part(sourceSize, src[0].Grid());
        double             time = 0.;

        // DeflateSources requires guesses on the same checkerboard as the sources
        for (unsigned int k = 0; k < sourceSize; ++k)
        {
            part[k].Checkerboard()  = src[k].Checkerboard();
            guess[k]                = Zero();
            guess[k].Checkerboard() = src[k].Checkerboard();
        }
        time -= usecond();
        for (unsigned int ev0 = 0; ev0 < epackSize; ev0 += chunkSize)
        {
            unsigned int n = std::min(chunkSize, epackSize - ev0);

            importChunk(ev0, n);
            defl.DeflateSources(srcCopy, part);
            for (unsigned int k = 0; k < sourceSize; ++k)
            {
                guess[k] += part[k];
            }
        }
        time += usecond();
        LOG(Message) << "BLASDeflatedGuesser: " << sourceSize << " sources, "
                     << epackSize << " eigenvectors (chunks of " << chunkSize
                     << "), total time " << time/1.e6 << " s" << std::endl;
    }
};

/******************************************************************************
 *                 TExactDeflationBLAS implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename EPack>
TExactDeflationBLAS<EPack>::TExactDeflationBLAS(const std::string name)
: Module<ExactDeflationBLASPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename EPack>
std::vector<std::string> TExactDeflationBLAS<EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};

    return in;
}

template <typename EPack>
std::vector<std::string> TExactDeflationBLAS<EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename EPack>
DependencyMap TExactDeflationBLAS<EPack>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename EPack>
void TExactDeflationBLAS<EPack>::setup(void)
{
    LOG(Message) << "Setting BLAS exact deflation guesser with eigenpack '"
                 << par().eigenPack << "' ("
                 << par().size << " modes, chunk size "
                 << par().chunkSize << ")" << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>, BLASDeflatedGuesser<Field>, getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.eval,
                     par().size, par().chunkSize);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack>
void TExactDeflationBLAS<EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // HADRONS_HAVE_MRHS_DEFLATION

#endif // Hadrons_MGuesser_ExactDeflationBLAS_hpp_
