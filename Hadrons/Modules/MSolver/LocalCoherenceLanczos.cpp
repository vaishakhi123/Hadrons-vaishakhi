/*
 * LocalCoherenceLanczos.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#include <Hadrons/Modules/MSolver/LocalCoherenceLanczos.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<FIMPL,600>;
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<STAGIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>;
#if HADRONS_DEFAULT_LANCZOS_NBASIS != 384  // otherwise the same class as the line above
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<STAGIMPL,384>;
#endif
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>;
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<FIMPL,600, FIMPLF>;
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<ZFIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, ZFIMPLF>;
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<STAGIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, STAGIMPLF>;
#if HADRONS_DEFAULT_LANCZOS_NBASIS != 384  // otherwise the same class as the line above
template class HADRONS_NAMESPACE::MSolver::TLocalCoherenceLanczos<STAGIMPL,384, STAGIMPLF>;
#endif
#endif
