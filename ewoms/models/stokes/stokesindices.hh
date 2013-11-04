// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Klaus Mosthaf                                     *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \copydoc Ewoms::StokesIndices
 */
#ifndef EWOMS_STOKES_INDICES_HH
#define EWOMS_STOKES_INDICES_HH

#include "stokesproperties.hh"
#include <ewoms/models/modules/energy/vcfvenergymodule.hh>

namespace Ewoms {
/*!
 * \ingroup VCFVStokesModel
 *
 * \brief The primary variable and equation indices of the
 *        (Navier-)Stokes model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class StokesIndices
    : public VcfvEnergyIndices<PVOffset
                                        + GET_PROP_VALUE(TypeTag, NumComponents)
                                        + GET_PROP_TYPE(TypeTag, GridView)::dimensionworld,
                                        GET_PROP_VALUE(TypeTag, EnableEnergy)>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { numComponents = FluidSystem::numComponents };


public:
    // Primary variable indices
    static const int conti0EqIdx = PVOffset + 0; //!< Index of the mass balance equation
    static const int momentum0EqIdx = PVOffset + numComponents; //!< Index of the first component of the momentum equation

    static const int pressureIdx = PVOffset + 0; //!< Index of the pressure in a solution vector
    static const int moleFrac1Idx = PVOffset + 1; //!< Index of the mole fraction of the _SECOND_ component
    static const int velocity0Idx = PVOffset + numComponents; //!< Index of the first component of the velocity
};
} // namespace Ewoms

#endif
