// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::ImmiscibleIndices
 */
#ifndef EWOMS_IMMISCIBLE_INDICES_HH
#define EWOMS_IMMISCIBLE_INDICES_HH

#include "immiscibleproperties.hh"
#include <ewoms/models/modules/energy/vcfvenergymodule.hh>

namespace Ewoms {

/*!
 * \ingroup ImmiscibleVcfvModel
 *
 * \brief The indices for the isothermal multi-phase model.
 */
template <class TypeTag, int PVOffset>
struct ImmiscibleIndices
    : public VcfvEnergyIndices<PVOffset + GET_PROP_VALUE(TypeTag, NumPhases), GET_PROP_VALUE(TypeTag, EnableEnergy)>

{
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef VcfvEnergyIndices<PVOffset + numPhases, enableEnergy> EnergyIndices;

public:
    // number of equations/primary variables
    static const int numEq = numPhases + EnergyIndices::numEq_;

    // Primary variable indices
    static const int pressure0Idx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturation0Idx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the equations
    static const int conti0EqIdx = PVOffset + 0; //!< Index of the continuity equation of the first phase
};
} // namespace Ewoms


#endif
