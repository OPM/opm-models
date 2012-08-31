// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Melanie Darcis                                    *
 *   Copyright (C) 2011 by Bernd Flemisch                                    *
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
 * \brief Defines the indices required for the two-phase box model.
 */
#ifndef DUMUX_BOX_IMMISCIBLE_INDICES_HH
#define DUMUX_BOX_IMMISCIBLE_INDICES_HH

#include "immiscibleproperties.hh"
#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>

namespace Dumux
{
/*!
 * \ingroup ImmiscibleBoxModel
 * \ingroup BoxIndices
 * \brief The indices for the isothermal two-phase model.
 */
template <class TypeTag, int PVOffset>
struct ImmiscibleIndices
    : public BoxMultiPhaseEnergyIndices<PVOffset + GET_PROP_VALUE(TypeTag, NumPhases), GET_PROP_VALUE(TypeTag, EnableEnergy)>

{
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef BoxMultiPhaseEnergyIndices<PVOffset + numPhases, enableEnergy> EnergyIndices;

public:
    // number of equations/primary variables
    static const int numEq = numPhases + EnergyIndices::numEq_;

    // Primary variable indices
    static const int pressure0Idx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturation0Idx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the equations
    static const int conti0EqIdx = PVOffset + 0; //!< Index of the continuity equation of the first phase
};
} // namespace Dumux


#endif
