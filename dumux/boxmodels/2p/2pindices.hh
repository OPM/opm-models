// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
#ifndef DUMUX_BOX_2P_INDICES_HH
#define DUMUX_BOX_2P_INDICES_HH

namespace Dumux
{
// \{

/*!
 * \ingroup TwoPBoxModel
 * \ingroup BoxIndices
 * \brief The indices for the isothermal two-phase model.
 */
template <int PVOffset=0>
struct TwoPIndices
{
    // Phase indices
    static const int wPhaseIdx = 0; //!< Index of the wetting phase in a phase vector
    static const int nPhaseIdx = 1; //!< Index of the non-wetting phase in a phase vector

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the primary variables
    static const int pwIdx = PVOffset + 0; //!< Pressure index of the wetting phase
    static const int SnIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    // indices of the equations
    static const int conti0EqIdx = PVOffset + 0; //!< Index of the continuity equation of the first phase
};

// \}
} // namespace Dumux


#endif
