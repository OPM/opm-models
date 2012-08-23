// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011 by Benjamin Faigle                                   *
 *   Copyright (C) 2009-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2010 by Melanie Darcis                                    *
 *   Copyright (C) 2010 by Philipp Nuske                                     *
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
 * \brief Defines the indices required for the black-oil box model.
 */
#ifndef DUMUX_BOX_BLACK_OIL_INDICES_HH
#define DUMUX_BOX_BLACK_OIL_INDICES_HH

namespace Dumux {
/*!
 * \ingroup BlackOilBoxModel
 * \ingroup BoxIndices
 * \brief The indices for the isothermal black-oil model.
 */
template <int PVOffset=0>
struct BlackOilIndices
{
    // Primary variable indices
    static const int pressure0Idx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturation0Idx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the equations
    static const int conti0EqIdx = PVOffset + 0; //!< Index of the continuity equation of the first phase
};

} // namespace Dumux

#endif
