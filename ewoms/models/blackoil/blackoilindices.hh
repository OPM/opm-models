// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::BlackOilIndices
 */
#ifndef EWOMS_BLACK_OIL_INDICES_HH
#define EWOMS_BLACK_OIL_INDICES_HH

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the black-oil model.
 */
template <int PVOffset = 0>
struct BlackOilIndices
{
    // Primary variable indices

    //! Index of pressure of the first fluid phase in a vector of primary variables
    static const int pressure0Idx  = PVOffset + 0;
    //! Index of the saturation of the first fluid phase in a vector of primary variables
    static const int saturation0Idx = PVOffset + 1;
    // numPhases - 2 saturations follow

    // indices of the equations

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = PVOffset + 0;
    // numPhases - 1 continuity equations follow

    //! The number of equations
    static const int numEq = 3;
};

} // namespace Ewoms

#endif
