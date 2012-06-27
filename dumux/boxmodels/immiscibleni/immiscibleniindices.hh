// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \brief Defines the indices used by the 2p2cni box model
 */
#ifndef DUMUX_IMMISCIBLE_NI_INDICES_HH
#define DUMUX_IMMISCIBLE_NI_INDICES_HH

#include <dumux/boxmodels/immiscible/immiscibleindices.hh>

namespace Dumux
{
// \{

/*!
 * \ingroup ImmiscibleNIModel
 * \ingroup BoxIndices
 * \brief Enumerations for the non-isothermal 2-phase 2-component model
 *
 * \tparam formulation The formulation, either pwSn or pnSw.
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class ImmiscibleNIIndices : public ImmiscibleIndices<PVOffset>
{
    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

public:
    static const int temperatureIdx = PVOffset + numPhases; //! The index for temperature in primary variable vectors.
    static const int energyEqIdx = PVOffset + numPhases; //! The index for energy in equation vectors.
};

// \}

}
#endif
