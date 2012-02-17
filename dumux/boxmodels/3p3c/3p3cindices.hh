// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief Defines the indices required for the 3p3c BOX model.
 */
#ifndef DUMUX_3P3C_INDICES_HH
#define DUMUX_3P3C_INDICES_HH

#include "3p3cproperties.hh"

namespace Dumux
{

/*!
 * \brief The indices for the isothermal ThreePThreeC model.
 *
 * \tparam formulation The formulation, only pgSwSn
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset = 0>
class ThreePThreeCIndices
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the water phase
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the NAPL phase
    static const int gPhaseIdx = FluidSystem::gPhaseIdx; //!< Index of the gas phase

    // Component indices
    static const int wCompIdx = FluidSystem::wCompIdx; //!< Index of the water component
    static const int cCompIdx = FluidSystem::cCompIdx; //!< Index of the NAPL component
    static const int aCompIdx = FluidSystem::aCompIdx; //!< Index of the air component

    // present phases (-> 'pseudo' primary variable)
    static const int wPhaseOnly = (1 << wPhaseIdx); //!< Only the water phase is present
    static const int nPhaseOnly = (1 << nPhaseIdx); //!< Only gas phase is present
    static const int gPhaseOnly = (1 << gPhaseIdx); //!< Only gas phase is present
    static const int threePhases = wPhaseOnly | nPhaseOnly | gPhaseOnly; //!< All three phases are present
    static const int gnPhaseOnly = nPhaseOnly | gPhaseOnly; //!< Only gas and NAPL phases are present
    static const int wnPhaseOnly = wPhaseOnly | nPhaseOnly; //!< Only water and NAPL phases are present
    static const int wgPhaseOnly = wPhaseOnly | gPhaseOnly; //!< Only water and gas phases are present

    // Primary variable indices
    static const int pressure0Idx = PVOffset + 0; //!< Index pressure of the phase with the lowest index
    static const int switch1Idx = PVOffset + 1; //!< Index 1 of saturation or mole fraction
    static const int switch2Idx = PVOffset + 2; //!< Index 2 of saturation or mole fraction

    // equation indices
    static const int conti0EqIdx = PVOffset; //!< Index of the mass conservation equation for the first component
    static const int contiWEqIdx = conti0EqIdx + wCompIdx; //!< Index of the mass conservation equation for the water component
    static const int contiCEqIdx = conti0EqIdx + cCompIdx; //!< Index of the mass conservation equation for the contaminant component
    static const int contiAEqIdx = conti0EqIdx + aCompIdx; //!< Index of the mass conservation equation for the air component
};

}

#endif
