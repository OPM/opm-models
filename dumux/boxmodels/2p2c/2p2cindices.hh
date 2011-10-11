// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
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
 * \brief Defines the indices required for the 2p2c BOX model.
 */
#ifndef DUMUX_2P2C_INDICES_HH
#define DUMUX_2P2C_INDICES_HH

#include "2p2cproperties.hh"

namespace Dumux
{
// \{

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxIndices
 * \brief Enumerates the formulations which the 2p2c model accepts.
 */
struct TwoPTwoCFormulation
{
    enum {
        plSg,
        pgSl
    };
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model.
 *
 * \tparam formulation The formulation, either pwSn or pnSw.
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag,
          int formulation = TwoPTwoCFormulation::plSg,
          int PVOffset = 0>
class TwoPTwoCIndices
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    static_assert(FluidSystem::numPhases == 2,
                  "Number of phases provided by the fluid system must be 2!");
    static_assert(FluidSystem::numComponents == 2,
                  "Number of components provided by the fluid system must be 2!");

public:
    // present phases (-> 'pseudo' primary variable)
    static const int lPhaseOnly = -10; //!< Only the non-wetting phase is present
    static const int gPhaseOnly = -20; //!< Only the wetting phase is present
    static const int bothPhases = -30; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pwIdx = pressureIdx; //!< Index for liquid phase pressure in a solution vector
    static const int SnOrXIdx = switchIdx; //!< Index of the either the saturation of the non-wetting phase or the mass fraction solute in the only phase

    // equation indices
    static const int conti0EqIdx = PVOffset; //!< Index of the mass conservation equation for the first component

    static const int DUMUX_DEPRECATED_MSG("use 0 instead") lCompIdx = 0;
    static const int DUMUX_DEPRECATED_MSG("use 1 instead") gCompIdx = 1;
    static const int DUMUX_DEPRECATED_MSG("use 0 instead") lPhaseIdx = 0;
    static const int DUMUX_DEPRECATED_MSG("use 1 instead") gPhaseIdx = 1;
    static const int DUMUX_DEPRECATED_MSG("use conti0EqIdx + componentIdx instead") contiLEqIdx = conti0EqIdx + 0;
    static const int DUMUX_DEPRECATED_MSG("use conti0EqIdx + componentIdx instead") contiGEqIdx = conti0EqIdx + 1;
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model in the pg-Sl
 *        formulation.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class TwoPTwoCIndices<TypeTag, TwoPTwoCFormulation::pgSl, PVOffset>
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    static_assert(FluidSystem::numPhases == 2,
                  "Number of phases provided by the fluid system must be 2!");
    static_assert(FluidSystem::numComponents == 2,
                  "Number of components provided by the fluid system must be 2!");
public:
    // present phases (-> 'pseudo' primary variable)
    static const int lPhaseOnly = -10; //!< Only the non-wetting phase is present
    static const int gPhaseOnly = -20; //!< Only the wetting phase is present
    static const int bothPhases = -30; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pnIdx = pressureIdx; //!< Index for gas phase pressure in a solution vector
    static const int SwOrXIdx = switchIdx; //!< Index of the either the saturation of the wetting phase or the mass fraction solute in the only phase

    // Equation indices
    static const int conti0EqIdx = PVOffset + 0; //!< Index of the mass conservation equation of the first component
};

// \}

}

#endif
