// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Dumux::BlackOilPrimaryVariables
 */
#ifndef DUMUX_BLACK_OIL_PRIMARY_VARIABLES_HH
#define DUMUX_BLACK_OIL_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/material/constraintsolvers/immiscibleflash.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

#include "blackoilproperties.hh"

namespace Dumux {

/*!
 * \ingroup BlackOilBoxModel
 *
 * \brief Represents the primary variables used by the black-oil model.
 */
template <class TypeTag>
class BlackOilPrimaryVariables
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    // phase indices from the fluid system
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { oPhaseIdx = FluidSystem::oPhaseIdx };

    // component indices from the fluid system
    enum { gCompIdx = FluidSystem::gCompIdx };
    enum { wCompIdx = FluidSystem::wCompIdx };
    enum { oCompIdx = FluidSystem::oCompIdx };

    static_assert(numPhases == 3, "The black-oil model has three phases!");
    static_assert(numComponents == 3, "The black-oil model has three components!");

public:
    BlackOilPrimaryVariables()
        : ParentType()
    { Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    BlackOilPrimaryVariables(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables &)
     */
    BlackOilPrimaryVariables(const BlackOilPrimaryVariables &value)
        : ParentType(value)
    { }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        (*this)[pressure0Idx] = fluidState.pressure(/*phaseIdx=*/0);

        Scalar saturation[numPhases];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            saturation[phaseIdx] = fluidState.saturation(phaseIdx);

        // check whether the oil is undersaturated
        Scalar po = fluidState.pressure(oPhaseIdx);
        Scalar poSat = FluidSystem::oilSaturationPressure(fluidState.massFraction(oPhaseIdx, gCompIdx));
        if (po > poSat) {
            // use a "negative saturation" of the gas phase to
            // compensate

            // first we calculate the "mass per cubic meter of pore
            // space" of gas which is missing to make the oil
            // saturated
            Scalar Bo = FluidSystem::oilFormationVolumeFactor(po);
            Scalar rhoo = FluidSystem::surfaceDensity(oPhaseIdx)/Bo;
            Scalar rhorefg = FluidSystem::surfaceDensity(gPhaseIdx);
            Scalar RsSat = FluidSystem::gasFormationFactor(po);
            Scalar XoGSat = RsSat*rhorefg / rhoo;

            Scalar rhogDef =
                fluidState.saturation(oPhaseIdx)
                * rhoo
                * (XoGSat - fluidState.massFraction(oPhaseIdx, gCompIdx));

            Scalar Bg = FluidSystem::gasFormationVolumeFactor(po);
            Scalar rhog = rhorefg/Bg;

            saturation[gPhaseIdx] = - rhogDef/rhog;
            saturation[oPhaseIdx] = 1 - saturation[wPhaseIdx] - saturation[gPhaseIdx];
        }

        (*this)[saturation0Idx] = saturation[/*phaseIdx=*/0];
        (*this)[saturation0Idx + 1] = saturation[/*phaseIdx=*/1];
    }
};

} // end namepace

#endif
