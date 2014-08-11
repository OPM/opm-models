/*
  Copyright (C) 2011-2014 by Andreas Lauser

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
 * \copydoc Ewoms::BlackOilPrimaryVariables
 */
#ifndef EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH
#define EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH

#include "blackoilproperties.hh"

#include <ewoms/disc/common/fvbaseprimaryvariables.hh>

#include <dune/common/fvector.hh>

#include <opm/material/constraintsolvers/ImmiscibleFlash.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 *
 * \brief Represents the primary variables used by the black-oil model.
 */
template <class TypeTag>
class BlackOilPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef FvBasePrimaryVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // primary variable indices
    enum { gasPressureIdx = Indices::gasPressureIdx };
    enum { waterSaturationIdx = Indices::waterSaturationIdx };
    enum { switchIdx = Indices::switchIdx };

    // phase indices from the fluid system
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    // component indices from the fluid system
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    static_assert(numPhases == 3,
                  "The black-oil model assumes three phases!");
    static_assert(numComponents == 3,
                  "The black-oil model assumes three components!");

public:
    BlackOilPrimaryVariables()
        : ParentType()
    { Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    BlackOilPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::SetUndefined(switchingVariableIsGasSaturation_);
        pvtRegionIdx_ = 0;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables &)
     */
    BlackOilPrimaryVariables(const BlackOilPrimaryVariables &value)
        : ParentType(value)
        , switchingVariableIsGasSaturation_(value.switchingVariableIsGasSaturation_)
        , pvtRegionIdx_(value.pvtRegionIdx_)
    { }

    void setPvtRegionIndex(int value)
    { pvtRegionIdx_ = value; }

    unsigned pvtRegionIndex() const
    { return pvtRegionIdx_; }

    bool switchingVariableIsGasSaturation() const
    { return switchingVariableIsGasSaturation_; }

    void setSwitchingVariableIsGasSaturation(bool yesno)
    { switchingVariableIsGasSaturation_ = yesno; }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        // the pressure of the first phase and the saturation of water
        // are always primary variables
        (*this)[gasPressureIdx] = fluidState.pressure(gasPhaseIdx);
        (*this)[waterSaturationIdx] = fluidState.saturation(waterPhaseIdx);

        bool gasPresent = (fluidState.saturation(gasPhaseIdx) > 0.0);
        bool oilPresent = (fluidState.saturation(oilPhaseIdx) > 0.0);
        switchingVariableIsGasSaturation_ = gasPresent || (!gasPresent && !oilPresent);

        // depending on the phases present, we use either Sg or x_o^G
        // as a primary variable
        if (switchingVariableIsGasSaturation())
            (*this)[switchIdx] = fluidState.saturation(gasPhaseIdx);
        else
            (*this)[switchIdx] = fluidState.moleFraction(oilPhaseIdx, gasPhaseIdx);
    }

    /*!
     * \brief Adapt the interpretation of the switching variable to a
     *        physically meaningful one.
     *
     * \return true Iff the interpretation of the switching variable
     *              was changed
     */
    bool adaptSwitchingVariable()
    {
        Scalar pg = (*this)[Indices::gasPressureIdx];

        if (switchingVariableIsGasSaturation()) {
            if ((*this)[Indices::waterSaturationIdx] < 1 &&
                (*this)[Indices::switchIdx] < 0.0)
            {
                // we switch to the gas mole fraction in the
                // oil phase if oil is present and if we would
                // encounter a negative gas saturation
                Scalar xoGsat = FluidSystem::saturatedOilGasMoleFraction(pg, pvtRegionIdx_);
                setSwitchingVariableIsGasSaturation(false);
                (*this)[Indices::switchIdx] = xoGsat;
                return true;
            }
        }
        else {
            // check if the amount of disolved gas in oil is
            // more that what's allowed
            Scalar xoGsat = FluidSystem::saturatedOilGasMoleFraction(pg, pvtRegionIdx_);
            if ((*this)[Indices::switchIdx] > xoGsat) {
                // yes, so we need to use gas saturation as
                // primary variable
                setSwitchingVariableIsGasSaturation(true);
                (*this)[Indices::switchIdx] = 0;
                return true;
            }
        }

        return false;
    }

private:
    bool switchingVariableIsGasSaturation_;
    unsigned char pvtRegionIdx_;
};

} // namespace Ewoms

#endif
