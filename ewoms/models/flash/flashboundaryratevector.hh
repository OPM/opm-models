// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Ewoms::FlashBoundaryRateVector
 */
#ifndef EWOMS_FLASH_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_FLASH_BOUNDARY_RATE_VECTOR_HH

#include "flashproperties.hh"

#include <ewoms/models/modules/energy/vcfvenergymodule.hh>
#include <opm/material/Valgrind.hpp>

namespace Ewoms {

/*!
 * \ingroup FlashModel
 *
 * \brief Implements a boundary vector for the fully implicit
 *        compositional multi-phase model which is based on flash
 *        calculations.
 */
template <class TypeTag>
class FlashBoundaryRateVector
    : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef VcfvEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    FlashBoundaryRateVector()
        : ParentType()
    { }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    FlashBoundaryRateVector(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector &)
     */
    FlashBoundaryRateVector(const FlashBoundaryRateVector &value)
        : ParentType(value)
    { }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setFreeFlow
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context &context,
                     int bfIdx,
                     int timeIdx,
                     const FluidState &fluidState)
    {
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        FluxVariables fluxVars;
        fluxVars.updateBoundary(context, bfIdx, timeIdx, fluidState, paramCache);
        const auto &insideVolVars = context.volVars(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar meanMBoundary = 0;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                meanMBoundary += fluidState.moleFraction(phaseIdx, compIdx)*FluidSystem::molarMass(compIdx);

            Scalar density;
            if (fluidState.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
                density = FluidSystem::density(fluidState, paramCache, phaseIdx);
            else
                density = insideVolVars.fluidState().density(phaseIdx);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                Scalar molarity;
                if (fluidState.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx)) {
                    molarity =
                        fluidState.moleFraction(phaseIdx, compIdx)
                        * density / meanMBoundary;
                }
                else  {
                    molarity = insideVolVars.fluidState().molarity(phaseIdx, compIdx);
                }


                // add advective flux of current component in current
                // phase
                (*this)[conti0EqIdx + compIdx] +=
                    fluxVars.volumeFlux(phaseIdx)
                    * molarity;
            }

            if (enableEnergy) {
                Scalar specificEnthalpy;
                if (fluidState.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
                    specificEnthalpy = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
                else
                    specificEnthalpy = insideVolVars.fluidState().enthalpy(phaseIdx);

                // currently we neglect heat conduction!
                Scalar enthalpyRate =
                    density
                    * fluxVars.volumeFlux(phaseIdx)
                    * specificEnthalpy;
                EnergyModule::addToEnthalpyRate(*this, enthalpyRate);
            }
        }

        EnergyModule::addToEnthalpyRate(*this, EnergyModule::heatConductionRate(fluxVars));

#ifndef NDEBUG
        for (int i = 0; i < numEq; ++ i) {
            Valgrind::CheckDefined((*this)[i]);
        };
#endif
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setInFlow
     */
    template <class Context, class FluidState>
    void setInFlow(const Context &context,
                   int bfIdx,
                   int timeIdx,
                   const FluidState &fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the direction opposite to the outer
        // unit normal
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            Scalar &val = this->operator[](eqIdx);
            val = std::min<Scalar>(0.0, val);
        };
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setOutFlow
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context &context,
                    int bfIdx,
                    int timeIdx,
                    const FluidState &fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the same direction as the outer
        // unit normal
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            Scalar &val = this->operator[](eqIdx);
            val = std::max<Scalar>(0.0, val);
        };
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setNoFlow
     */
    void setNoFlow()
    { (*this) = 0.0; }
};

} // namespace Ewoms

#endif
