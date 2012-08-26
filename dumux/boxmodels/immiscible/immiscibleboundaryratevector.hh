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
 * \brief Implements a boundary vector for the fully implicit two-phase model.
 */
#ifndef DUMUX_BOX_IMMISCIBLE_BOUNDARY_RATE_VECTOR_HH
#define DUMUX_BOX_IMMISCIBLE_BOUNDARY_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/material/constraintsolvers/ncpflash.hh>

#include "immisciblevolumevariables.hh"

namespace Dumux
{
/*!
 * \ingroup ImmiscibleModel
 *
 * \brief Implements a boundary vector for the fully implicit two-phase model.
 */
template <class TypeTag>
class ImmiscibleBoundaryRateVector
    : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef BoxMultiPhaseEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Default constructor
     */
    ImmiscibleBoundaryRateVector()
        : ParentType()
    { }

    /*!
     * \brief Constructor with assignment from scalar
     */
    ImmiscibleBoundaryRateVector(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \brief Copy constructor
     */
    ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector &value)
        : ParentType(value)
    { }

    /*!
     * \brief Specify a free-flow boundary
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context &context, 
                     int bfIdx, 
                     int timeIdx,
                     const FluidState &fs)
    {
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fs);

        FluxVariables fluxVars;
        fluxVars.updateBoundary(context, bfIdx, timeIdx, fs, paramCache);      
        const auto &insideVolVars = context.volVars(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            Scalar density;
            if (fs.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
                density = FluidSystem::density(fs, paramCache, phaseIdx);
            else 
                density = insideVolVars.fluidState().density(phaseIdx);

            // add advective flux of current component in current
            // phase
            (*this)[conti0EqIdx + phaseIdx] +=
                fluxVars.volumeFlux(phaseIdx)
                * density;

            if (enableEnergy) {
                Scalar specificEnthalpy;
                if (fs.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
                    specificEnthalpy = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
                else 
                    specificEnthalpy = insideVolVars.fluidState().enthalpy(phaseIdx);

                // currently we neglect heat conduction!
                Scalar enthalpyRate = 
                    density
                    * fluxVars.volumeFlux(phaseIdx)
                    * specificEnthalpy;
                EnergyModule::setEnthalpyRate(*this, enthalpyRate);
            }
        }

#ifndef NDEBUG
        for (int i = 0; i < numEq; ++ i) {
            Valgrind::CheckDefined((*this)[i]);
        };
        Valgrind::CheckDefined(*this);
#endif
    }

    /*!
     * \brief Specify an inflow boundary
     */
    template <class Context, class FluidState>
    void setInFlow(const Context &context, 
                   int bfIdx, 
                   int timeIdx,
                   const FluidState &fs)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fs);
        
        // we only allow fluxes in the direction opposite to the outer
        // unit normal
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            Scalar &val = this->operator[](eqIdx);
            val = std::min<Scalar>(0.0, val);
        };
    }

    /*!
     * \brief Specify an outflow boundary
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context &context, 
                    int bfIdx, 
                    int timeIdx,
                    const FluidState &fs)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fs);
        
        // we only allow fluxes in the same direction as the outer
        // unit normal
        for (int eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            Scalar &val = this->operator[](eqIdx);
            val = std::max<Scalar>(0.0, val);
        };
    }
    
    /*!
     * \brief Specify a no-flow boundary.
     */
    void setNoFlow()
    { (*this) = 0.0; }
};

} // end namepace

#endif
