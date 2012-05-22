// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \brief Implements a boundary vector for the fully implicit single-phase model.
 */
#ifndef DUMUX_BOX_1P2C_BOUNDARY_RATE_VECTOR_HH
#define DUMUX_BOX_1P2C_BOUNDARY_RATE_VECTOR_HH

#include "1p2cproperties.hh"

#include <dumux/common/valgrind.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup 2PModel
 *
 * \brief Implements a boundary vector for the fully implicit single-phase model.
 */
template <class TypeTag>
class OnePTwoCBoundaryRateVector
    : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { contiEqIdx = Indices::contiEqIdx };
    enum { transEqIdx = Indices::transEqIdx };

public:
    /*!
     * \brief Default constructor
     */
    OnePTwoCBoundaryRateVector()
        : ParentType()
    { }

    /*!
     * \brief Constructor with assignment from scalar
     */
    OnePTwoCBoundaryRateVector(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \brief Copy constructor
     */
    OnePTwoCBoundaryRateVector(const OnePTwoCBoundaryRateVector &value)
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
            Scalar molarDensity, x1;
            if (fs.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx)) {
                x1 = fs.moleFraction(phaseIdx, /*compIdx=*/1);
                Scalar meanM =
                    (1 - x1) * FluidSystem::molarMass(0) 
                    + x1 * FluidSystem::molarMass(1);

                molarDensity = FluidSystem::density(fs, paramCache, phaseIdx)/meanM;
            }
            else  {
                molarDensity = insideVolVars.fluidState().molarDensity(phaseIdx);
                x1 = insideVolVars.fluidState().moleFraction(phaseIdx, /*compIdx=*/1);
            }

            // add advective flux of current component in current
            // phase
            (*this)[contiEqIdx] +=
                fluxVars.filterVelocityNormal(phaseIdx)
                * molarDensity;

            (*this)[transEqIdx] +=
                fluxVars.filterVelocityNormal(phaseIdx)
                * molarDensity
                * x1;
        }
        
#ifndef NDEBUG
        for (int i = 0; i < numEq; ++ i) {
            Valgrind::CheckDefined((*this)[i]);
        };
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
