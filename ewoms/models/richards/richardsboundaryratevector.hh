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
 * \copydoc Ewoms::RichardsBoundaryRateVector
 */
#ifndef EWOMS_RICHARDS_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_RICHARDS_BOUNDARY_RATE_VECTOR_HH

#include <opm/common/valgrind.hh>
#include <opm/material/constraintsolvers/ncpflash.hh>

#include "richardsvolumevariables.hh"

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 *
 * \brief Implements a boundary vector for the fully implicit Richards model.
 */
template <class TypeTag>
class RichardsBoundaryRateVector
    : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { contiWEqIdx = Indices::contiWEqIdx };
    enum { wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };

public:
    RichardsBoundaryRateVector()
        : ParentType()
    { }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    RichardsBoundaryRateVector(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector &)
     */
    RichardsBoundaryRateVector(const RichardsBoundaryRateVector &value)
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

        int phaseIdx = wPhaseIdx;
        Scalar density;
        if (fluidState.pressure(phaseIdx) > insideVolVars.fluidState().pressure(phaseIdx))
            density = FluidSystem::density(fluidState, paramCache, phaseIdx);
        else
            density = insideVolVars.fluidState().density(phaseIdx);

        // add advective flux of current component in current
        // phase
        (*this)[contiWEqIdx] +=
            fluxVars.volumeFlux(phaseIdx)
            * density;

#ifndef NDEBUG
        for (int i = 0; i < numEq; ++ i) {
            Valgrind::CheckDefined((*this)[i]);
        };
        Valgrind::CheckDefined(*this);
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

} // end namepace

#endif
