// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \copydoc Ewoms::RichardsBoundaryRateVector
 */
#ifndef EWOMS_RICHARDS_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_RICHARDS_BOUNDARY_RATE_VECTOR_HH

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include "richardsintensivequantities.hh"

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 *
 * \brief Implements a boundary vector for the fully implicit Richards model.
 */
template <class TypeTag>
class RichardsBoundaryRateVector : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { contiEqIdx = Indices::contiEqIdx };
    enum { liquidPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };

    typedef Opm::MathToolbox<Evaluation> Toolbox;

public:
    RichardsBoundaryRateVector() : ParentType()
    {}

    /*!
     * \copydoc
     * ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    RichardsBoundaryRateVector(const Evaluation& value)
        : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const
     * ImmiscibleBoundaryRateVector &)
     */
    RichardsBoundaryRateVector(const RichardsBoundaryRateVector &value)
        : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setFreeFlow
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context &context, int bfIdx, int timeIdx, const FluidState &fluidState)
    {
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);

        ExtensiveQuantities extQuants;
        extQuants.updateBoundary(context, bfIdx, timeIdx, fluidState, paramCache);
        const auto &insideIntQuants = context.intensiveQuantities(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = Toolbox::createConstant(0.0);

        int phaseIdx = liquidPhaseIdx;
        Evaluation density;
        if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx))
            density = FluidSystem::density(fluidState, paramCache, phaseIdx);
        else
            density = insideIntQuants.fluidState().density(phaseIdx);

        // add advective flux of current component in current
        // phase
        (*this)[contiEqIdx] += extQuants.volumeFlux(phaseIdx) * density;

#ifndef NDEBUG
        for (int i = 0; i < numEq; ++i) {
            Valgrind::CheckDefined((*this)[i]);
        }
        Valgrind::CheckDefined(*this);
#endif
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setInFlow
     */
    template <class Context, class FluidState>
    void setInFlow(const Context &context, int bfIdx, int timeIdx,
                   const FluidState &fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the direction opposite to the outer unit normal
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::min(0.0, val);
        }
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setOutFlow
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context &context, int bfIdx, int timeIdx,
                    const FluidState &fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the same direction as the outer unit normal
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::max(0.0, val);
        }
    }

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::setNoFlow
     */
    void setNoFlow()
    { (*this) = Toolbox::createConstant(0.0); }
};

} // namespace Ewoms

#endif
