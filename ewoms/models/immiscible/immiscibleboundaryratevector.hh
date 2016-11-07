// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::ImmiscibleBoundaryRateVector
 */
#ifndef EWOMS_IMMISCIBLE_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_IMMISCIBLE_BOUNDARY_RATE_VECTOR_HH

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include "immiscibleintensivequantities.hh"

namespace Ewoms {

/*!
 * \ingroup ImmiscibleModel
 *
 * \brief Implements a boundary vector for the fully implicit
 *        multi-phase model which assumes immiscibility.
 */
template <class TypeTag>
class ImmiscibleBoundaryRateVector : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    ImmiscibleBoundaryRateVector()
        : ParentType()
    {}

    /*!
     * \brief Constructor that assigns all entries to a scalar value.
     *
     * \param value The scalar value to which all components of the
     *              boundary rate vector will be set.
     */
    ImmiscibleBoundaryRateVector(const Evaluation& value)
        : ParentType(value)
    {}

    /*!
     * \brief Copy constructor
     *
     * \param value The boundary rate vector to be duplicated.
     */
    ImmiscibleBoundaryRateVector(const ImmiscibleBoundaryRateVector& value) = default;

    ImmiscibleBoundaryRateVector& operator=(const ImmiscibleBoundaryRateVector& value) = default;

    /*!
     * \brief Specify a free-flow boundary
     *
     * \param context The execution context for which the boundary rate should
     *                be specified.
     * \param bfIdx The local space index of the boundary segment.
     * \param timeIdx The index used by the time discretization.
     * \param fluidState The repesentation of the thermodynamic state
     *                   of the system on the integration point of the
     *                   boundary segment.
     */
    template <class Context, class FluidState>
    void setFreeFlow(const Context& context, unsigned bfIdx, unsigned timeIdx, const FluidState& fluidState)
    {
        typename FluidSystem::template ParameterCache<typename FluidState::Scalar> paramCache;
        paramCache.updateAll(fluidState);

        ExtensiveQuantities extQuants;
        extQuants.updateBoundary(context, bfIdx, timeIdx, fluidState, paramCache);
        const auto& insideIntQuants = context.intensiveQuantities(bfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        (*this) = Evaluation(0.0);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Evaluation density;
            if (fluidState.pressure(phaseIdx) > insideIntQuants.fluidState().pressure(phaseIdx))
                density = FluidSystem::density(fluidState, paramCache, phaseIdx);
            else
                density = insideIntQuants.fluidState().density(phaseIdx);

            Valgrind::CheckDefined(density);
            Valgrind::CheckDefined(extQuants.volumeFlux(phaseIdx));

            // add advective flux of current component in current
            // phase
            (*this)[conti0EqIdx + phaseIdx] += extQuants.volumeFlux(phaseIdx)*density;

            if (enableEnergy) {
                Evaluation specificEnthalpy;
                Scalar pBoundary = fluidState.pressure(phaseIdx);
                const Evaluation& pElement = insideIntQuants.fluidState().pressure(phaseIdx);
                if (pBoundary > pElement)
                    specificEnthalpy = FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
                else
                    specificEnthalpy = insideIntQuants.fluidState().enthalpy(phaseIdx);

                // currently we neglect heat conduction!
                Evaluation enthalpyRate = density*extQuants.volumeFlux(phaseIdx)*specificEnthalpy;
                EnergyModule::addToEnthalpyRate(*this, enthalpyRate);
            }
        }

        EnergyModule::addToEnthalpyRate(*this, EnergyModule::heatConductionRate(extQuants));

#ifndef NDEBUG
        for (unsigned i = 0; i < numEq; ++i)
            Valgrind::CheckDefined((*this)[i]);
        Valgrind::CheckDefined(*this);
#endif
    }

    /*!
     * \brief Specify an inflow boundary
     *
     * An inflow boundary condition is basically a free flow boundary
     * condition that is not prevented from specifying a flow out of
     * the domain.
     *
     * \copydetails setFreeFlow
     */
    template <class Context, class FluidState>
    void setInFlow(const Context& context,
                   unsigned bfIdx,
                   unsigned timeIdx,
                   const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the direction opposite to the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::min(0.0, val);
        }
    }

    /*!
     * \brief Specify an outflow boundary
     *
     * An outflow boundary condition is basically a free flow boundary
     * condition that is not prevented from specifying a flow into
     * the domain.
     *
     * \copydetails setFreeFlow
     */
    template <class Context, class FluidState>
    void setOutFlow(const Context& context,
                    unsigned bfIdx,
                    unsigned timeIdx,
                    const FluidState& fluidState)
    {
        this->setFreeFlow(context, bfIdx, timeIdx, fluidState);

        // we only allow fluxes in the same direction as the outer unit normal
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            Evaluation& val = this->operator[](eqIdx);
            val = Toolbox::max(0.0, val);
        }
    }

    /*!
     * \brief Specify a no-flow boundary for all conserved quantities.
     */
    void setNoFlow()
    { (*this) = Evaluation(0.0); }
};

} // namespace Ewoms

#endif
