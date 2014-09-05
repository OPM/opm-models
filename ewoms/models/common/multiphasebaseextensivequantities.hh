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
 * \copydoc Ewoms::MultiPhaseBaseExtensiveQuantities
 */
#ifndef EWOMS_MULTI_PHASE_BASE_EXTENSIVE_QUANTITIES_HH
#define EWOMS_MULTI_PHASE_BASE_EXTENSIVE_QUANTITIES_HH

#include "multiphasebaseproperties.hh"

#include <ewoms/models/common/quantitycallbacks.hh>
#include <ewoms/disc/common/fvbaseextensivequantities.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/common/fvector.hh>

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief This class calculates the pressure potential gradients and
 *        the filter velocities for multi-phase flow in porous media
 */
template <class TypeTag>
class MultiPhaseBaseExtensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscExtensiveQuantities)
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscExtensiveQuantities) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;
    typedef typename VelocityModule::VelocityExtensiveQuantities VelocityExtensiveQuantities;

public:
    /*!
     * \brief Register all run-time parameters for the extensive quantities.
     */
    static void registerParameters()
    {
        VelocityModule::registerParameters();
    }

    /*!
     * \brief Update the extensive quantities for a given sub-control-volume-face.
     *
     * \param elemCtx Reference to the current element context.
     * \param scvfIdx The local index of the sub-control-volume face for
     *                which the extensive quantities should be calculated.
     * \param timeIdx The index used by the time discretization.
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx);

        const auto &scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);

        // compute the pressure potential gradients
        calculateGradients_(elemCtx, scvfIdx, timeIdx);

        // determine the upstream indices. since this is a semi-smooth non-linear solver,
        // make upstream only look at the evaluation point for the upstream decision
        const auto &evalExtQuants = elemCtx.evalPointExtensiveQuantities(scvfIdx, timeIdx);
        if (&evalExtQuants == this) {
            // we _are_ the evaluation point. Check whether the pressure potential is in
            // the same direction as the face normal or in the opposite one
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                    Valgrind::SetUndefined(upstreamScvIdx_[phaseIdx]);
                    Valgrind::SetUndefined(downstreamScvIdx_[phaseIdx]);
                    continue;
                }

                Scalar tmp = 0;
                for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                    tmp += potentialGrad(phaseIdx)[dimIdx] * scvf.normal()[dimIdx];

                if (tmp > 0) {
                    upstreamScvIdx_[phaseIdx] = this->exteriorIndex();
                    downstreamScvIdx_[phaseIdx] = this->interiorIndex();
                }
                else {
                    upstreamScvIdx_[phaseIdx] = this->interiorIndex();
                    downstreamScvIdx_[phaseIdx] = this->exteriorIndex();
                }
            }
        }
        else {
            // we are *not* the evaluation point. in this case, we just take the
            // up-/downstream indices from the evaluation point.
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                upstreamScvIdx_[phaseIdx] = evalExtQuants.upstreamIndex(phaseIdx);
                downstreamScvIdx_[phaseIdx] = evalExtQuants.downstreamIndex(phaseIdx);
            }
        }

        VelocityExtensiveQuantities::calculateVelocities_(elemCtx, scvfIdx, timeIdx);
    }


    /*!
     * \brief Update the extensive quantities for a given boundary face.
     *
     * \param context Reference to the current execution context.
     * \param bfIdx The local index of the boundary face for which
     *              the extensive quantities should be calculated.
     * \param timeIdx The index used by the time discretization.
     * \param fluidState The FluidState on the domain boundary.
     * \param paramCache The FluidSystem's parameter cache.
     */
    template <class Context, class FluidState>
    void updateBoundary(const Context &context,
                        int bfIdx,
                        int timeIdx,
                        const FluidState &fluidState,
                        typename FluidSystem::ParameterCache &paramCache)
    {
        ParentType::updateBoundary(context, bfIdx, timeIdx, fluidState, paramCache);

        calculateBoundaryGradients_(context, bfIdx, timeIdx, fluidState, paramCache);
        VelocityExtensiveQuantities::calculateBoundaryVelocities_(context,
                                                            bfIdx,
                                                            timeIdx,
                                                            fluidState,
                                                            paramCache);
    }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase
     *        at the face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const DimVector &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    /*!
     * \brief Return the local index of the upstream control volume for a given phase as
     *        a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    short upstreamIndex(int phaseIdx) const
    { return upstreamScvIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    short downstreamIndex(int phaseIdx) const
    { return downstreamScvIdx_[phaseIdx]; }

    /*!
     * \brief Return the weight of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar upstreamWeight(int phaseIdx) const
    { return 1.0; }

    /*!
     * \brief Return the weight of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar downstreamWeight(int phaseIdx) const
    { return 1.0 - upstreamWeight(phaseIdx); }

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int faceIdx,
                             int timeIdx)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        Ewoms::PressureCallback<TypeTag> pressureCallback(elemCtx);

        // calculate the pressure gradient
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                Valgrind::SetUndefined(potentialGrad_[phaseIdx]);
                continue;
            }

            pressureCallback.setPhaseIndex(phaseIdx);
            gradCalc.calculateGradient(potentialGrad_[phaseIdx],
                                       elemCtx,
                                       faceIdx,
                                       pressureCallback);
        }

        // correct the pressure gradients by the gravitational acceleration
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            auto g = elemCtx.problem().gravity(elemCtx, this->interiorIndex(), timeIdx);
            const auto& gEx = elemCtx.problem().gravity(elemCtx, this->exteriorIndex(), timeIdx);
            g += gEx;
            g /= 2;
            Valgrind::CheckDefined(g);

            const auto &intQuantsIn = elemCtx.intensiveQuantities(this->interiorIndex(), timeIdx);
            const auto &intQuantsEx = elemCtx.intensiveQuantities(this->exteriorIndex(), timeIdx);

            const auto &posIn = elemCtx.pos(this->interiorIndex(), timeIdx);
            const auto &posEx = elemCtx.pos(this->exteriorIndex(), timeIdx);

            // the distance between the centers of the control volumes
            DimVector distVec(posEx);
            distVec -= posIn;
            Scalar absDistSquared = distVec.two_norm2();

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate the hydrostatic pressure at the reference height, i.e., at
                // the surface...
                Scalar rhoIn = intQuantsIn.fluidState().density(phaseIdx);
                Scalar rhoEx = intQuantsEx.fluidState().density(phaseIdx);
                Scalar rho = (rhoIn + rhoEx)/2;

                // calculate the difference in height of the exterior SCV compared to
                // interior one.
                Scalar x = distVec * g;
                x /= g*g;

                DimVector deltaHeightVec = g;
                deltaHeightVec *= x;

                // compute the hydrostatic pressure for the exterior control volume
                // assuming constant density...
                Scalar pStatIn = 0;
                Scalar pStatEx = pStatIn - rho*(g*deltaHeightVec);

                // compute the hydrostatic gradient between the two control volumes (this
                // gradient exhibitis the same direction as the vector between the two
                // control volume centers and the length (pStaticExterior -
                // pStaticInterior)/distanceInteriorToExterior
                auto f(distVec);
                f *= (pStatEx - pStatIn)/(absDistSquared);

                // calculate the final potential gradient
                potentialGrad_[phaseIdx] += f;
                if (!std::isfinite(potentialGrad_[phaseIdx].two_norm())) {
                    OPM_THROW(Opm::NumericalProblem,
                               "Non finite potential gradient for phase '"
                               << FluidSystem::phaseName(phaseIdx) << "'");
                }
            }
        }
    }

    template <class Context, class FluidState>
    void calculateBoundaryGradients_(const Context &context,
                                     int bfIdx,
                                     int timeIdx,
                                     const FluidState &fluidState,
                                     const typename FluidSystem::ParameterCache &paramCache)
    {
        const auto& gradCalc = context.gradientCalculator();
        Ewoms::BoundaryPressureCallback<TypeTag, FluidState>
            pressureCallback(context.elementContext(), fluidState);

        // calculate the pressure gradient using two-point gradient
        // appoximation
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!context.model().phaseIsConsidered(phaseIdx)) {
                potentialGrad_[phaseIdx] = 0;
                continue;
            }

            pressureCallback.setPhaseIndex(phaseIdx);
            gradCalc.calculateBoundaryGradient(potentialGrad_[phaseIdx],
                                               context.elementContext(),
                                               bfIdx,
                                               pressureCallback);
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            DimVector g(context.problem().gravity(context.elementContext(),
                                                  this->interiorIndex(),
                                                  timeIdx));

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                if (!context.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate volumetric gravity acceleration force
                DimVector f(g);
                f *= context.intensiveQuantities(bfIdx, timeIdx).fluidState().density(phaseIdx);

                // calculate the final potential gradient
                potentialGrad_[phaseIdx] -= f;
            }
        }
    }

    short upstreamScvIdx_[numPhases];
    short downstreamScvIdx_[numPhases];

    // pressure potential gradients of all phases
    DimVector potentialGrad_[numPhases];
};

} // namespace Ewoms

#endif
