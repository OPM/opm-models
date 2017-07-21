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
 * \copydoc Ewoms::StokesBoundaryRateVector
 */
#ifndef EWOMS_STOKES_BOUNDARY_RATE_VECTOR_HH
#define EWOMS_STOKES_BOUNDARY_RATE_VECTOR_HH

#include <opm/material/densead/Math.hpp>
#include <opm/common/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include <dune/common/fvector.hh>

#include "stokesintensivequantities.hh"

namespace Ewoms {

/*!
 * \ingroup StokesModel
 *
 * \brief Implements a boundary vector for the fully implicit (Navier-)Stokes
 *        model.
 */
template <class TypeTag>
class StokesBoundaryRateVector : public GET_PROP_TYPE(TypeTag, RateVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numComponents = FluidSystem::numComponents };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { dimWorld = GridView::dimensionworld };

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { momentum0EqIdx = Indices::momentum0EqIdx };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;

public:
    StokesBoundaryRateVector() : ParentType()
    {}

    /*!
     * \copydoc
     * ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(Scalar)
     */
    StokesBoundaryRateVector(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleBoundaryRateVector::ImmiscibleBoundaryRateVector(const
     * ImmiscibleBoundaryRateVector& )
     */
    StokesBoundaryRateVector(const StokesBoundaryRateVector& value) = default;
    StokesBoundaryRateVector& operator=(const StokesBoundaryRateVector& value) = default;

    /*!
     * \param context The execution context for which the boundary rate should be specified.
     * \param bfIdx The local index of the boundary segment (-> local space index).
     * \param timeIdx The index used by the time discretization.
     * \param velocity The velocity vector [m/s] at the boundary.
     * \param fluidState The repesentation of the thermodynamic state
     *                   of the system on the integration point of the
     *                   boundary segment.
     */
    template <class VelocityEval, class Context, class FluidState>
    void setFreeFlow(const Context& context,
                     unsigned bfIdx,
                     unsigned timeIdx,
                     const Dune::FieldVector<VelocityEval, dimWorld>& velocity,
                     const FluidState& fluidState)
    {
        const auto& stencil = context.stencil(timeIdx);
        const auto& scvf = stencil.boundaryFace(bfIdx);

        unsigned insideScvIdx = context.interiorScvIndex(bfIdx, timeIdx);
        unsigned focusDofIdx = context.focusDofIndex();
        //const auto& insideScv = stencil.subControlVolume(insideScvIdx);
        const auto& insideIntQuants = context.intensiveQuantities(bfIdx, timeIdx);

        // the outer unit normal
        const auto& normal = scvf.normal();

        // distance between the center of the SCV and center of the boundary face
        DimVector distVec = context.element().geometry().corner(static_cast<int>(insideScvIdx));
        distVec -= stencil.subControlVolume(insideScvIdx).geometry().center();

        Scalar dist = 0.0;
        Evaluation volumeFlux = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            dist += distVec[dimIdx]*normal[dimIdx];
            volumeFlux -= velocity[dimIdx]*normal[dimIdx];
        }

        // mass fluxes over the boundary
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.updatePhase(fluidState, phaseIdx);

        const auto& density = FluidSystem::density(fluidState, paramCache, phaseIdx);
        const auto& molarDensity = density / Opm::scalarValue(fluidState.averageMolarMass(phaseIdx));
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            if (volumeFlux > 0.0) {
                // outflow
                if (insideScvIdx == focusDofIdx)
                    (*this)[conti0EqIdx + compIdx] =
                        volumeFlux
                        * insideIntQuants.fluidState().molarity(phaseIdx, compIdx);
                else
                    (*this)[conti0EqIdx + compIdx] =
                        Opm::scalarValue(volumeFlux)
                        * Opm::scalarValue(insideIntQuants.fluidState().molarity(phaseIdx, compIdx));
            }
            else {
                // inflow
                (*this)[conti0EqIdx + compIdx] =
                    volumeFlux
                    * molarDensity
                    * fluidState.moleFraction(phaseIdx, compIdx);
            }
        }

        // momentum flux over the boundary
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            if (volumeFlux > 0.0) {
                // outflow
                if (insideScvIdx == focusDofIdx)
                    (*this)[momentum0EqIdx + dimIdx] =
                        - insideIntQuants.fluidState().density(phaseIdx)
                        * normal[dimIdx]
                        * velocity[dimIdx];
                else
                    (*this)[momentum0EqIdx + dimIdx] =
                        - Opm::scalarValue(insideIntQuants.fluidState().density(phaseIdx))
                        * normal[dimIdx]
                        * Opm::scalarValue(velocity[dimIdx]);
            }
            else {
                // inflow
                (*this)[momentum0EqIdx + dimIdx] =
                    - density
                    * normal[dimIdx]
                    * velocity[dimIdx];
            }
        }

        EnergyModule::setEnthalpyRate(*this, fluidState, phaseIdx, volumeFlux);
    }

    /*!
     * \brief Set a in-flow boundary in the (Navier-)Stoke model
     *
     * \param context The execution context for which the boundary rate should
     *                be specified.
     * \param bfIdx The local space index of the boundary segment.
     * \param timeIdx The index used by the time discretization.
     * \param velocity The velocity vector [m/s] at the boundary.
     * \param fluidState The repesentation of the thermodynamic state
     *                   of the system on the integration point of the
     *                   boundary segment.
     */
    template <class Context, class FluidState>
    void setInFlow(const Context& context,
                   unsigned bfIdx,
                   unsigned timeIdx,
                   const DimVector& velocity,
                   const FluidState& fluidState)
    {
        const auto& stencil = context.stencil(timeIdx);
        const auto& scvf = stencil.boundaryFace(bfIdx);

        // the outer unit normal
        const auto& normal = scvf.normal();

        const auto& intQuants = context.intensiveQuantities(bfIdx, timeIdx);

        Scalar velDir = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            velDir += Toolbox::value(velocity[dimIdx])*normal[dimIdx];

        // don't allow velocities to exhibit the same direction as the face normal
        if (velDir > 0.0)
            (*this) = 0.0;
        else {
            const auto& fluidState = intQuants.fluidState();
            setFreeFlow(context, bfIdx, timeIdx, velocity, fluidState);
        }
    }

    /*!
     * \brief Set a out-flow boundary in the (Navier-)Stoke model
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    void setOutFlow(const Context& context, unsigned bfIdx, unsigned timeIdx)
    {
        const auto& stencil = context.stencil(timeIdx);
        const auto& scvf = stencil.boundaryFace(bfIdx);

        // the outer unit normal
        const auto& normal = scvf.normal();

        const auto& intQuants = context.intensiveQuantities(bfIdx, timeIdx);
        const EvalDimVector& velocity = intQuants.velocity();

        Scalar velDir = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
            velDir += Toolbox::value(velocity[dimIdx])*normal[dimIdx];
        }

        // don't allow velocities to exhibit the opposite direction as the face normal
        if (velDir < 0.0)
            (*this) = 0.0;
        else {
            const auto& fluidState = intQuants.fluidState();
            setFreeFlow(context, bfIdx, timeIdx, velocity, fluidState);
        }

        const auto& fluidState = intQuants.fluidState();
        setFreeFlow(context, bfIdx, timeIdx, velocity, fluidState);
    }

    /*!
     * \brief Set a no-flow boundary in the (Navier-)Stoke model
     *
     * \copydoc Doxygen::contextParams
     */
    template <class Context>
    void setNoFlow(const Context& context,
                   unsigned spaceIdx,
                   unsigned timeIdx)
    {
        // due to shear stresses, there usually is a flux of momentum over the boundary
        // even if there is no mass flux. we can thus treat the boundary segment as if it
        // was free-flow but only use the part of the velocity that is perpendicular to
        // the face normal

        const auto& intQuants = context.intensiveQuantities(spaceIdx, timeIdx);
        const auto& stencil = context.stencil(timeIdx);
        const auto& scvf = stencil.boundaryFace(spaceIdx);
        const auto& normal = scvf.normal();
        const auto& velocity = intQuants.velocity();

        auto beta = velocity[0];
        beta = 0.0;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            beta += velocity[dimIdx]*normal[dimIdx];

        auto projectedVelocity = velocity;
        for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
            projectedVelocity[dimIdx] = velocity[dimIdx] - beta*normal[dimIdx];

        const auto& fluidState = intQuants.fluidState();
        setFreeFlow(context, spaceIdx, timeIdx, projectedVelocity, fluidState);
    }
};

} // namespace Ewoms

#endif
