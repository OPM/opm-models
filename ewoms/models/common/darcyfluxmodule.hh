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
 * \brief This file contains the necessary classes to calculate the
 *        volumetric fluxes out of a pressure potential gradient using the
 *        Darcy relation.
 */
#ifndef EWOMS_DARCY_FLUX_MODULE_HH
#define EWOMS_DARCY_FLUX_MODULE_HH

#include "multiphasebaseproperties.hh"
#include <ewoms/models/common/quantitycallbacks.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(MaterialLaw);
}

template <class TypeTag>
class DarcyIntensiveQuantities;

template <class TypeTag>
class DarcyExtensiveQuantities;

template <class TypeTag>
class DarcyBaseProblem;

/*!
 * \ingroup FluxModules
 * \brief Specifies a flux module which uses the Darcy relation.
 */
template <class TypeTag>
struct DarcyFluxModule
{
    typedef DarcyIntensiveQuantities<TypeTag> FluxIntensiveQuantities;
    typedef DarcyExtensiveQuantities<TypeTag> FluxExtensiveQuantities;
    typedef DarcyBaseProblem<TypeTag> FluxBaseProblem;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    { }
};

/*!
 * \ingroup FluxModules
 * \brief Provides the defaults for the parameters required by the
 *        Darcy velocity approach.
 */
template <class TypeTag>
class DarcyBaseProblem
{ };

/*!
 * \ingroup FluxModules
 * \brief Provides the intensive quantities for the Darcy flux module
 */
template <class TypeTag>
class DarcyIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
protected:
    void update_(const ElementContext &elemCtx, int dofIdx, int timeIdx)
    { }
};

/*!
 * \ingroup FluxModules
 * \brief Provides the Darcy flux module
 *
 * The commonly used Darcy relation looses its validity for Reynolds numbers \f$ Re <
 * 1\f$.  If one encounters flow velocities in porous media above this threshold, the
 * Forchheimer approach can be used.
 *
 * The Darcy equation is given by the following relation:
 *
 * \f[
  \vec{v}_\alpha =
  \left( \nabla p_\alpha - \rho_\alpha \vec{g}\right)
  \frac{\mu_\alpha}{k_{r,\alpha} K}
 \f]
 */
template <class TypeTag>
class DarcyExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename Opm::MathToolbox<Evaluation> Toolbox;
    typedef typename FluidSystem::template ParameterCache<Evaluation> ParameterCache;
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \brief Returns the intrinsic permeability tensor for a given
     *        sub-control volume face.
     */
    const DimMatrix &intrinsicPermability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase
     *        at the face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    /*!
     * \brief Return the filter velocity of a fluid phase at the
     *        face's integration point [m/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector &filterVelocity(int phaseIdx) const
    { return filterVelocity_[phaseIdx]; }

    /*!
     * \brief Return the volume flux of a fluid phase at the face's integration point
     *        \f$[m^3/s / m^2]\f$
     *
     * This is the fluid volume of a phase per second and per square meter of face
     * area.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Evaluation& volumeFlux(int phaseIdx) const
    { return volumeFlux_[phaseIdx]; }

protected:
    short upstreamIndex_(int phaseIdx) const
    { return upstreamDofIdx_[phaseIdx]; }

    short downstreamIndex_(int phaseIdx) const
    { return downstreamDofIdx_[phaseIdx]; }

    /*!
     * \brief Calculate the gradients which are required to determine the volumetric fluxes
     *
     * The the upwind directions is also determined by method.
     */
    void calculateGradients_(const ElementContext &elemCtx,
                             int faceIdx,
                             int timeIdx)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        Ewoms::PressureCallback<TypeTag> pressureCallback(elemCtx);

        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(faceIdx);
        const auto &faceNormal = scvf.normal();

        interiorDofIdx_ = scvf.interiorIndex();
        exteriorDofIdx_ = scvf.exteriorIndex();

        // calculate the "raw" pressure gradient
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
            Valgrind::CheckDefined(potentialGrad_[phaseIdx]);
        }

        // correct the pressure gradients by the gravitational acceleration
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            const auto& gIn = elemCtx.problem().gravity(elemCtx, interiorDofIdx_, timeIdx);
            const auto& gEx = elemCtx.problem().gravity(elemCtx, exteriorDofIdx_, timeIdx);

            const auto &intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);
            const auto &intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx_, timeIdx);

            const auto &posIn = elemCtx.pos(interiorDofIdx_, timeIdx);
            const auto &posEx = elemCtx.pos(exteriorDofIdx_, timeIdx);
            const auto &posFace = scvf.integrationPos();

            // the distance between the centers of the control volumes
            DimVector distVecIn(posIn);
            DimVector distVecEx(posEx);
            DimVector distVecTotal(posEx);

            distVecIn -= posFace;
            distVecEx -= posFace;
            distVecTotal -= posIn;
            Scalar absDistTotalSquared = distVecTotal.two_norm2();
            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate the hydrostatic pressure at the integration point of the face
                auto rhoIn = intQuantsIn.fluidState().density(phaseIdx);
                auto pStatIn = - rhoIn*(gIn*distVecIn);

                // the quantities on the exterior side of the face do not influence the
                // result for the TPFA scheme, so they can be treated as scalar values.
                Scalar rhoEx = Toolbox::value(intQuantsEx.fluidState().density(phaseIdx));
                Scalar pStatEx = - rhoEx*(gEx*distVecEx);

                // compute the hydrostatic gradient between the two control volumes (this
                // gradient exhibitis the same direction as the vector between the two
                // control volume centers and the length (pStaticExterior -
                // pStaticInterior)/distanceInteriorToExterior
                EvalDimVector f(distVecTotal);
                f *= (pStatEx - pStatIn)/absDistTotalSquared;

                // calculate the final potential gradient
                for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                    potentialGrad_[phaseIdx][dimIdx] += f[dimIdx];

                for (unsigned i = 0; i < potentialGrad_[phaseIdx].size(); ++i) {
                    if (!std::isfinite(Toolbox::value(potentialGrad_[phaseIdx][i]))) {
                        OPM_THROW(Opm::NumericalProblem,
                                  "Non-finite potential gradient for phase '"
                                  << FluidSystem::phaseName(phaseIdx) << "'");
                    }
                }
            }
        }

        Valgrind::SetUndefined(K_);
        elemCtx.problem().intersectionIntrinsicPermeability(K_, elemCtx, faceIdx, timeIdx);
        Valgrind::CheckDefined(K_);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                Valgrind::SetUndefined(potentialGrad_[phaseIdx]);
                continue;
            }

            // determine the upstream and downstream DOFs
            Evaluation tmp = 0.0;
            for (unsigned i = 0; i < faceNormal.size(); ++i)
                tmp += potentialGrad_[phaseIdx][i]*faceNormal[i];

            if (tmp > 0) {
                upstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
                downstreamDofIdx_[phaseIdx] = interiorDofIdx_;
            }
            else {
                upstreamDofIdx_[phaseIdx] = interiorDofIdx_;
                downstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
            }

            const auto &up = elemCtx.intensiveQuantities(upstreamDofIdx_[phaseIdx], timeIdx);
            // this is also slightly hacky because it assumes that the derivative of the
            // flux between two DOFs only depends on the primary variables in the
            // upstream direction. For non-TPFA flux approximation schemes, this is not
            // true...
            if (upstreamDofIdx_[phaseIdx] == interiorDofIdx_)
                mobility_[phaseIdx] = up.mobility(phaseIdx);
            else
                mobility_[phaseIdx] = Toolbox::value(up.mobility(phaseIdx));
        }
    }

    /*!
     * \brief Calculate the gradients at the grid boundary which are required to
     *        determine the volumetric fluxes
     *
     * The the upwind directions is also determined by method.
     */
    template <class FluidState>
    void calculateBoundaryGradients_(const ElementContext &elemCtx,
                                     int boundaryFaceIdx,
                                     int timeIdx,
                                     const FluidState& fluidState,
                                     const typename FluidSystem::template ParameterCache<typename FluidState::Scalar>& paramCache)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        Ewoms::BoundaryPressureCallback<TypeTag, FluidState> pressureCallback(elemCtx, fluidState);

        // calculate the pressure gradient
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                Valgrind::SetUndefined(potentialGrad_[phaseIdx]);
                continue;
            }

            pressureCallback.setPhaseIndex(phaseIdx);
            gradCalc.calculateBoundaryGradient(potentialGrad_[phaseIdx],
                                               elemCtx,
                                               boundaryFaceIdx,
                                               pressureCallback);
            Valgrind::CheckDefined(potentialGrad_[phaseIdx]);
        }

        const auto& scvf = elemCtx.stencil(timeIdx).boundaryFace(boundaryFaceIdx);
        interiorDofIdx_ = scvf.interiorIndex();
        exteriorDofIdx_ = -1;

        // calculate the intrinsic permeability
        const auto &intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);
        K_ = intQuantsIn.intrinsicPermeability();

        // correct the pressure gradients by the gravitational acceleration
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            const auto& gIn = elemCtx.problem().gravity(elemCtx, interiorDofIdx_, timeIdx);
            const auto& posIn = elemCtx.pos(interiorDofIdx_, timeIdx);
            const auto& posFace = scvf.integrationPos();

            // the distance between the face center and the center of the control volume
            DimVector distVecIn(posIn);
            distVecIn -= posFace;
            Scalar absDist = distVecIn.two_norm();
            Scalar gTimesDist = gIn*distVecIn;

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate the hydrostatic pressure at the integration point of the face
                Evaluation rhoIn = intQuantsIn.fluidState().density(phaseIdx);
                Evaluation pStatIn = - gTimesDist*rhoIn;

                Valgrind::CheckDefined(pStatIn);

                // compute the hydrostatic gradient between the two control volumes (this
                // gradient exhibitis the same direction as the vector between the two
                // control volume centers and the length (pStaticExterior -
                // pStaticInterior)/distanceInteriorToExterior. Note that for the
                // boundary, 'pStaticExterior' is zero as the boundary pressure is
                // defined on boundary face's integration point...
                EvalDimVector f(distVecIn);
                f *= - pStatIn/absDist;

                // calculate the final potential gradient
                for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                    potentialGrad_[phaseIdx][dimIdx] += f[dimIdx];

                Valgrind::CheckDefined(potentialGrad_[phaseIdx]);
                for (unsigned i = 0; i < potentialGrad_[phaseIdx].size(); ++i) {
                    if (!std::isfinite(Toolbox::value(potentialGrad_[phaseIdx][i]))) {
                        OPM_THROW(Opm::NumericalProblem,
                                  "Non finite potential gradient for phase '"
                                  << FluidSystem::phaseName(phaseIdx) << "'");
                    }
                }
            }
        }

        // determine the upstream and downstream DOFs
        const auto &faceNormal = scvf.normal();

        const auto &matParams = elemCtx.problem().materialLawParams(elemCtx, interiorDofIdx_, timeIdx);

        Scalar kr[numPhases];
        MaterialLaw::relativePermeabilities(kr, matParams, fluidState);
        Valgrind::CheckDefined(kr);

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            Evaluation tmp = 0.0;
            for (unsigned i = 0; i < faceNormal.size(); ++i)
                tmp += potentialGrad_[phaseIdx][i]*faceNormal[i];

            if (tmp > 0) {
                upstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
                downstreamDofIdx_[phaseIdx] = interiorDofIdx_;
            }
            else {
                upstreamDofIdx_[phaseIdx] = interiorDofIdx_;
                downstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
            }

            // take the phase mobility from the DOF in upstream direction
            if (upstreamDofIdx_[phaseIdx] < 0)
                mobility_[phaseIdx] =
                    kr[phaseIdx] / FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            else
                mobility_[phaseIdx] = intQuantsIn.mobility(phaseIdx);
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }
    }

    /*!
     * \brief Calculate the volumetric fluxes of all phases
     *
     * The pressure potentials and upwind directions must already be
     * determined before calling this method!
     */
    void calculateFluxes_(const ElementContext& elemCtx, int scvfIdx, int timeIdx)
    {
        const auto &scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const DimVector &normal = scvf.normal();
        Valgrind::CheckDefined(normal);

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            filterVelocity_[phaseIdx] = 0.0;
            volumeFlux_[phaseIdx] = 0.0;
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                continue;
            }

            asImp_().calculateFilterVelocity_(phaseIdx);
            Valgrind::CheckDefined(filterVelocity_[phaseIdx]);
            volumeFlux_[phaseIdx] = 0.0;
            for (unsigned i = 0; i < normal.size(); ++i)
                volumeFlux_[phaseIdx] += (filterVelocity_[phaseIdx][i] * normal[i]);
        }
    }

    /*!
     * \brief Calculate the volumetric fluxes at a boundary face of all fluid phases
     *
     * The pressure potentials and upwind directions must already be determined before
     * calling this method!
     */
    void calculateBoundaryFluxes_(const ElementContext& elemCtx,
                                      int boundaryFaceIdx,
                                      int timeIdx)
    {
        const auto &scvf = elemCtx.stencil(timeIdx).boundaryFace(boundaryFaceIdx);
        const DimVector &normal = scvf.normal();
        Valgrind::CheckDefined(normal);

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                filterVelocity_[phaseIdx] = 0.0;
                volumeFlux_[phaseIdx] = 0.0;
                continue;
            }

            asImp_().calculateFilterVelocity_(phaseIdx);
            Valgrind::CheckDefined(filterVelocity_[phaseIdx]);
            volumeFlux_[phaseIdx] = 0.0;
            for (unsigned i = 0; i < normal.size(); ++i)
                volumeFlux_[phaseIdx] += (filterVelocity_[phaseIdx][i] * normal[i]);
        }
    }

    void calculateFilterVelocity_(int phaseIdx)
    {
#ifndef NDEBUG
        assert(std::isfinite(Toolbox::value(mobility_[phaseIdx])));
        for (unsigned i = 0; i < K_.M(); ++ i)
            for (unsigned j = 0; j < K_.N(); ++ j)
                assert(std::isfinite(K_[i][j]));
#endif

        K_.mv(potentialGrad_[phaseIdx], filterVelocity_[phaseIdx]);
        filterVelocity_[phaseIdx] *= - mobility_[phaseIdx];

#ifndef NDEBUG
        for (unsigned i = 0; i < filterVelocity_[phaseIdx].size(); ++ i)
            assert(std::isfinite(Toolbox::value(filterVelocity_[phaseIdx][i])));
#endif
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

protected:
    // intrinsic permeability tensor and its square root
    DimMatrix K_;

    // interior, exterior, upstream and downstream DOFs
    short interiorDofIdx_;
    short exteriorDofIdx_;
    short upstreamDofIdx_[numPhases];
    short downstreamDofIdx_[numPhases];

    // mobilities of all fluid phases [1 / (Pa s)]
    Evaluation mobility_[numPhases];

    // filter velocities of all phases [m/s]
    EvalDimVector filterVelocity_[numPhases];

    // the volumetric flux of all fluid phases over the control
    // volume's face [m^3/s / m^2]
    Evaluation volumeFlux_[numPhases];

    // pressure potential gradients of all phases [Pa / m]
    EvalDimVector potentialGrad_[numPhases];
};

} // namespace Ewoms

#endif
