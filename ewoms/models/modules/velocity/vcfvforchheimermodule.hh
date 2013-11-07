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
 * \brief This file contains the necessary classes to calculate the
 *        velocity out of a pressure potential gradient using the
 *        Forchhheimer approach.
 */
#ifndef EWOMS_VCFV_FORCHHEIMER_MODULE_HH
#define EWOMS_VCFV_FORCHHEIMER_MODULE_HH

#include <ewoms/disc/vcfv/vcfvproperties.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(MaterialLaw);
}}

namespace Ewoms {
template <class TypeTag>
class VcfvForchheimerVolumeVariables;

template <class TypeTag>
class VcfvForchheimerFluxVariables;

template <class TypeTag>
class VcfvForchheimerBaseProblem;

/*!
 * \ingroup VcfvForchheimerVelocity
 * \brief Specifies a velocity module which uses the Forchheimer relation.
 */
template <class TypeTag>
struct VcfvForchheimerVelocityModule
{
    typedef VcfvForchheimerVolumeVariables<TypeTag> VelocityVolumeVariables;
    typedef VcfvForchheimerFluxVariables<TypeTag> VelocityFluxVariables;
    typedef VcfvForchheimerBaseProblem<TypeTag> VelocityBaseProblem;

    /*!
     * \brief Register all run-time parameters for the velocity module.
     */
    static void registerParameters()
    { }
};

/*!
 * \ingroup VcfvForchheimerVelocity
 * \brief Provides the defaults for the parameters required by the
 *        Forchheimer velocity approach.
 */
template <class TypeTag>
class VcfvForchheimerBaseProblem
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:

    /*!
     * \brief Returns the Ergun coefficient.
     *
     * The Ergun coefficient is a measure how much the velocity is
     * reduced by turbolence. It is a quantity that does not depend on
     * the fluid phase but only on the porous medium in question. A
     * value of 0 means that the velocity is not influenced by
     * turbolence.
     */
    template <class Context>
    Scalar ergunCoefficient(Context &context, int spaceIdx, int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                  "Not implemented: Problem::ergunCoefficient()");
    }

    /*!
     * \brief Returns the ratio between the phase mobility \f$k_{r,\alpha}\f$ and its passability
     *        \f$\eta_{r,\alpha}\f$ for a given fluid phase \f$\alpha\f$.
     *
     * The passability coefficient specifies the influence of the
     * other fluid phases on the turbolent behaviour of a given fluid
     * phase. By default it is equal to the relative permeability. The
     * mobility to passability ratio is the inverse of phase' the viscosity.
     */
    template <class Context>
    Scalar mobilityPassabilityRatio(Context &context, int spaceIdx, int timeIdx, int phaseIdx) const
    { return 1.0/context.volVars(spaceIdx, timeIdx).fluidState().viscosity(phaseIdx); }
};

/*!
 * \ingroup VcfvForchheimerVelocity
 * \brief Provides the volume variables for the Forchheimer module
 */
template <class TypeTag>
class VcfvForchheimerVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

public:
    /*!
     * \brief Returns the Ergun coefficient.
     *
     * The Ergun coefficient is a measure how much the velocity is
     * reduced by turbolence. A value of 0 means that it is not
     * influenced.
     */
    Scalar ergunCoefficient() const
    { return ergunCoefficient_; };

    /*!
     * \brief Returns the passability of a phase.
     */
    Scalar mobilityPassabilityRatio(int phaseIdx) const
    { return mobilityPassabilityRatio_[phaseIdx]; };

protected:
    void update_(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        const auto &problem = elemCtx.problem();
        ergunCoefficient_ = problem.ergunCoefficient(elemCtx, scvIdx, timeIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mobilityPassabilityRatio_[phaseIdx] = problem.mobilityPassabilityRatio(elemCtx, scvIdx, timeIdx, phaseIdx);
    }

private:
    Scalar ergunCoefficient_;
    Scalar mobilityPassabilityRatio_[numPhases];
};

/*!
 * \ingroup VcfvForchheimerVelocity
 * \brief Provides the Forchheimer velocity module
 *
 * The commonly used Darcy relation looses its validity for Reynolds
 * numbers \f$ Re > 1\f$.  If one encounters flow velocities in porous
 * media above this threshold, the Forchheimer relation can be
 * used. Like the Darcy relation, it relates the gradient in potential
 * to velocity.  However, this relation is not linear (as in the Darcy
 * case) any more.
 *
 * Therefore, the Newton scheme is used to solve the non-linear
 * relation. This velocity is then used like the Darcy velocity
 * e.g. by the local residual.
 *
 * For Reynolds numbers above \f$\approx 500\f$ the standard Forchheimer
 * relation also looses it's validity.
 *
 * The Forchheimer equation is given by the following relation:
 *
 * \f[
  \nabla p_\alpha - \rho_\alpha \vec{g} =
  - \frac{\mu_\alpha}{k_{r,\alpha}} K^{-1}\vec{v}_\alpha
  - \frac{\rho_\alpha C_E}{\eta_{r,\alpha}} \sqrt{K}^{-1}
  \left| \vec{v}_\alpha \right| \vec{v}_\alpha
 \f]
 *
 * Where \f$C_E\f$ is the modified Ergun parameter and
 * \f$\eta_{r,\alpha}\f$ is the passability which is given by a
 * closure relation (usually it is assumed to be identical to the
 * relative permeability). To avoid numerical problems, the relation
 * implemented by this class multiplies both sides with
 * \f$\frac{k_{r_alpha}}{mu} K\f$, so we get
 *
 * \f[
  \frac{k_{r_alpha}}{mu} K \left( \nabla p_\alpha - \rho_\alpha \vec{g}\right) =
  - \vec{v}_\alpha
  - \frac{\rho_\alpha C_E}{\eta_{r,\alpha}}  \frac{k_{r_alpha}}{mu} \sqrt{K}
  \left| \vec{v}_\alpha \right| \vec{v}_\alpha
 \f]

 */
template <class TypeTag>
class VcfvForchheimerFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

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
     * \brief Return the Ergun coefficent at the face's integration
     *        point.
     */
    Scalar ergunCoefficient() const
    { return ergunCoefficient_; }

    /*!
     * \brief Return the ratio of the mobility divided by the
     *        passability at the face's integration point for a given
     *        fluid phase.
     *
     * Usually, that's the inverse of the viscosity.
     */
    Scalar mobilityPassabilityRatio(int phaseIdx) const
    { return mobilityPassabilityRatio_[phaseIdx]; }

    /*!
     * \brief Return the filter velocity of a fluid phase at the
     *        face's integration point [m/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const DimVector &filterVelocity(int phaseIdx) const
    { return filterVelocity_[phaseIdx]; }

    /*!
     * \brief Return the volume flux of a fluid phase at the
     *        face's integration point \f$[m^3/s]\f$
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar volumeFlux(int phaseIdx) const
    { return volumeFlux_[phaseIdx]; }

protected:
    /*!
     * \brief Calculate the filter velocities of all phases
     *
     * The pressure potentials and upwind directions must already be
     * determined before calling this method!
     */
    void calculateVelocities_(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        const auto &problem = elemCtx.problem();

        const auto &volVarsI = elemCtx.volVars(asImp_().insideIndex(), timeIdx);
        const auto &volVarsJ = elemCtx.volVars(asImp_().outsideIndex(), timeIdx);

        // calculate the intrinsic permeability
        const auto &Ki = volVarsI.intrinsicPermeability();
        const auto &Kj = volVarsJ.intrinsicPermeability();
        problem.meanK(K_, Ki, Kj);
        Valgrind::CheckDefined(K_);

        assert(isDiagonal_(K_));
        sqrtK_ = 0.0;
        for (int i = 0; i < dimWorld; ++i)
            sqrtK_[i][i] = std::sqrt(K_[i][i]);

        const DimVector &normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;
        Valgrind::CheckDefined(normal);

        // obtain the Ergun coefficient from the volume
        // variables. Until a better method comes along, we use
        // arithmetic averaging.
        ergunCoefficient_ =
            (volVarsI.ergunCoefficient()
             + volVarsJ.ergunCoefficient())
            / 2;

        ///////////////
        // calculate the weights of the upstream and the downstream
        // control volumes
        ///////////////
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                filterVelocity_[phaseIdx] = 0;
                volumeFlux_[phaseIdx] = 0;
                continue;
            }

            const auto &up = elemCtx.volVars(asImp_().upstreamIndex(phaseIdx), timeIdx);
            density_[phaseIdx] = up.fluidState().density(phaseIdx);
            mobility_[phaseIdx] = up.mobility(phaseIdx);
            mobilityPassabilityRatio_[phaseIdx] = up.mobilityPassabilityRatio(phaseIdx);

            calculateForchheimerVelocity_(phaseIdx);
            volumeFlux_[phaseIdx] = filterVelocity_[phaseIdx] * normal;
        }
    }

    /*!
     * \brief Calculate the filter velocities of all phases at the domain boundary
     */
    template <class Context, class FluidState>
    void calculateBoundaryVelocities_(const Context &context,
                                      int bfIdx,
                                      int timeIdx,
                                      const FluidState &fluidState,
                                      const typename FluidSystem::ParameterCache &paramCache)
    {
        const auto &elemCtx = context.elemContext();
        const auto &problem = elemCtx.problem();
        int insideScvIdx = asImp_().insideIndex();
        const auto &volVarsInside = elemCtx.volVars(insideScvIdx, timeIdx);
        const auto &fsInside = volVarsInside.fluidState();

        // calculate the intrinsic permeability
        K_ = volVarsInside.intrinsicPermeability();

        assert(isDiagonal_(K_));
        sqrtK_ = 0.0;
        for (int i = 0; i < dimWorld; ++i)
            sqrtK_[i][i] = std::sqrt(K_[i][i]);

        DimVector normal = context.fvElemGeom(timeIdx).boundaryFace[bfIdx].normal;

        const auto &matParams = problem.materialLawParams(elemCtx, insideScvIdx, timeIdx);

        Scalar kr[numPhases];
        MaterialLaw::relativePermeabilities(kr, matParams, fluidState);

        // obtain the Ergun coefficient. Because we are on the
        // boundary here, we will take the Ergun coefficient of
        // the interior sub-control volume
        ergunCoefficient_ = volVarsInside.ergunCoefficient();

        ///////////////
        // calculate the weights of the upstream and the downstream
        // control volumes
        ///////////////
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                filterVelocity_[phaseIdx] = 0;
                volumeFlux_[phaseIdx] = 0;
                continue;
            }

            // calculate the actual darcy velocities by multiplying
            // the current "filter velocity" with the upstream mobility
            if (fsInside.pressure(phaseIdx) < fluidState.pressure(phaseIdx)) {
                // the outside of the domain has higher pressure. we
                // need to calculate the mobility
                mobility_[phaseIdx] = kr[phaseIdx]/FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
                density_[phaseIdx] = FluidSystem::density(fluidState, paramCache, phaseIdx);
            }
            else {
                mobility_[phaseIdx] = volVarsInside.mobility(phaseIdx);
                density_[phaseIdx] = fsInside.density(phaseIdx);
            }

            // take the mobility/passability ratio of the phase from the inside SCV
            mobilityPassabilityRatio_[phaseIdx] = volVarsInside.mobilityPassabilityRatio(phaseIdx);

            calculateForchheimerVelocity_(phaseIdx);
            volumeFlux_[phaseIdx] = filterVelocity_[phaseIdx] * normal;
        }
    }

    void calculateForchheimerVelocity_(int phaseIdx)
    {
        // initial guess: filter velocity is zero
        DimVector &velocity = filterVelocity_[phaseIdx];
        velocity = 0;

        DimVector deltaV(1e5); // the change of velocity between two consecutive Newton iterations
        DimVector residual; // the residual (function value that is to be minimized ) of the equation that is to be fulfilled
        DimMatrix gradResid; // slope of equation that is to be solved

        // search by means of the Newton method for a root of Forchheimer equation
        int newtonIter = 0;
        while (deltaV.two_norm() > 1e-11) {
            if (newtonIter >= 30)
                OPM_THROW(Opm::NumericalProblem,
                           "Could not determine Forchheimer velocity within "
                           << newtonIter <<" iterations");
            ++newtonIter;

            // calculate the residual and its Jacobian matrix
            gradForchheimerResid_(residual, gradResid, phaseIdx);

            // newton method
            gradResid.solve(deltaV, residual);
            velocity -= deltaV;
        }

    };

    void forchheimerResid_(DimVector &residual,
                           int phaseIdx) const
    {
        const auto &imp = asImp_();

        const DimVector &velocity = filterVelocity_[phaseIdx];

        // Obtaining the upstreamed quantities
        Scalar mobility = mobility_[phaseIdx];
        Scalar density = density_[phaseIdx];
        Scalar mobilityPassabilityRatio = mobilityPassabilityRatio_[phaseIdx];

        // optain the quantites for the integration point
        const auto &pGrad = imp.potentialGrad(phaseIdx);

        // residual = v_\alpha
        residual = velocity;

        // residual += mobility_\alpha K(\grad p_\alpha - \rho_\alpha g)
        K_.usmv(mobility, pGrad, residual);

        // Forchheimer turbulence correction:
        //
        // residual += \rho_\alpha * mobility_\alpha * C_E / \eta_{r,\alpha} * abs(v_\alpha) * sqrt(K)*v_\alpha
        sqrtK_.usmv(density
                    * mobilityPassabilityRatio
                    * ergunCoefficient_
                    * velocity.two_norm(),
                    velocity,
                    residual);

        Valgrind::CheckDefined(residual);
    }

    void gradForchheimerResid_(DimVector &residual,
                               DimMatrix &gradResid,
                               int phaseIdx)
    {
        DimVector &velocity = filterVelocity_[phaseIdx];
        forchheimerResid_(residual, phaseIdx);

        Scalar eps = 1e-11;
        DimVector tmp;
        for (int i = 0; i < dimWorld; ++i) {
            Scalar coordEps = std::max(eps, velocity[i]*(1 + eps));
            velocity[i] += coordEps;
            forchheimerResid_(tmp, phaseIdx);
            tmp -= residual;
            tmp /= coordEps;
            gradResid[i] = tmp;
            velocity[i] -= coordEps;
        }
    }

    /*!
     * \brief Check whether all off-diagonal entries of a tensor are zero.
     *
     * \param K the tensor that is to be checked.
     * \return True iff all off-diagonals are zero.
     *
     */
    bool isDiagonal_(const DimMatrix &K) const
    {
        for (int i = 0; i < dimWorld; i++) {
            for (int j = 0; j < dimWorld; j++){
                if (i==j)
                    continue;

                if (std::abs(K[i][j]) > 1e-25)
                    return false;
            }
        }
        return true;
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

protected:
    // intrinsic permeability tensor and its square root
    DimMatrix K_;
    DimMatrix sqrtK_;

    // Ergun coefficient of all phases at the integration point
    Scalar ergunCoefficient_;

    // Passability of all phases at the integration point
    Scalar mobilityPassabilityRatio_[numPhases];

    // Density of all phases at the integration point
    Scalar density_[numPhases];

    // Mobility of all phases at the integration point
    Scalar mobility_[numPhases];

    // filter velocities of all phases [m/s]
    DimVector filterVelocity_[numPhases];

    // the volumetric flux of all fluid phases over the control
    // volume's face [m^3/s]
    Scalar volumeFlux_[numPhases];
};

} // namespace Ewoms

#endif
