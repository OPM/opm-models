/*
  Copyright (C) 2008-2013 by Andreas Lauser

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
 * \brief This file contains the necessary classes to calculate the
 *        velocity out of a pressure potential gradient using the
 *        Forchhheimer approach.
 */
#ifndef EWOMS_FORCHHEIMER_MODULE_HH
#define EWOMS_FORCHHEIMER_MODULE_HH

#include <ewoms/disc/common/fvbaseproperties.hh>
#include "darcyvelocitymodule.hh"


#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(MaterialLaw);
}
}

namespace Ewoms {
template <class TypeTag>
class ForchheimerIntensiveQuantities;

template <class TypeTag>
class ForchheimerExtensiveQuantities;

template <class TypeTag>
class ForchheimerBaseProblem;

/*!
 * \ingroup ForchheimerVelocity
 * \brief Specifies a velocity module which uses the Forchheimer relation.
 */
template <class TypeTag>
struct ForchheimerVelocityModule
{
    typedef ForchheimerIntensiveQuantities<TypeTag> VelocityIntensiveQuantities;
    typedef ForchheimerExtensiveQuantities<TypeTag> VelocityExtensiveQuantities;
    typedef ForchheimerBaseProblem<TypeTag> VelocityBaseProblem;

    /*!
     * \brief Register all run-time parameters for the velocity module.
     */
    static void registerParameters()
    {}
};

/*!
 * \ingroup ForchheimerVelocity
 * \brief Provides the defaults for the parameters required by the
 *        Forchheimer velocity approach.
 */
template <class TypeTag>
class ForchheimerBaseProblem
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
     * \brief Returns the ratio between the phase mobility
     *        \f$k_{r,\alpha}\f$ and its passability
     *        \f$\eta_{r,\alpha}\f$ for a given fluid phase
     *        \f$\alpha\f$.
     *
     * The passability coefficient specifies the influence of the
     * other fluid phases on the turbolent behaviour of a given fluid
     * phase. By default it is equal to the relative permeability. The
     * mobility to passability ratio is the inverse of phase' the viscosity.
     */
    template <class Context>
    Scalar mobilityPassabilityRatio(Context &context, int spaceIdx, int timeIdx,
                                    int phaseIdx) const
    {
        return 1.0 / context.intensiveQuantities(spaceIdx, timeIdx).fluidState().viscosity(
                         phaseIdx);
    }
};

/*!
 * \ingroup ForchheimerVelocity
 * \brief Provides the intensive quantities for the Forchheimer module
 */
template <class TypeTag>
class ForchheimerIntensiveQuantities
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
    {
        return ergunCoefficient_;
    }

    /*!
     * \brief Returns the passability of a phase.
     */
    Scalar mobilityPassabilityRatio(int phaseIdx) const
    {
        return mobilityPassabilityRatio_[phaseIdx];
    }

protected:
    void update_(const ElementContext &elemCtx, int dofIdx, int timeIdx)
    {
        const auto &problem = elemCtx.problem();
        ergunCoefficient_ = problem.ergunCoefficient(elemCtx, dofIdx, timeIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            mobilityPassabilityRatio_[phaseIdx] =
                problem.mobilityPassabilityRatio(elemCtx,
                                                 dofIdx,
                                                 timeIdx,
                                                 phaseIdx);
    }

private:
    Scalar ergunCoefficient_;
    Scalar mobilityPassabilityRatio_[numPhases];
};

/*!
 * \ingroup ForchheimerVelocity
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
class ForchheimerExtensiveQuantities
    : public DarcyExtensiveQuantities<TypeTag>
{
    typedef DarcyExtensiveQuantities<TypeTag> DarcyExtQuants;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;


    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \brief Return the Ergun coefficent at the face's integration point.
     */
    Scalar ergunCoefficient() const
    { return ergunCoefficient_; }

    /*!
     * \brief Return the ratio of the mobility divided by the passability at the face's
     *        integration point for a given fluid phase.
     *
     * Usually, that's the inverse of the viscosity.
     */
    Scalar mobilityPassabilityRatio(int phaseIdx) const
    { return mobilityPassabilityRatio_[phaseIdx]; }

protected:
    void calculateGradients_(const ElementContext &elemCtx,
                             int faceIdx,
                             int timeIdx)
    {
        DarcyExtQuants::calculateGradients_(elemCtx, faceIdx, timeIdx);

        const auto &intQuantsIn = elemCtx.intensiveQuantities(this->interiorDofIdx_, timeIdx);
        const auto &intQuantsEx = elemCtx.intensiveQuantities(this->exteriorDofIdx_, timeIdx);

        // calculate the square root of the intrinsic permeability
        assert(isDiagonal_(this->K_));
        sqrtK_ = 0.0;
        for (int i = 0; i < dimWorld; ++i)
            sqrtK_[i][i] = std::sqrt(this->K_[i][i]);

        // obtain the Ergun coefficient. Lacking better ideas, we use its the arithmetic mean.
        ergunCoefficient_ = (intQuantsIn.ergunCoefficient() + intQuantsEx.ergunCoefficient())/2;

        // obtain the mobility to passability ratio for each phase.
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            const auto &up = elemCtx.intensiveQuantities(this->upstreamIndex_(phaseIdx), timeIdx);

            density_[phaseIdx] = up.fluidState().density(phaseIdx);
            mobilityPassabilityRatio_[phaseIdx] = up.mobilityPassabilityRatio(phaseIdx);
        }
    }

    template <class FluidState>
    void calculateBoundaryGradients_(const ElementContext &elemCtx,
                                     int boundaryFaceIdx,
                                     int timeIdx,
                                     const FluidState& fluidState,
                                     const typename FluidSystem::ParameterCache& paramCache)
    {
        DarcyExtQuants::calculateBoundaryGradients_(elemCtx,
                                                    boundaryFaceIdx,
                                                    timeIdx,
                                                    fluidState,
                                                    paramCache);

        const auto &intQuantsIn = elemCtx.intensiveQuantities(this->interiorDofIdx_, timeIdx);

        // obtain the Ergun coefficient. Because we are on the boundary here, we will
        // take the Ergun coefficient of the interior
        ergunCoefficient_ = intQuantsIn.ergunCoefficient();

        // calculate the square root of the intrinsic permeability
        assert(isDiagonal_(this->K_));
        sqrtK_ = 0.0;
        for (int i = 0; i < dimWorld; ++i)
            sqrtK_[i][i] = std::sqrt(this->K_[i][i]);

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            density_[phaseIdx] = intQuantsIn.fluidState().density(phaseIdx);
            mobilityPassabilityRatio_[phaseIdx] = intQuantsIn.mobilityPassabilityRatio(phaseIdx);
        }
    }

    /*!
     * \brief Calculate the filter velocities of all phases
     *
     * The pressure potentials and upwind directions must already be
     * determined before calling this method!
     */
    void calculateVelocities_(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        const auto &intQuantsI = elemCtx.intensiveQuantities(asImp_().interiorIndex(), timeIdx);
        const auto &intQuantsJ = elemCtx.intensiveQuantities(asImp_().exteriorIndex(), timeIdx);

        const auto &scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const auto &normal = scvf.normal();
        Valgrind::CheckDefined(normal);

        // obtain the Ergun coefficient from the intensive quantity object. Until a
        // better method comes along, we use arithmetic averaging.
        ergunCoefficient_ = (intQuantsI.ergunCoefficient() + intQuantsJ.ergunCoefficient()) / 2;

        ///////////////
        // calculate the weights of the upstream and the downstream control volumes
        ///////////////
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                this->filterVelocity_[phaseIdx] = 0;
                this->volumeFlux_[phaseIdx] = 0;
                continue;
            }

            calculateForchheimerVelocity_(phaseIdx);
            this->volumeFlux_[phaseIdx] = (this->filterVelocity_[phaseIdx] * normal);
        }
    }

    /*!
     * \brief Calculate the filter velocities of all phases at the domain
     * boundary
     */
    void calculateBoundaryVelocities_(const ElementContext& elemCtx,
                                      int bfIdx,
                                      int timeIdx)
    {
        const auto &boundaryFace = elemCtx.stencil(timeIdx).boundaryFace(bfIdx);
        const auto &normal = boundaryFace.normal();

        ///////////////
        // calculate the weights of the upstream and the downstream degrees of freedom
        ///////////////
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                this->filterVelocity_[phaseIdx] = 0;
                this->volumeFlux_[phaseIdx] = 0;
                continue;
            }

            calculateForchheimerVelocity_(phaseIdx);
            this->volumeFlux_[phaseIdx] = (this->filterVelocity_[phaseIdx] * normal);
        }
    }

    void calculateForchheimerVelocity_(int phaseIdx)
    {
        // initial guess: filter velocity is zero
        DimVector &velocity = this->filterVelocity_[phaseIdx];
        velocity = 0;

        // the change of velocity between two consecutive Newton iterations
        DimVector deltaV(1e5);
        // the function value that is to be minimized of the equation that is to be
        // fulfilled
        DimVector residual;
        // derivative of equation that is to be solved
        DimMatrix gradResid;

        // search by means of the Newton method for a root of Forchheimer equation
        int newtonIter = 0;
        while (deltaV.two_norm() > 1e-11) {
            if (newtonIter >= 50)
                OPM_THROW(Opm::NumericalProblem,
                          "Could not determine Forchheimer velocity within "
                          << newtonIter << " iterations");
            ++newtonIter;

            // calculate the residual and its Jacobian matrix
            gradForchheimerResid_(residual, gradResid, phaseIdx);

            // newton method
            gradResid.solve(deltaV, residual);
            velocity -= deltaV;
        }
    }

    void forchheimerResid_(DimVector &residual, int phaseIdx) const
    {
        const DimVector &velocity = this->filterVelocity_[phaseIdx];

        // Obtaining the upstreamed quantities
        Scalar mobility = this->mobility_[phaseIdx];
        Scalar density = density_[phaseIdx];
        Scalar mobilityPassabilityRatio = mobilityPassabilityRatio_[phaseIdx];

        // optain the quantites for the integration point
        const auto &pGrad = this->potentialGrad_[phaseIdx];

        // residual = v_\alpha
        residual = velocity;

        // residual += mobility_\alpha K(\grad p_\alpha - \rho_\alpha g)
        this->K_.usmv(mobility, pGrad, residual);

        // Forchheimer turbulence correction:
        //
        // residual +=
        //   \rho_\alpha
        //   * mobility_\alpha
        //   * C_E / \eta_{r,\alpha}
        //   * abs(v_\alpha) * sqrt(K)*v_\alpha
        sqrtK_.usmv(density*mobilityPassabilityRatio*ergunCoefficient_*velocity.two_norm(),
                    velocity,
                    residual);

        Valgrind::CheckDefined(residual);
    }

    void gradForchheimerResid_(DimVector &residual, DimMatrix &gradResid,
                               int phaseIdx)
    {
        DimVector &velocity = this->filterVelocity_[phaseIdx];
        forchheimerResid_(residual, phaseIdx);

        Scalar eps = 1e-11;
        DimVector tmp;
        for (int i = 0; i < dimWorld; ++i) {
            Scalar coordEps = std::max(eps, velocity[i] * (1 + eps));
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
            for (int j = 0; j < dimWorld; j++) {
                if (i == j)
                    continue;

                if (std::abs(K[i][j]) > 1e-25)
                    return false;
            }
        }
        return true;
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

protected:
    // intrinsic permeability tensor and its square root
    DimMatrix sqrtK_;

    // Ergun coefficient of all phases at the integration point
    Scalar ergunCoefficient_;

    // Passability of all phases at the integration point
    Scalar mobilityPassabilityRatio_[numPhases];

    // Density of all phases at the integration point
    Scalar density_[numPhases];
};

} // namespace Ewoms

#endif
