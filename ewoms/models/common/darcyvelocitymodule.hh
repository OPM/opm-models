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
 *        Darcy relation.
 */
#ifndef EWOMS_DARCY_MODULE_HH
#define EWOMS_DARCY_MODULE_HH

#include <ewoms/disc/common/fvbaseproperties.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(MaterialLaw);
}}

namespace Ewoms {
template <class TypeTag>
class DarcyVolumeVariables;

template <class TypeTag>
class DarcyFluxVariables;

template <class TypeTag>
class DarcyBaseProblem;

/*!
 * \ingroup DarcyVelocity
 * \brief Specifies a velocity module which uses the Darcy relation.
 */
template <class TypeTag>
struct DarcyVelocityModule
{
    typedef DarcyVolumeVariables<TypeTag> VelocityVolumeVariables;
    typedef DarcyFluxVariables<TypeTag> VelocityFluxVariables;
    typedef DarcyBaseProblem<TypeTag> VelocityBaseProblem;

    /*!
     * \brief Register all run-time parameters for the velocity module.
     */
    static void registerParameters()
    { }
};

/*!
 * \ingroup DarcyVelocity
 * \brief Provides the defaults for the parameters required by the
 *        Darcy velocity approach.
 */
template <class TypeTag>
class DarcyBaseProblem
{ };

/*!
 * \ingroup DarcyVelocity
 * \brief Provides the volume variables for the Darcy velocity module
 */
template <class TypeTag>
class DarcyVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
protected:
    void update_(const ElementContext &elemCtx, int dofIdx, int timeIdx)
    { }
};

/*!
 * \ingroup DarcyVelocity
 * \brief Provides the Darcy velocity module
 *
 * The commonly used Darcy relation looses its validity for Reynolds
 * numbers \f$ Re > 1\f$.  If one encounters flow velocities in porous
 * media above this threshold, the Darcy relation can be
 * used. Like the Darcy relation, it relates the gradient in potential
 * to velocity.  However, this relation is not linear (as in the Darcy
 * case) any more.
 *
 * Therefore, the Newton scheme is used to solve the non-linear
 * relation. This velocity is then used like the Darcy velocity
 * e.g. by the local residual.
 *
 * For Reynolds numbers above \f$\approx 500\f$ the standard Darcy
 * relation also looses it's validity.
 *
 * The Darcy equation is given by the following relation:
 *
 * \f[
  \nabla p_\alpha - \rho_\alpha \vec{g} =
  - \frac{\mu_\alpha}{k_{r,\alpha} K}\vec{v}_\alpha
  - \frac{\rho_\alpha C_E}{\eta_{r,\alpha} \sqrt{K}}
  \left| \vec{v}_\alpha \right| \vec{v}_\alpha
 \f]
 *
 * Where \f$C_E\f$ is the modified Ergun parameter and
 * \f$\eta_{r,\alpha}\f$ is the passability which is given by a
 * closure relation.
 */
template <class TypeTag>
class DarcyFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

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

        const auto &volVarsI = elemCtx.volVars(asImp_().interiorIndex(), timeIdx);
        const auto &volVarsJ = elemCtx.volVars(asImp_().exteriorIndex(), timeIdx);

        // calculate the intrinsic permeability
        const auto &Ki = volVarsI.intrinsicPermeability();
        const auto &Kj = volVarsJ.intrinsicPermeability();
        problem.meanK(K_, Ki, Kj);
        Valgrind::CheckDefined(K_);

        const auto &scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const DimVector &normal = scvf.normal();
        Valgrind::CheckDefined(normal);

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
            mobility_[phaseIdx] = up.mobility(phaseIdx);

            calculateDarcyVelocity_(phaseIdx);
            volumeFlux_[phaseIdx] = (filterVelocity_[phaseIdx] * normal);
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
        const auto &elemCtx = context.elementContext();
        const auto &problem = elemCtx.problem();
        int interiorScvIdx = asImp_().interiorIndex();
        const auto &volVarsInterior = elemCtx.volVars(interiorScvIdx, timeIdx);
        const auto &fsInterior = volVarsInterior.fluidState();

        // calculate the intrinsic permeability
        K_ = volVarsInterior.intrinsicPermeability();

        const auto &boundaryFace = context.stencil(timeIdx).boundaryFace(bfIdx);
        const auto &normal = boundaryFace.normal();

        const auto &matParams = problem.materialLawParams(elemCtx, interiorScvIdx, timeIdx);

        Scalar kr[numPhases];
        MaterialLaw::relativePermeabilities(kr, matParams, fluidState);

        ///////////////
        // calculate the weights of the upstream and the downstream
        // control volumes
        ///////////////
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            if (!context.model().phaseIsConsidered(phaseIdx)) {
                filterVelocity_[phaseIdx] = 0;
                volumeFlux_[phaseIdx] = 0;
                continue;
            }

            // calculate the actual darcy velocities by multiplying
            // the current "filter velocity" with the upstream mobility
            if (fsInterior.pressure(phaseIdx) < fluidState.pressure(phaseIdx)) {
                // the exterior of the domain has higher pressure. we
                // need to calculate the mobility
                mobility_[phaseIdx] =
                    kr[phaseIdx]
                    / FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            }
            else {
                mobility_[phaseIdx] = volVarsInterior.mobility(phaseIdx);
            }

            calculateDarcyVelocity_(phaseIdx);

            Scalar tmp = 0;
            for (int dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                tmp += filterVelocity_[phaseIdx][dimIdx] * normal[dimIdx];

            volumeFlux_[phaseIdx] = tmp;
        }
    }

    void calculateDarcyVelocity_(int phaseIdx)
    {
        // calculate the "prelimanary" filter velocity not
        // taking the mobility into account
        K_.mv(asImp_().potentialGrad(phaseIdx), filterVelocity_[phaseIdx]);
        filterVelocity_[phaseIdx] *= - mobility_[phaseIdx];
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

protected:
    // intrinsic permeability tensor and its square root
    DimMatrix K_;

    // filter velocities of all phases [m/s]
    Scalar mobility_[numPhases];

    // filter velocities of all phases [m/s]
    DimVector filterVelocity_[numPhases];

    // the volumetric flux of all fluid phases over the control
    // volume's face [m^3/s]
    Scalar volumeFlux_[numPhases];
};

} // namespace Ewoms

#endif
