// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012-2013 by Andreas Lauser                               *
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
 * \brief Classes required for molecular diffusion.
 */
#ifndef EWOMS_VCFV_DIFFUSION_MODULE_HH
#define EWOMS_VCFV_DIFFUSION_MODULE_HH

#include <ewoms/disc/vcfv/vcfvproperties.hh>
#include <dune/common/fvector.hh>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Indices);
}}

namespace Ewoms {
/*!
 * \ingroup VcfvDiffusion
 * \class Ewoms::VcfvDiffusionModule
 * \brief Provides the auxiliary methods required for consideration of the diffusion equation.
 */
template <class TypeTag, bool enableDiffusion>
class VcfvDiffusionModule;

/*!
 * \copydoc Ewoms::VcfvDiffusionModule
 */
template <class TypeTag>
class VcfvDiffusionModule<TypeTag, /*enableDiffusion=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    /*!
     * \brief Register all run-time parameters for the diffusion module.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Adds the diffusive mass flux flux to the flux vector
     *        over the face of a sub-control volume.
      */
    template <class Context>
    static void addDiffusiveFlux(RateVector &flux,
                                 const Context &context,
                                 int spaceIdx,
                                 int timeIdx)
    { }
};

/*!
 * \copydoc Ewoms::VcfvDiffusionModule
 */
template <class TypeTag>
class VcfvDiffusionModule<TypeTag, /*enableDiffusion=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { conti0EqIdx = Indices::conti0EqIdx };

public:
    /*!
     * \brief Register all run-time parameters for the diffusion module.
     */
    static void registerParameters()
    { }

    /*!
     * \brief Adds the mass flux due to molecular diffusion to the
     *        flux vector over the face of a sub-control volume.
     */
    template <class Context>
    static void addDiffusiveFlux(RateVector &flux,
                                 const Context &context,
                                 int spaceIdx,
                                 int timeIdx)
    {
        const auto &fluxVars = context.fluxVars(spaceIdx, timeIdx);

        const auto &fluidStateI = context.volVars(fluxVars.insideIndex(), timeIdx).fluidState();
        const auto &fluidStateJ = context.volVars(fluxVars.outsideIndex(), timeIdx).fluidState();

        RateVector molarRate(0.0);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // arithmetic mean of the phase's molar density
            Scalar rhoMolar =
                (fluidStateI.molarDensity(phaseIdx)
                 + fluidStateJ.molarDensity(phaseIdx))/2;

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                // mass flux due to molecular diffusion
                molarRate[conti0EqIdx + compIdx] +=
                    - rhoMolar
                    * fluxVars.moleFractionGradientNormal(phaseIdx, compIdx)
                    * fluxVars.effectiveDiffusionCoefficient(phaseIdx, compIdx);
        }

        flux += molarRate;
    }
};

/*!
 * \ingroup VcfvDiffusion
 * \class Ewoms::VcfvDiffusionVolumeVariables
 *
 * \brief Provides the volumetric quantities required for the
 *        calculation of molecular diffusive fluxes.
 */
template <class TypeTag, bool enableDiffusion>
class VcfvDiffusionVolumeVariables;

/*!
 * \copydoc Ewoms::VcfvDiffusionVolumeVariables
 */
template <class TypeTag>
class VcfvDiffusionVolumeVariables<TypeTag, /*enableDiffusion=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    /*!
     * \brief Returns the tortuousity of the sub-domain of a fluid
     *        phase in the porous medium.
     */
    Scalar tortuosity(int phaseIdx) const
    { OPM_THROW(std::logic_error, "Method tortuosity() does not make sense if diffusion is disabled"); }

    /*!
     * \brief Returns the molecular diffusion coefficient for a
     *        component in a phase.
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    { OPM_THROW(std::logic_error, "Method diffusionCoefficient() does not make sense if diffusion is disabled"); }

    /*!
     * \brief Returns the effective molecular diffusion coefficient of
     *        the porous medium for a component in a phase.
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIdx) const
    { OPM_THROW(std::logic_error, "Method effectiveDiffusionCoefficient() does not make sense if diffusion is disabled"); }

protected:
    /*!
     * \brief Update the quantities required to calculate diffusive
     *        mass fluxes.
     */
    template <class FluidState>
    void update_(FluidState &fs,
                 ParameterCache &paramCache,
                 const ElementContext &elemCtx,
                 int scvIdx,
                 int timeIdx)
    { }
};

/*!
 * \copydoc Ewoms::VcfvDiffusionVolumeVariables
 */
template <class TypeTag>
class VcfvDiffusionVolumeVariables<TypeTag, /*enableDiffusion=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename FluidSystem::ParameterCache ParameterCache;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

public:
    /*!
     * \brief Returns the molecular diffusion coefficient for a
     *        component in a phase.
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    { return diffusionCoefficient_[phaseIdx][compIdx]; }

    /*!
     * \brief Returns the tortuousity of the sub-domain of a fluid
     *        phase in the porous medium.
     */
    Scalar tortuosity(int phaseIdx) const
    { return tortuosity_[phaseIdx]; }

    /*!
     * \brief Returns the effective molecular diffusion coefficient of
     *        the porous medium for a component in a phase.
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIdx) const
    { return tortuosity_[phaseIdx]*diffusionCoefficient_[phaseIdx][compIdx]; }

protected:
    /*!
     * \brief Update the quantities required to calculate diffusive
     *        mass fluxes.
     */
    template <class FluidState>
    void update_(FluidState &fluidState,
                 ParameterCache &paramCache,
                 const ElementContext &elemCtx,
                 int scvIdx,
                 int timeIdx)
    {
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            // TODO: let the problem do this (this is a constitutive
            // relation of which the model should be free of from the
            // abstraction POV!)
            tortuosity_[phaseIdx] =
                1.0/(volVars.porosity() * volVars.porosity())
                * std::pow(std::max(0.0001, volVars.porosity() * volVars.fluidState().saturation(phaseIdx)), 7.0/3);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                diffusionCoefficient_[phaseIdx][compIdx] =
                    FluidSystem::diffusionCoefficient(fluidState, paramCache, phaseIdx, compIdx);
            }
        }
    }

private:
    Scalar tortuosity_[numPhases];
    Scalar diffusionCoefficient_[numPhases][numComponents];
};

/*!
 * \ingroup VcfvDiffusion
 * \class Ewoms::VcfvDiffusionFluxVariables
 *
 * \brief Provides the quantities required to calculate diffusive mass fluxes.
 */
template <class TypeTag, bool enableDiffusion>
class VcfvDiffusionFluxVariables;

/*!
 * \copydoc Ewoms::VcfvDiffusionFluxVariables
 */
template <class TypeTag>
class VcfvDiffusionFluxVariables<TypeTag, /*enableDiffusion=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

protected:
    /*!
     * \brief Update the quantities required to calculate
     *        the diffusive mass fluxes.
     */
    void update_(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    { }

    template <class Context, class FluidState>
    void updateBoundary_(const Context &context,
                         int bfIdx,
                         int timeIdx,
                         const FluidState &fluidState)
    { }

public:
    /*!
     * \brief The the gradient of the mole fraction times the face normal.
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    Scalar moleFractionGradientNormal(int phaseIdx, int compIdx) const
    { OPM_THROW(std::logic_error, "Method moleFractionGradient() does not make sense if diffusion is disabled."); }

    /*!
     * \brief The effective diffusion coeffcient of a component in a
     *        fluid phase at the face's integration point
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIdx) const
    { OPM_THROW(std::logic_error, "Method effectiveDiffusionCoefficient() does not make sense if diffusion is disabled."); }
};

/*!
 * \copydoc Ewoms::VcfvDiffusionFluxVariables
 */
template <class TypeTag>
class VcfvDiffusionFluxVariables<TypeTag, /*enableDiffusion=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

protected:
    /*!
     * \brief Update the quantities required to calculate
     *        the diffusive mass fluxes.
     */
    void update_(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        const auto &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);

        const auto &volVarsInside = elemCtx.volVars(fluxVars.insideIndex(), timeIdx);
        const auto &volVarsOutside = elemCtx.volVars(fluxVars.outsideIndex(), timeIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                DimVector moleFractionGradient(0.0);
                DimVector tmp;

                // calculate mole fraction gradient using finite element
                // gradients
                for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++)
                {
                    const auto &feGrad = scvf.grad[scvIdx];
                    const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);

                    tmp = feGrad;
                    tmp *= volVars.fluidState().moleFraction(phaseIdx, compIdx);
                    moleFractionGradient += tmp;
                }

                moleFractionGradientNormal_[phaseIdx][compIdx] =
                    scvf.normal*moleFractionGradient;
                Valgrind::CheckDefined(moleFractionGradientNormal_[phaseIdx][compIdx]);

                // use the arithmetic average for the effective
                // diffusion coefficients.
                effectiveDiffusionCoefficient_[phaseIdx][compIdx] =
                    (volVarsInside.effectiveDiffusionCoefficient(phaseIdx, compIdx)
                     + volVarsOutside.effectiveDiffusionCoefficient(phaseIdx, compIdx))/2;

                Valgrind::CheckDefined(effectiveDiffusionCoefficient_[phaseIdx][compIdx]);
            }
        }
    }

    template <class Context, class FluidState>
    void updateBoundary_(const Context &context,
                         int bfIdx,
                         int timeIdx,
                         const FluidState &fluidState)
    {
        const auto &fvElemGeom = context.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.boundaryFace[bfIdx];

        const auto &elemCtx = context.elemContext();
        int insideScvIdx = scvf.i;
        const auto &insideScv = elemCtx.fvElemGeom(timeIdx).subContVol[insideScvIdx];

        const auto &volVarsInside = elemCtx.volVars(insideScvIdx, timeIdx);
        const auto &fluidStateInside = volVarsInside.fluidState();

        // distance between the center of the SCV and center of the boundary face
        DimVector distVec = scvf.ipGlobal;
        distVec -= context.element().geometry().global(insideScv.localGeometry->center());

        const DimVector &n = scvf.normal;
        Scalar dist = (distVec * n)/n.two_norm();

        // if the following assertation triggers, the center of the
        // center of the interior SCV was not inside the element!
        assert(dist > 0);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                DimVector moleFractionGradient(0.0);

                // calculate mole fraction gradient using two-point
                // gradients
                moleFractionGradientNormal_[phaseIdx][compIdx] =
                    (fluidState.moleFraction(phaseIdx, compIdx) -
                     fluidStateInside.moleFraction(phaseIdx, compIdx))/dist;
                Valgrind::CheckDefined(moleFractionGradientNormal_[phaseIdx][compIdx]);

                // use effective diffusion coefficients of the interior finite volume.
                effectiveDiffusionCoefficient_[phaseIdx][compIdx] =
                    volVarsInside.effectiveDiffusionCoefficient(phaseIdx, compIdx);

                Valgrind::CheckDefined(effectiveDiffusionCoefficient_[phaseIdx][compIdx]);
            }
        }
    }

public:
    /*!
     * \brief The the gradient of the mole fraction times the face normal.
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    Scalar moleFractionGradientNormal(int phaseIdx, int compIdx) const
    { return moleFractionGradientNormal_[phaseIdx][compIdx]; }

    /*!
     * \brief The effective diffusion coeffcient of a component in a
     *        fluid phase at the face's integration point
     *
     * \copydoc Doxygen::phaseIdxParam
     * \copydoc Doxygen::compIdxParam
     */
    Scalar effectiveDiffusionCoefficient(int phaseIdx, int compIdx) const
    { return effectiveDiffusionCoefficient_[phaseIdx][compIdx]; }

private:
    Scalar moleFractionGradientNormal_[numPhases][numComponents];
    Scalar effectiveDiffusionCoefficient_[numPhases][numComponents];
};

} // namespace Ewoms

#endif
