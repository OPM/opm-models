// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2011 by Andreas Lauser                               *
 *   Copyright (C) 2011 by Philipp Nuske                                     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
#ifndef DUMUX_NCP_LOCAL_RESIDUAL_MASS_HH
#define DUMUX_NCP_LOCAL_RESIDUAL_MASS_HH

#include <dune/common/fvector.hh>

#include <dumux/boxmodels/ncp/ncpproperties.hh>
#include <dumux/material/constraintsolvers/compositionfromfugacities.hh>

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>

#include "../diffusion/ncpdiffusion.hh"
#include "../energy/ncplocalresidualenergy.hh"

namespace Dumux
{
/*!
 * \brief The mass conservation part of the compositional NCP model.
 *
 * This is the class represents methods which are shared amongst all
 * mass conservation modules.
 */
template<class TypeTag>
class NcpLocalResidualMassCommon
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef NcpDiffusion<TypeTag, enableDiffusion> Diffusion;

    typedef NcpLocalResidualEnergy<TypeTag, enableEnergy> EnergyResid;

public:
    /*!
     * \brief Evaluate the amount moles within a sub-control volume in
     *        a phase.
     *
     * The result should be averaged over the volume.
     */
    static void computePhaseStorage(ComponentVector &compStorage,
                                    const VolumeVariables &volVars,
                                    int phaseIdx)
    {
        // compute storage term of all components within all phases
        compStorage = 0;
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            compStorage[compIdx] +=
                volVars.fluidState().saturation(phaseIdx)*
                volVars.fluidState().molarity(phaseIdx, compIdx);
        }

        compStorage *= volVars.porosity();
    }

    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     */
    static void computeAdvectiveFlux(ComponentVector &compFlux,
                                     const ElementContext &elemCtx,
                                     int scvfIdx,
                                     int timeIdx,
                                     int phaseIdx)
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const FluxVariables &evalFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);
        const VolumeVariables &up = elemCtx.volVars(evalFluxVars.upstreamIdx(phaseIdx), timeIdx);

        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            compFlux[compIdx] =
                up.fluidState().molarity(phaseIdx, compIdx)
                * fluxVars.filterVelocityNormal(phaseIdx);
        }
    }

    /*!
     * \brief Evaluates the advective flux of all conservation
     *        quantities over a face of a subcontrol volume via a
     *        fluid phase.
     */
    static void computeDiffusiveFlux(ComponentVector &flux,
                                     const ElementContext &elemCtx,
                                     int scvfIdx,
                                     int timeIdx,
                                     int phaseIdx)
    {
        if (!enableDiffusion) {
            flux = 0.0;
            return;
        }

        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const VolumeVariables &volVarsI = elemCtx.volVars(fluxVars.insideIdx(), timeIdx);
        const VolumeVariables &volVarsJ = elemCtx.volVars(fluxVars.outsideIdx(), timeIdx);
        if (volVarsI.fluidState().saturation(phaseIdx) < 1e-4 ||
            volVarsJ.fluidState().saturation(phaseIdx) < 1e-4)
        {
            return; // phase is not present in one of the finite volumes
        }

        // approximate the total concentration of the phase at the
        // integration point by the arithmetic mean of the
        // concentration of the sub-control volumes
        Scalar molarDensity;
        molarDensity = volVarsI.fluidState().molarDensity(phaseIdx);
        molarDensity += volVarsJ.fluidState().molarDensity(phaseIdx);
        molarDensity /= 2;

        Diffusion::flux(flux, elemCtx, scvfIdx, timeIdx, phaseIdx, molarDensity);
    }
};

/*!
 * \brief The mass conservation part of the compositional NCP model.
 *
 * This is the specialization for the case where kinetic mass transfer
 * is not considered.
 */
template<class TypeTag>
class NcpLocalResidualMass
{
    typedef NcpLocalResidualMassCommon<TypeTag> MassCommon;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { conti0EqIdx = Indices::conti0EqIdx };

    typedef typename Dune::FieldVector<Scalar, numComponents> ComponentVector;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef NcpLocalResidualEnergy<TypeTag, enableEnergy> EnergyResid;

public:
    /*!
     * \brief Calculate the storage for all mass balance equations
     */
    static void computeStorage(EqVector &storage,
                               const VolumeVariables &volVars)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] = 0.0;

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            addPhaseStorage(storage, volVars, phaseIdx);
        }
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     *        within a single fluid phase
     */
    static void addPhaseStorage(EqVector &storage,
                                const VolumeVariables &volVars,
                                int phaseIdx)
    {
        // calculate the component-wise mass storage
        ComponentVector phaseComponentValues;
        MassCommon::computePhaseStorage(phaseComponentValues,
                                        volVars,
                                        phaseIdx);

        // copy to the primary variables
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            storage[conti0EqIdx + compIdx] += phaseComponentValues[compIdx];
    }

    /*!
     * \brief Calculate the storage for all mass balance equations
     */
    static void computeFlux(EqVector &flux,
                            const ElementContext &elemCtx,
                            int scvfIdx,
                            int timeIdx)
    {
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            flux[conti0EqIdx + compIdx] = 0.0;

        ComponentVector diffusiveFlux(0.);
        ComponentVector advectiveFlux(0.);
        ComponentVector totalFluxes[numPhases]; // goes into the energy module

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            MassCommon::computeAdvectiveFlux(advectiveFlux,
                                             elemCtx,
                                             scvfIdx,
                                             timeIdx,
                                             phaseIdx);
            MassCommon::computeDiffusiveFlux(diffusiveFlux,
                                             elemCtx,
                                             scvfIdx,
                                             timeIdx,
                                             phaseIdx);

            // add the fluxes to the primary variables vector
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                flux[conti0EqIdx + compIdx] +=
                    advectiveFlux[compIdx] + diffusiveFlux[compIdx];

            // Right now I think that adding the two contributions
            // individually into the flux is best for debugging and
            // understanding. The Energy module needs both
            // contributions.
            totalFluxes[phaseIdx] = advectiveFlux + diffusiveFlux;

            Valgrind::CheckDefined(advectiveFlux);
            Valgrind::CheckDefined(diffusiveFlux);
            Valgrind::CheckDefined(flux);
        }

        // \todo
        //
        // The computeFlux() of the energy module needs a
        // component-wise flux for the enthalpy transport. It makes
        // some sense calling into the energy module from here,
        // because energy is carried by mass. However, it is not
        // really a clean solution.

        // energy transport in fluid phases
        EnergyResid::computeFlux(flux,
                                 elemCtx,
                                 scvfIdx,
                                 timeIdx,
                                 totalFluxes);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Calculate the source terms for all mass balance
     *        equations
     */
    static void computeSource(EqVector &source,
                              const ElementContext &elemCtx,
                              int scvIdx,
                              int timeIdx)
    {
        // set source terms for mass to 0 for all components
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            source[conti0EqIdx + compIdx] = 0.0;
        Valgrind::CheckDefined(source);
    }
};

} // end namepace

#endif
