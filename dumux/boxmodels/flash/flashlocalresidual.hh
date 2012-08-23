// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2011-2012 by Klaus Mosthaf                                *
 *   Copyright (C) 2010-2011 by Melanie Darcis                               *
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
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component box model.
 */
#ifndef DUMUX_FLASH_LOCAL_RESIDUAL_HH
#define DUMUX_FLASH_LOCAL_RESIDUAL_HH

#include "flashproperties.hh"
#include "flashvolumevariables.hh"
#include "flashfluxvariables.hh"

#include <dumux/boxmodels/modules/energy/multiphaseenergymodule.hh>
#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

namespace Dumux {
/*!
 * \ingroup FlashModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase two-component box model.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the 2P-2C flow.
 */
template<class TypeTag>
class FlashLocalResidual: public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
        conti0EqIdx = Indices::conti0EqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum  { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef BoxMultiPhaseEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Evaluate the storage term [kg/m^3] in a single phase.
     *
     * \param element The element
     * \param phaseIdx The index of the fluid phase
     */
    void addPhaseStorage(EqVector &storage,
                         const ElementContext &elemCtx,
                         int scvIdx,
                         int timeIdx,
                         int phaseIdx) const
    {
        const VolumeVariables &volVars = elemCtx.volVars(scvIdx, timeIdx);
        const auto &fs = volVars.fluidState();
        
        // compute storage term of all components within all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            int eqIdx = conti0EqIdx + compIdx;
            storage[eqIdx] +=
                fs.molarity(phaseIdx, compIdx)
                * fs.saturation(phaseIdx)
                * volVars.porosity();
        }

        EnergyModule::addPhaseStorage(storage, elemCtx.volVars(scvIdx, timeIdx), phaseIdx);
    }

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param result The number of moles of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        storage = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            addPhaseStorage(storage, elemCtx, scvIdx, timeIdx, phaseIdx);

        EnergyModule::addSolidHeatStorage(storage, elemCtx.volVars(scvIdx, timeIdx));
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        flux = 0.0;
        computeAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);

        computeDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    {
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VolumeVariables &up = elemCtx.volVars(evalPointFluxVars.upstreamIndex(phaseIdx), timeIdx);
            const VolumeVariables &dn = elemCtx.volVars(evalPointFluxVars.downstreamIndex(phaseIdx), timeIdx);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = conti0EqIdx + compIdx;

                flux[eqIdx] +=
                    fluxVars.volumeFlux(phaseIdx)
                    *(fluxVars.upstreamWeight(phaseIdx)
                      * up.fluidState().molarity(phaseIdx, compIdx)
                      +
                      fluxVars.downstreamWeight(phaseIdx)
                      * dn.fluidState().molarity(phaseIdx, compIdx));

                Valgrind::CheckDefined(flux[eqIdx]);
            }
        }
        
        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param vars The flux variables at the current sub control volume face
     */
    void computeDiffusiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    {
#if 0
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;

        // add diffusive flux of gas component in liquid phase
        Scalar tmp = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            int compIdx = 1; // HACK
            tmp = - (fluxVars.moleFracGrad(phaseIdx, compIdx)*normal);
            tmp *=
                fluxVars.porousDiffCoeff(phaseIdx, compIdx) *
                fluxVars.molarDensity(phaseIdx);

            flux[conti0EqIdx + compIdx] += tmp;
            flux[conti0EqIdx + (1 - compIdx)] -= tmp;
        }
#endif

        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source The source/sink term [kg/m^3] in the sub control
     *               volume for each component
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, scvIdx, timeIdx);
        Valgrind::CheckDefined(source);
    }
};

} // end namepace

#endif
