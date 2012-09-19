// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
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
 * \copydoc Dumux::PvsLocalResidual
 */
#ifndef DUMUX_PVS_LOCAL_RESIDUAL_HH
#define DUMUX_PVS_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/common/math.hh>

#include "pvsproperties.hh"
#include "pvsvolumevariables.hh"
#include "pvsfluxvariables.hh"
#include "pvsnewtoncontroller.hh"

namespace Dumux {

/*!
 * \ingroup PvsModel
 *
 * \brief Element-wise calculation of the local residual for the
 *        compositional multi-phase primary variable switching box model.
 */
template<class TypeTag>
class PvsLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
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
     * \copydoc ImmiscibleLocalResidual::addPhaseStorage
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
     * \copydoc ImmiscibleLocalResidual::computeStorage
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
     * \copydoc ImmiscibleLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        flux = 0.0;
        addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);

        addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::addAdvectiveFlux
     */
    void addAdvectiveFlux(RateVector &flux,
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
     * \copydoc ImmiscibleLocalResidual::addDiffusiveFlux
     */
    void addDiffusiveFlux(RateVector &flux,
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
     * \copydoc ImmiscibleLocalResidual::computeSource
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
