// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Melanie Darcis                                    *
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
 * \copydoc Dumux::ImmiscibleLocalResidual
 */
#ifndef DUMUX_IMMISCIBLE_LOCAL_RESIDUAL_BASE_HH
#define DUMUX_IMMISCIBLE_LOCAL_RESIDUAL_BASE_HH

#include "immiscibleproperties.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/boxmodels/common/boxmodel.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup ImmiscibleBoxModel
 *
 * \brief Calculates the local residual of the immiscible multi-phase
 *        box model.
 */
template<class TypeTag>
class ImmiscibleLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy)
    };

    typedef BoxMultiPhaseEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Adds the amount all conservation quantities (e.g. phase
     *        mass) within a single fluid phase
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::boxScvCtxParams
     * \copydetails Doxygen::phaseIdxParam
     */
    void addPhaseStorage(EqVector &storage,
                         const ElementContext &elemCtx,
                         int scvIdx,
                         int timeIdx,
                         int phaseIdx) const
    {
        // retrieve the volume variables for the SCV at the specified
        // point in time
        const VolumeVariables &volVars = elemCtx.volVars(scvIdx, timeIdx);

        storage[conti0EqIdx + phaseIdx] +=
            volVars.porosity()
            * volVars.fluidState().saturation(phaseIdx)
            * volVars.fluidState().density(phaseIdx);
        
        EnergyModule::addPhaseStorage(storage, elemCtx.volVars(scvIdx, timeIdx), phaseIdx);
    }

    /*!
     * \copydoc BoxLocalResidual::computeStorage
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
     * \copydoc BoxLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        flux = 0;
        asImp_()->addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        asImp_()->addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Add the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \copydetails computeFlux
     */
    void addAdvectiveFlux(RateVector &flux,
                          const ElementContext &elemCtx,
                          int scvfIdx,
                          int timeIdx) const
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const FluxVariables &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase. The upstream decision has to be
            // made using the evaluation point and *not* the current
            // solution. (although the actual secondary variables must
            // obviously come from the current solution.)
            int upIdx = evalPointFluxVars.upstreamIndex(phaseIdx);

            const VolumeVariables &up = elemCtx.volVars(upIdx, /*timeIdx=*/0);

            // add advective flux of current component in current
            // phase
            flux[conti0EqIdx + phaseIdx] +=
                fluxVars.volumeFlux(phaseIdx)
                * up.fluidState().density(phaseIdx);
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Adds the diffusive flux of all conservation quantitis to
     *        the flux vector over the face of a sub-control volume.
     *
     * For the immiscible model, this is a no-op for mass fluxes. For
     * energy it adds the contribution of heat conduction to the
     * enthalpy flux.
     *
     * \copydetails computeFlux
     */
    void addDiffusiveFlux(RateVector &flux,
                          const ElementContext &elemCtx,
                          int scvfIdx,
                          int timeIdx) const
    {
        // no diffusive mass fluxes for the immiscible model
        
        // heat conduction
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc BoxLocalResidual::computeSource
     *
     * By default, this method only asks the problem to specify a
     * source term.
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

protected:
    Implementation *asImp_()
    { return static_cast<Implementation *> (this);  }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *> (this); }
};

}

#endif
