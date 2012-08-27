// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Melanie Darcis                                    *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
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
 * \brief Element-wise calculation of the residual for the two-phase box model.
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
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase box model.
 *
 * This class is also used for the non-isothermal model, which means
 * that it uses static polymorphism.
 */
template<class TypeTag>
class ImmiscibleLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
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
     *  \param result The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
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
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a finite sub-control volume.
     *
     *  \param result The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        // retrieve the volume variables for the SCV at the specified
        // point in time
        const VolumeVariables &volVars = elemCtx.volVars(scvIdx, timeIdx);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            storage[conti0EqIdx + phaseIdx] =
                volVars.porosity()
                * volVars.fluidState().saturation(phaseIdx)
                * volVars.fluidState().density(phaseIdx);
        }

        EnergyModule::addSolidHeatStorage(storage, elemCtx.volVars(scvIdx, timeIdx));
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub-control
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face for each phase
     * \param faceIdx The index of the SCV face
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        asImp_()->computeDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each phase
     * \param fluxVars The flux variables at the current SCV
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
     */
    void computeAdvectiveFlux(RateVector &flux,
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
            int dnIdx = evalPointFluxVars.downstreamIndex(phaseIdx);

            const VolumeVariables &up = elemCtx.volVars(upIdx, /*timeIdx=*/0);
            const VolumeVariables &dn = elemCtx.volVars(dnIdx, /*timeIdx=*/0);

            // add advective flux of current component in current
            // phase
            flux[conti0EqIdx + phaseIdx] +=
                fluxVars.volumeFlux(phaseIdx)
                * (fluxVars.upstreamWeight(phaseIdx)
                   * up.fluidState().density(phaseIdx)
                   +
                   fluxVars.downstreamWeight(phaseIdx)
                   * dn.fluidState().density(phaseIdx));
        }

        EnergyModule::addAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each phase
     * \param fluxData The flux variables at the current SCV
     *
     * It doesn't do anything in two-phase model but is used by the
     * non-isothermal two-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    {
        // no diffusive mass fluxes for the immiscible model
        
        // heat conduction
        EnergyModule::addDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q The source/sink in the SCV for each phase
     * \param localVertexIdx The index of the SCV
     *
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
