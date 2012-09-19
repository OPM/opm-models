// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \copydoc Dumux::BlackOilLocalResidual
 */
#ifndef DUMUX_BLACK_OIL_LOCAL_RESIDUAL_HH
#define DUMUX_BLACK_OIL_LOCAL_RESIDUAL_HH

#include "blackoilproperties.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

/*!
 * \ingroup BlackOilBoxModel
 *
 * \brief Calculates the local residual of the black oil
 *        box model.
 */
template<class TypeTag>
class BlackOilLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

public:
    /*!
     * \copydoc BoxLocalResidual::computeStorage
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        // retrieve the volume variables for the SCV at the specified
        // point in time
        const VolumeVariables &volVars = elemCtx.volVars(scvIdx, timeIdx);

        storage = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                storage[conti0EqIdx + compIdx] +=
                    volVars.porosity()
                    * volVars.fluidState().saturation(phaseIdx)
                    * volVars.fluidState().molarity(phaseIdx, compIdx);
            }
        }
    }

    /*!
     * \copydoc BoxLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const FluxVariables &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);

        flux = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int upIdx = evalPointFluxVars.upstreamIndex(phaseIdx);
                const VolumeVariables &up = elemCtx.volVars(upIdx, /*timeIdx=*/0);

                // add advective flux of current component in current
                // phase
                flux[conti0EqIdx + compIdx] +=
                    fluxVars.volumeFlux(phaseIdx)
                    * up.fluidState().molarity(phaseIdx, compIdx);
            }
        }
    }

    /*!
     * \copydoc BoxLocalResidual::computeSource
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        // retrieve the source term intrinsic to the problem
        elemCtx.problem().source(source, elemCtx, scvIdx, timeIdx);
    }
};

}

#endif
