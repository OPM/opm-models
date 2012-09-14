// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2010 by Felix Bode                                        *
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
 * \brief Element-wise calculation of the residual for the Richards box model.
 */
#ifndef DUMUX_RICHARDS_LOCAL_RESIDUAL_HH
#define DUMUX_RICHARDS_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxlocalresidual.hh>

#include "richardsvolumevariables.hh"

#include "richardsfluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsModel
 * \brief Element-wise calculation of the residual for the Richards box model.
 */
template<class TypeTag>
class RichardsLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        contiEqIdx = Indices::contiEqIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

public:
    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the Richards
     *        model.
     *
     * This function should not include the source and sink terms.
     *
     * \param result Stores the average mass per unit volume for each phase \f$\mathrm{[kg/m^3]}\f$
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        const VolumeVariables &volVars = elemCtx.volVars(scvIdx, timeIdx);

        // partial time derivative of the wetting phase mass
        storage[contiEqIdx] =
            volVars.fluidState().density(wPhaseIdx)
            * volVars.fluidState().saturation(wPhaseIdx)
            * volVars.porosity();
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        const auto &fluxVarsEval = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);
        //const auto &fluxVarsEval = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);

        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up = elemCtx.volVars(fluxVarsEval.upstreamIndex(wPhaseIdx), timeIdx);
        const VolumeVariables &dn = elemCtx.volVars(fluxVarsEval.downstreamIndex(wPhaseIdx), timeIdx);

        flux[contiEqIdx] =
            fluxVars.volumeFlux(wPhaseIdx)
            *( fluxVars.upstreamWeight(wPhaseIdx)*up.fluidState().density(wPhaseIdx)
               +
               fluxVars.downstreamWeight(wPhaseIdx)*dn.fluidState().density(wPhaseIdx));
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        elemCtx.problem().source(source,
                                 elemCtx,
                                 scvIdx,
                                 timeIdx);
    }
};

}

#endif
