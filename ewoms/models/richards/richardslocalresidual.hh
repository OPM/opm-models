/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \copydoc Ewoms::RichardsLocalResidual
 */
#ifndef EWOMS_RICHARDS_LOCAL_RESIDUAL_HH
#define EWOMS_RICHARDS_LOCAL_RESIDUAL_HH

#include "richardsintensivequantities.hh"

#include "richardsextensivequantities.hh"

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 * \brief Element-wise calculation of the residual for the Richards model.
 */
template <class TypeTag>
class RichardsLocalResidual : public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { contiEqIdx = Indices::contiEqIdx };
    enum { liquidPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };

public:
    /*!
     * \copydoc ImmiscibleLocalResidual::computeStorage
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int dofIdx,
                        int timeIdx) const
    {
        const IntensiveQuantities &intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);

        // partial time derivative of the wetting phase mass
        storage[contiEqIdx] =
            intQuants.fluidState().density(liquidPhaseIdx)
            * intQuants.fluidState().saturation(liquidPhaseIdx)
            * intQuants.porosity();
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux, const ElementContext &elemCtx,
                     int scvfIdx, int timeIdx) const
    {
        const auto &extQuantsEval = elemCtx.evalPointExtensiveQuantities(scvfIdx, timeIdx);
        // const auto &extQuantsEval = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const auto &extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        // data attached to upstream and the downstream DOFs
        // of the current phase
        const IntensiveQuantities &up =
            elemCtx.intensiveQuantities(extQuantsEval.upstreamIndex(liquidPhaseIdx), timeIdx);
        const IntensiveQuantities &dn =
            elemCtx.intensiveQuantities(extQuantsEval.downstreamIndex(liquidPhaseIdx), timeIdx);

        flux[contiEqIdx] =
            extQuants.volumeFlux(liquidPhaseIdx)
            * (extQuants.upstreamWeight(liquidPhaseIdx)
               * up.fluidState().density(liquidPhaseIdx)
               +
               extQuants.downstreamWeight(liquidPhaseIdx)
               * dn.fluidState().density(liquidPhaseIdx));
    }

    /*!
     * \copydoc ImmiscibleLocalResidual::computeSource
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int dofIdx,
                       int timeIdx) const
    {
        elemCtx.problem().source(source,
                                 elemCtx,
                                 dofIdx,
                                 timeIdx);
    }
};

} // namespace Ewoms

#endif
