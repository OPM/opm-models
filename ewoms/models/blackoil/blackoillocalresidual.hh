/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Ewoms::BlackOilLocalResidual
 */
#ifndef EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH
#define EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH

#include "blackoilproperties.hh"

namespace Ewoms {
/*!
 * \ingroup BlackOilModel
 *
 * \brief Calculates the local residual of the black oil model.
 */
template <class TypeTag>
class BlackOilLocalResidual : public GET_PROP_TYPE(TypeTag, DiscLocalResidual)
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

public:
    /*!
     * \copydoc FvBaseLocalResidual::computeStorage
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int dofIdx,
                        int timeIdx) const
    {
        // retrieve the intensive quantities for the SCV at the specified point in time
        const IntensiveQuantities &intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);

        storage = 0.0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                storage[conti0EqIdx + compIdx] +=
                    intQuants.porosity()
                    * intQuants.fluidState().saturation(phaseIdx)
                    * intQuants.fluidState().molarity(phaseIdx, compIdx);
                assert(std::isfinite(storage[conti0EqIdx + compIdx]));
            }
        }
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        //std::cout << "------------\n";
        const ExtensiveQuantities &extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const ExtensiveQuantities &evalPointExtQuants =
            elemCtx.evalPointExtensiveQuantities(scvfIdx, timeIdx);

        flux = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            int upIdx = evalPointExtQuants.upstreamIndex(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                const IntensiveQuantities &up = elemCtx.intensiveQuantities(upIdx, /*timeIdx=*/0);

                // add advective flux of current component in current phase
                flux[conti0EqIdx + compIdx] +=
                    extQuants.volumeFlux(phaseIdx) * up.fluidState().molarity(phaseIdx, compIdx);

                assert(std::isfinite(flux[conti0EqIdx + compIdx]));
            }
        }
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int dofIdx,
                       int timeIdx) const
    {
        // retrieve the source term intrinsic to the problem
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
    }
};

} // namespace Ewoms

#endif
