// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    typedef Opm::MathToolbox<Evaluation> Toolbox;

public:
    /*!
     * \copydoc FvBaseLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq> &storage,
                        const ElementContext &elemCtx,
                        int dofIdx,
                        int timeIdx) const
    {
        // retrieve the intensive quantities for the SCV at the specified point in time
        const IntensiveQuantities &intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);

        LhsEval waterSurfaceVolume =
            Toolbox::template toLhs<LhsEval>(intQuants.fluidState().saturation(waterPhaseIdx))
            * Toolbox::template toLhs<LhsEval>(intQuants.fluidState().invB(waterPhaseIdx))
            * Toolbox::template toLhs<LhsEval>(intQuants.porosity());
        LhsEval oilSurfaceVolume =
            Toolbox::template toLhs<LhsEval>(intQuants.fluidState().saturation(oilPhaseIdx))
            * Toolbox::template toLhs<LhsEval>(intQuants.fluidState().invB(oilPhaseIdx))
            * Toolbox::template toLhs<LhsEval>(intQuants.porosity());
        LhsEval gasSurfaceVolume =
            Toolbox::template toLhs<LhsEval>(intQuants.fluidState().saturation(gasPhaseIdx))
            * Toolbox::template toLhs<LhsEval>(intQuants.fluidState().invB(gasPhaseIdx))
            * Toolbox::template toLhs<LhsEval>(intQuants.porosity());

        storage[conti0EqIdx + waterCompIdx] = waterSurfaceVolume;
        storage[conti0EqIdx + oilCompIdx] = oilSurfaceVolume;
        storage[conti0EqIdx + gasCompIdx] = gasSurfaceVolume;

        // account for dissolved gas and vaporized oil
        if (FluidSystem::enableDissolvedGas()) {
            storage[conti0EqIdx + gasCompIdx] +=
                Toolbox::template toLhs<LhsEval>(intQuants.fluidState().Rs())
                * oilSurfaceVolume;
        }

        if (FluidSystem::enableVaporizedOil()) {
            storage[conti0EqIdx + oilCompIdx] +=
                Toolbox::template toLhs<LhsEval>(intQuants.fluidState().Rv())
                * gasSurfaceVolume;
        }

        // convert surface volumes to component masses
        unsigned pvtRegionIdx = intQuants.pvtRegionIndex();
        storage[conti0EqIdx + waterCompIdx] *=
            FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
        storage[conti0EqIdx + gasCompIdx] *=
            FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
        storage[conti0EqIdx + oilCompIdx] *=
            FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        assert(timeIdx == 0);

        const ExtensiveQuantities &extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        for (int eqIdx=0; eqIdx < numEq; eqIdx++)
            flux[eqIdx] = 0;

        Evaluation b[numPhases];
        unsigned interiorIdx = extQuants.interiorIndex();

        unsigned upIdx[numPhases];
        const IntensiveQuantities* up[numPhases];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)  {
            upIdx[phaseIdx] = extQuants.upstreamIndex(phaseIdx);
            up[phaseIdx] = &elemCtx.intensiveQuantities(upIdx[phaseIdx], /*timeIdx=*/0);

            // this is a bit hacky because it is specific to the element-centered
            // finite volume scheme. (N.B. that if finite differences are used to
            // linearize the system of equations, it does not matter.)
            if (upIdx[phaseIdx] == interiorIdx)
                b[phaseIdx] = up[phaseIdx]->fluidState().invB(phaseIdx);
            else
                b[phaseIdx] = Toolbox::value(up[phaseIdx]->fluidState().invB(phaseIdx));
        }

        // deal with the "main" components of each phase
        flux[conti0EqIdx + waterCompIdx] = b[waterPhaseIdx];
        flux[conti0EqIdx + waterCompIdx] *= extQuants.volumeFlux(waterPhaseIdx);
        flux[conti0EqIdx + waterCompIdx] *= FluidSystem::referenceDensity(waterPhaseIdx, up[waterPhaseIdx]->pvtRegionIndex());

        flux[conti0EqIdx + gasCompIdx] = b[gasPhaseIdx];
        flux[conti0EqIdx + gasCompIdx] *= extQuants.volumeFlux(gasPhaseIdx);
        flux[conti0EqIdx + gasCompIdx] *= FluidSystem::referenceDensity(gasPhaseIdx, up[gasPhaseIdx]->pvtRegionIndex());

        flux[conti0EqIdx + oilCompIdx] = b[oilPhaseIdx];
        flux[conti0EqIdx + oilCompIdx] *= extQuants.volumeFlux(oilPhaseIdx);
        flux[conti0EqIdx + oilCompIdx] *= FluidSystem::referenceDensity(oilPhaseIdx, up[oilPhaseIdx]->pvtRegionIndex());

        // dissolved gas (in the oil phase).
        if (FluidSystem::enableDissolvedGas()) {
            if (upIdx[oilPhaseIdx] == interiorIdx) {
                flux[conti0EqIdx + gasCompIdx] +=
                    up[oilPhaseIdx]->fluidState().Rs()
                    * up[oilPhaseIdx]->fluidState().invB(oilPhaseIdx)
                    * extQuants.volumeFlux(oilPhaseIdx)
                    * FluidSystem::referenceDensity(gasPhaseIdx, up[oilPhaseIdx]->pvtRegionIndex());
            }
            else  {
                flux[conti0EqIdx + gasCompIdx] +=
                    Toolbox::value(up[oilPhaseIdx]->fluidState().Rs())
                    * Toolbox::value(up[oilPhaseIdx]->fluidState().invB(oilPhaseIdx))
                    * extQuants.volumeFlux(oilPhaseIdx)
                    * FluidSystem::referenceDensity(gasPhaseIdx, up[oilPhaseIdx]->pvtRegionIndex());
            }
        }

        // vaporized oil (in the gas phase).
        if (FluidSystem::enableVaporizedOil()) {
            if (upIdx[gasPhaseIdx] == interiorIdx) {
                flux[conti0EqIdx + oilCompIdx] +=
                    up[gasPhaseIdx]->fluidState().Rv()
                    * up[gasPhaseIdx]->fluidState().invB(gasPhaseIdx)
                    * extQuants.volumeFlux(gasPhaseIdx)
                    * FluidSystem::referenceDensity(oilPhaseIdx, up[gasPhaseIdx]->pvtRegionIndex());
            }
            else  {
                flux[conti0EqIdx + gasCompIdx] +=
                    Toolbox::value(up[gasPhaseIdx]->fluidState().Rv())
                    * Toolbox::value(up[gasPhaseIdx]->fluidState().invB(gasPhaseIdx))
                    * extQuants.volumeFlux(gasPhaseIdx)
                    * FluidSystem::referenceDensity(oilPhaseIdx, up[gasPhaseIdx]->pvtRegionIndex());
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
