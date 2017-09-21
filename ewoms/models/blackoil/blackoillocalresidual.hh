// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::BlackOilLocalResidual
 */
#ifndef EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH
#define EWOMS_BLACK_OIL_LOCAL_RESIDUAL_HH

#include "blackoilproperties.hh"
#include "blackoilsolventmodules.hh"
#include "blackoilpolymermodules.hh"


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
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

    static const bool compositionSwitchEnabled = Indices::gasEnabled;
    static const bool waterEnabled = Indices::waterEnabled;

    static constexpr bool blackoilConserveSurfaceVolume = GET_PROP_VALUE(TypeTag, BlackoilConserveSurfaceVolume);

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef BlackOilSolventModule<TypeTag> SolventModule;
    typedef BlackOilPolymerModule<TypeTag> PolymerModule;


public:
    /*!
     * \copydoc FvBaseLocalResidual::computeStorage
     */
    template <class LhsEval>
    void computeStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                        const ElementContext& elemCtx,
                        unsigned dofIdx,
                        unsigned timeIdx) const
    {
        // retrieve the intensive quantities for the SCV at the specified point in time
        const IntensiveQuantities& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();

        storage = 0.0;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            LhsEval surfaceVolume =
                Toolbox::template decay<LhsEval>(fs.saturation(phaseIdx))
                * Toolbox::template decay<LhsEval>(fs.invB(phaseIdx))
                * Toolbox::template decay<LhsEval>(intQuants.porosity());

            storage[conti0EqIdx + activeCompIdx] += surfaceVolume;

            // account for dissolved gas
            if (phaseIdx == oilPhaseIdx && FluidSystem::enableDissolvedGas()) {
                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                storage[conti0EqIdx + activeGasCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rs())
                    * surfaceVolume;
            }

            // account for vaporized oil
            if (phaseIdx == gasPhaseIdx && FluidSystem::enableVaporizedOil()) {
                unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
                storage[conti0EqIdx + activeOilCompIdx] +=
                    Toolbox::template decay<LhsEval>(intQuants.fluidState().Rv())
                    * surfaceVolume;
            }
        }

        // convert surface volumes to component masses
        if (!blackoilConserveSurfaceVolume) {
            unsigned pvtRegionIdx = intQuants.pvtRegionIndex();
            if (waterEnabled) {
                unsigned activeWaterCompIdx = Indices::canonicalToActiveComponentIndex(waterCompIdx);
                storage[conti0EqIdx +  activeWaterCompIdx] *=
                        FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            }
            if (compositionSwitchEnabled) {
                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                storage[conti0EqIdx + activeGasCompIdx] *=
                        FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
            unsigned activeOilCompIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
            storage[conti0EqIdx + activeOilCompIdx] *=
                    FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
        }
        // deal with solvents (if present)
        SolventModule::addStorage(storage, intQuants);

        // deal with polymer (if present)
        PolymerModule::addStorage(storage, intQuants);

    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        assert(timeIdx == 0);

        flux = 0.0;

        const ExtensiveQuantities& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned focusDofIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
            if (upIdx == focusDofIdx)
                evalPhaseFluxes_<Evaluation>(flux, phaseIdx, extQuants, up);
            else
                evalPhaseFluxes_<Scalar>(flux, phaseIdx, extQuants, up);
        }

        if (!blackoilConserveSurfaceVolume) {
           for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
               unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
               unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
               const IntensiveQuantities& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
               unsigned pvtRegionIdx = up.pvtRegionIndex();
               Scalar refDensity = FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);
               flux[conti0EqIdx + activeCompIdx] *= refDensity;
           }
        }

        // deal with solvents (if present)
        SolventModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        // deal with polymer (if present)
        PolymerModule::computeFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeSource
     */
    void computeSource(RateVector& source,
                       const ElementContext& elemCtx,
                       unsigned dofIdx,
                       unsigned timeIdx) const
    {
        // retrieve the source term intrinsic to the problem
        elemCtx.problem().source(source, elemCtx, dofIdx, timeIdx);
    }

protected:
    template <class UpEval>
    void evalPhaseFluxes_(RateVector& flux,
                          unsigned phaseIdx,
                          const ExtensiveQuantities& extQuants,
                          const IntensiveQuantities& up) const
    {
        unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
        const auto& fs = up.fluidState();

        Evaluation surfaceVolumeFlux =
                Toolbox::template decay<UpEval>(fs.invB(phaseIdx))
                * extQuants.volumeFlux(phaseIdx);

        flux[conti0EqIdx + compIdx] +=
                surfaceVolumeFlux;

        // dissolved gas (in the oil phase).
        if (phaseIdx == oilPhaseIdx && FluidSystem::enableDissolvedGas()) {
            flux[conti0EqIdx + Indices::canonicalToActiveComponentIndex(gasCompIdx)] +=
                    Toolbox::template decay<UpEval>(fs.Rs())
                    * surfaceVolumeFlux;

        }

        // vaporized oil (in the gas phase).
        if (phaseIdx == gasPhaseIdx && FluidSystem::enableVaporizedOil()) {
            flux[conti0EqIdx + Indices::canonicalToActiveComponentIndex(oilCompIdx)] +=
                    Toolbox::template decay<UpEval>(fs.Rv())
                    * surfaceVolumeFlux;

        }
    }

};

} // namespace Ewoms

#endif
