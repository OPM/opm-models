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
 * \copydoc Ewoms::BlackOilIntensiveQuantities
 */
#ifndef EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_HH
#define EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_HH

#include "blackoilproperties.hh"

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <dune/common/fmatrix.hh>

namespace Ewoms {
/*!
 * \ingroup BlackOilModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model.
 */
template <class TypeTag>
class BlackOilIntensiveQuantities
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BlackOilFluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;

public:
    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext &elemCtx,
                int dofIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        fluidState_.setTemperature(elemCtx.problem().temperature(elemCtx, dofIdx, timeIdx));

        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        int pvtRegionIdx = priVars.pvtRegionIndex();

        // extract the water and the gas saturations for convenience
        Scalar Sw = priVars[Indices::waterSaturationIdx];

        Scalar Sg = 0.0;
        if (priVars.switchingVarMeaning() == PrimaryVariables::GasSaturation)
            Sg = priVars[Indices::switchIdx];

        fluidState_.setSaturation(waterPhaseIdx, Sw);
        fluidState_.setSaturation(gasPhaseIdx, Sg);
        fluidState_.setSaturation(oilPhaseIdx, 1 - Sw - Sg);

        // reference phase (-> gas) pressure
        Scalar pg = priVars[Indices::gasPressureIdx];

        // now we compute all phase pressures
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        Scalar pC[numPhases];
        const auto &materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
            if (fluidState_.pressure(phaseIdx) < 1e5) {
                OPM_THROW(Opm::NumericalProblem,
                          "All pressures must be at least 1 bar.");
            }
        }

        // oil phase temperature and pressure
        Scalar T = fluidState_.temperature(oilPhaseIdx);

        // update phase compositions. first, set everything to 0...
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState_.setMoleFraction(phaseIdx, compIdx, 0.0);

        // ... then set the default composition of all phases (default = assume immscibility)
        fluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, 1.0);
        fluidState_.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);
        fluidState_.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);

        // take the meaning of the switiching primary variable into account for the gas
        // and oil phase compositions
        if (priVars.switchingVarMeaning() == PrimaryVariables::GasSaturation) {
            // gas is present, i.e. we use the compositions of the gas-saturated oil and
            // oil-saturated gas phases if dissolved gas and vaporized oil are
            // enabled. if they are disabled, continue to assume immiscibility for the
            // respective phase.
            Scalar xgO = 0.0;
            Scalar xoG = 0.0;

            if (FluidSystem::enableVaporizedOil())
                xgO = FluidSystem::saturatedGasOilMoleFraction(T, pg, /*regionIdx=*/0);

            else if (FluidSystem::enableDissolvedGas())
                // TODO (?): use oil instead of gas pressure
                xoG = FluidSystem::saturatedOilGasMoleFraction(T, pg, /*regionIdx=*/0);

            fluidState_.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
            fluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, 1 - xoG);

            fluidState_.setMoleFraction(gasPhaseIdx, gasCompIdx, 1 - xgO);
            fluidState_.setMoleFraction(gasPhaseIdx, oilCompIdx, xgO);
        }
        else if (priVars.switchingVarMeaning() == PrimaryVariables::GasMoleFractionInOil) {
            // if the switching variable is the mole fraction of the gas component in the
            // oil phase, we can directly set the composition of the oil phase
            Scalar xoG = priVars[Indices::switchIdx];

            fluidState_.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
            fluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, 1 - xoG);
        }
        else {
            assert(priVars.switchingVarMeaning() == PrimaryVariables::GasMoleFractionInOil);

            // if the switching variable is the mole fraction of the oil component in the
            // gas phase, we can directly set the composition of the gas phase
            Scalar xgO = priVars[Indices::switchIdx];

            fluidState_.setMoleFraction(gasPhaseIdx, gasCompIdx, 1 - xgO);
            fluidState_.setMoleFraction(gasPhaseIdx, oilCompIdx, xgO);
        }

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::ParameterCache paramCache;
        paramCache.setRegionIndex(pvtRegionIdx);
        paramCache.updateAll(fluidState_);

        // set the phase densities and viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fluidState_, paramCache, phaseIdx);
            fluidState_.setDensity(phaseIdx, rho);

            Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);
        }

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // retrieve the porosity from the problem
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);

        // the porosity must be modified by the compressibility of the
        // rock...
        Scalar rockCompressibility = problem.rockCompressibility(elemCtx, dofIdx, timeIdx);
        if (rockCompressibility > 0.0) {
            Scalar rockRefPressure = problem.rockReferencePressure(elemCtx, dofIdx, timeIdx);
            Scalar x = rockCompressibility*(pg - rockRefPressure);
            porosity_ *= 1.0 + x + 0.5*x*x;
        }

        // now get the intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities which are required by the chosen
        // velocity model
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

#ifndef NDEBUG
        // some safety checks in debug mode
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            assert(std::isfinite(fluidState_.density(phaseIdx)));
            assert(std::isfinite(fluidState_.saturation(phaseIdx)));
            assert(std::isfinite(fluidState_.temperature(phaseIdx)));
            assert(std::isfinite(fluidState_.pressure(phaseIdx)));
            assert(std::isfinite(fluidState_.viscosity(phaseIdx)));
            assert(std::isfinite(relativePermeability_[phaseIdx]));
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
                assert(std::isfinite(fluidState_.moleFraction(phaseIdx, compIdx)));
            }
        }
        assert(std::isfinite(intrinsicPerm_.frobenius_norm()));
        assert(std::isfinite(porosity_));
#endif
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::intrinsicPermeability
     */
    const DimMatrix &intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    Scalar mobility(int phaseIdx) const
    {
        return relativePermeability(phaseIdx) / fluidState().viscosity(phaseIdx);
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    Scalar porosity() const
    { return porosity_; }

private:
    FluidState fluidState_;
    Scalar porosity_;
    DimMatrix intrinsicPerm_;
    Scalar relativePermeability_[numPhases];
};

} // namespace Ewoms

#endif
