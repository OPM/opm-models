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
 * \copydoc Ewoms::BlackOilVolumeVariables
 */
#ifndef EWOMS_BLACK_OIL_VOLUME_VARIABLES_HH
#define EWOMS_BLACK_OIL_VOLUME_VARIABLES_HH

#include "blackoilproperties.hh"

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <dune/common/fmatrix.hh>

namespace Ewoms {
/*!
 * \ingroup BlackOilModel
 * \ingroup VolumeVariables
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model.
 */
template <class TypeTag>
class BlackOilVolumeVariables
    : public GET_PROP_TYPE(TypeTag, DiscVolumeVariables)
    , public GET_PROP_TYPE(TypeTag, VelocityModule)::VelocityVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscVolumeVariables) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, BlackOilFluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, VelocityModule) VelocityModule;

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

    typedef typename VelocityModule::VelocityVolumeVariables VelocityVolumeVariables;

public:
    /*!
     * \copydoc VolumeVariables::update
     */
    void update(const ElementContext &elemCtx,
                int dofIdx,
                int timeIdx)
    {
        ParentType::update(elemCtx,
                           dofIdx,
                           timeIdx);

        fluidState_.setTemperature(elemCtx.problem().temperature(elemCtx, dofIdx, timeIdx));

        const auto &problem = elemCtx.problem();
        const auto &priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // extract the water and the gas saturations for convenience
        Scalar Sw = priVars[Indices::waterSaturationIdx];

        Scalar Sg = 0.0;
        if (priVars.switchingVariableIsGasSaturation())
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

        // oil phase pressure
        Scalar po = fluidState_.pressure(oilPhaseIdx);

        // update phase compositions. first, set everything to 0, then
        // make the gas/water phases consist of only the gas/water
        // components
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState_.setMoleFraction(phaseIdx, compIdx, 0.0);
        // set composition of gas and water phases
        fluidState_.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);
        fluidState_.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);

        // retrieve the relevant parameters from the fluid
        // system.
        Scalar Bo = FluidSystem::oilFormationVolumeFactor(po);
        Scalar Bg = FluidSystem::gasFormationVolumeFactor(pg);
        Scalar Rs = FluidSystem::gasDissolutionFactor(pg);

        Scalar rhorefo = FluidSystem::surfaceDensity(oilPhaseIdx);
        Scalar rhorefg = FluidSystem::surfaceDensity(gasPhaseIdx);

        Scalar rhoo = rhorefo / Bo;
        Scalar rhog = rhorefg / Bg;

        Scalar MG = FluidSystem::molarMass(gasPhaseIdx);
        Scalar MO = FluidSystem::molarMass(oilPhaseIdx);

        if (priVars.switchingVariableIsGasSaturation()) {
            // we need to calculate the composition of the
            // gas-saturated oil phase
            //
            // for this, we first calculate composition of the
            // saturated oil phase in terms of mass fractions.
            Scalar XoG = Rs * rhorefg / rhoo;

            // Then we convert these to mole fractions. This assumes
            // that the gas and water phases are immiscible and only
            // consist of their respective components.
            Scalar avgMolarMass = MO / (1 + XoG * (MO / MG - 1));
            Scalar xoG = XoG * avgMolarMass / MG;
            Scalar xoO = 1 - xoG;

            // now we can finally set the oil-phase composition
            fluidState_.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
            fluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, xoO);
        }
        else {
            // if the switching variable is not the gas saturation, it
            // is the mole fraction of the gas component in the oil
            // phase.
            Scalar xoG = priVars[Indices::switchIdx];

            // this makes setting the oil-phase composition simple
            fluidState_.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
            fluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, 1 - xoG);
        }

        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        // set the phase densities
        fluidState_.setDensity(waterPhaseIdx,
                               FluidSystem::density(fluidState_, paramCache, waterPhaseIdx));
        fluidState_.setDensity(gasPhaseIdx, rhog);
        fluidState_.setDensity(oilPhaseIdx,
                               FluidSystem::density(fluidState_, paramCache, oilPhaseIdx));

        // compute the phase viscosities
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            fluidState_.setViscosity(phaseIdx, mu);
        }

        // calculate relative permeabilities
        MaterialLaw::relativePermeabilities(relativePermeability_, materialParams, fluidState_);
        Valgrind::CheckDefined(relativePermeability_);

        // retrieve the porosity from the problem ...
        porosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);

        // ... as well as the intrinsic permeability
        intrinsicPerm_ = problem.intrinsicPermeability(elemCtx, dofIdx, timeIdx);

        // update the quantities which are required by the chosen
        // velocity model
        VelocityVolumeVariables::update_(elemCtx, dofIdx, timeIdx);

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
     * \copydoc ImmiscibleVolumeVariables::fluidState
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::intrinsicPermeability
     */
    const DimMatrix &intrinsicPermeability() const
    { return intrinsicPerm_; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::relativePermeability
     */
    Scalar relativePermeability(int phaseIdx) const
    { return relativePermeability_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleVolumeVariables::mobility
     */
    Scalar mobility(int phaseIdx) const
    {
        return relativePermeability(phaseIdx) / fluidState().viscosity(phaseIdx);
    }

    /*!
     * \copydoc ImmiscibleVolumeVariables::porosity
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
