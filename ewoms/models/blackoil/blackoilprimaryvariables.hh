// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2014 by Andreas Lauser

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
 * \copydoc Ewoms::BlackOilPrimaryVariables
 */
#ifndef EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH
#define EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH

#include "blackoilproperties.hh"

#include <ewoms/disc/common/fvbaseprimaryvariables.hh>

#include <dune/common/fvector.hh>

#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

namespace Ewoms {
/*!
 * \ingroup BlackOilModel
 *
 * \brief Represents the primary variables used by the black-oil model.
 */
template <class TypeTag>
class BlackOilPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    typedef FvBasePrimaryVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // number of equations
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    // primary variable indices
    enum { gasPressureIdx = Indices::gasPressureIdx };
    enum { waterSaturationIdx = Indices::waterSaturationIdx };
    enum { switchIdx = Indices::switchIdx };

    // phase indices from the fluid system
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    // component indices from the fluid system
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

    typedef typename Opm::MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    static_assert(numPhases == 3, "The black-oil model assumes three phases!");
    static_assert(numComponents == 3, "The black-oil model assumes three components!");

public:
    enum SwitchingVarMeaning {
        GasSaturation,
        OilMoleFractionInGas,
        GasMoleFractionInOil
    };

    BlackOilPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetUndefined(*this);
        pvtRegionIdx_ = 0;
        temperature_ = FluidSystem::surfaceTemperature;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    BlackOilPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::SetUndefined(switchingVarMeaning_);
        pvtRegionIdx_ = 0;
        temperature_ = FluidSystem::surfaceTemperature;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables &)
     */
    BlackOilPrimaryVariables(const BlackOilPrimaryVariables &value)
        : ParentType(value)
        , temperature_(FluidSystem::surfaceTemperature)
        , switchingVarMeaning_(value.switchingVarMeaning_)
        , pvtRegionIdx_(value.pvtRegionIdx_)
    { }

    void setPvtRegionIndex(int value)
    { pvtRegionIdx_ = value; }

    unsigned pvtRegionIndex() const
    { return pvtRegionIdx_; }

    SwitchingVarMeaning switchingVarMeaning() const
    { return switchingVarMeaning_; }

    void setSwitchingVarMeaning(SwitchingVarMeaning newMeaning)
    { switchingVarMeaning_ = newMeaning; }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState &fluidState,
                                const MaterialLawParams &matParams,
                                bool isInEquilibrium = false)
    {
        typedef typename std::remove_reference<typename FluidState::Scalar>::type ConstEvaluation;
        typedef typename std::remove_const<ConstEvaluation>::type FsEvaluation;
        typedef typename Opm::MathToolbox<FsEvaluation> FsToolbox;

#ifndef NDEBUG
        // make sure the temperature is the same in all fluid phases
        for (int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx) {
            Valgrind::CheckDefined(fluidState.temperature(0));
            Valgrind::CheckDefined(fluidState.temperature(phaseIdx));

            assert(fluidState.temperature(0) == fluidState.temperature(phaseIdx));
        }
#endif // NDEBUG

        // for the equilibrium case, we don't need complicated
        // computations.
        if (isInEquilibrium) {
            assignNaive(fluidState);
            return;
        }

        // If your compiler bails out here, you're probably not using a suitable black
        // oil fluid system.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.setRegionIndex(pvtRegionIdx_);

        // create a mutable fluid state with well defined densities based on the input
        typedef Opm::NcpFlash<Scalar, FluidSystem> NcpFlash;
        typedef Opm::CompositionalFluidState<Scalar, FluidSystem> FlashFluidState;
        FlashFluidState fsFlash;
        fsFlash.setTemperature(FsToolbox::value(fluidState.temperature(/*phaseIdx=*/0)));
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fsFlash.setPressure(phaseIdx, FsToolbox::value(fluidState.pressure(phaseIdx)));
            fsFlash.setSaturation(phaseIdx, FsToolbox::value(fluidState.saturation(phaseIdx)));
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fsFlash.setMoleFraction(phaseIdx, compIdx, FsToolbox::value(fluidState.moleFraction(phaseIdx, compIdx)));
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::template density<FlashFluidState, Scalar>(fsFlash, paramCache, phaseIdx);
            fsFlash.setDensity(phaseIdx, rho);
        }

        // calculate the "global molarities"
        ComponentVector globalMolarities(0.0);
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx] +=
                    fsFlash.saturation(phaseIdx) * fsFlash.molarity(phaseIdx, compIdx);
            }
        }

        // use a flash calculation to calculate a fluid state in
        // thermodynamic equilibrium

        // run the flash calculation
        NcpFlash::template solve<MaterialLaw>(fsFlash, paramCache, matParams, globalMolarities);

        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        typedef typename std::remove_reference<typename FluidState::Scalar>::type ConstEvaluation;
        typedef typename std::remove_const<ConstEvaluation>::type FsEvaluation;
        typedef typename Opm::MathToolbox<FsEvaluation> FsToolbox;

        temperature_ = FsToolbox::value(fluidState.temperature(gasPhaseIdx));

        // the pressure of the first phase and the saturation of water
        // are always primary variables
        (*this)[gasPressureIdx] = FsToolbox::value(fluidState.pressure(gasPhaseIdx));
        (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));

        bool gasPresent = (fluidState.saturation(gasPhaseIdx) > 0.0);
        bool oilPresent = (fluidState.saturation(oilPhaseIdx) > 0.0);

        if ((gasPresent && oilPresent) || (!gasPresent && !oilPresent))
            // gas and oil
            switchingVarMeaning_ = GasSaturation;
        else if (oilPresent) {
            // only oil
            if (FluidSystem::enableDissolvedGas())
                switchingVarMeaning_ = GasMoleFractionInOil;
            else
                switchingVarMeaning_ = GasSaturation;
        }
        else {
            assert(gasPresent);
            // only gas
            if (FluidSystem::enableVaporizedOil())
                switchingVarMeaning_ = OilMoleFractionInGas;
            else
                switchingVarMeaning_ = GasSaturation;
        }

        // depending on the phases present, we use either Sg, x_o^G or x_g^O as the
        // switching primary variable
        if (switchingVarMeaning() == GasSaturation)
            (*this)[switchIdx] = FsToolbox::value(fluidState.saturation(gasPhaseIdx));
        else if (switchingVarMeaning() == GasMoleFractionInOil)
            (*this)[switchIdx] = FsToolbox::value(fluidState.moleFraction(oilPhaseIdx, gasCompIdx));
        else {
            assert(switchingVarMeaning() == OilMoleFractionInGas);
            (*this)[switchIdx] = FsToolbox::value(fluidState.moleFraction(gasPhaseIdx, oilCompIdx));
        }
    }

    /*!
     * \brief Adapt the interpretation of the switching variable to a
     *        physically meaningful one.
     *
     * \return true Iff the interpretation of the switching variable was changed
     */
    bool adaptSwitchingVariable()
    {
        Scalar pg = (*this)[Indices::gasPressureIdx];

        if (switchingVarMeaning() == GasSaturation) {
            // both hydrocarbon phases
            Scalar Sw = (*this)[Indices::waterSaturationIdx];
            Scalar Sg = (*this)[Indices::switchIdx];
            Scalar So = 1 - Sw - Sg;

            if (Sg < 0.0 && So > 0.0 && FluidSystem::enableDissolvedGas()) {
                // we switch to the gas mole fraction in the oil phase if some oil phase
                // is present and dissolved gas is enabled
                Scalar Rssat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temperature_, pg);
                Scalar XoGsat = FluidSystem::convertRsToXoG(Rssat, pvtRegionIdx_);
                Scalar xoGsat = FluidSystem::convertXoGToxoG(XoGsat, pvtRegionIdx_);

                setSwitchingVarMeaning(GasMoleFractionInOil);
                (*this)[Indices::switchIdx] = xoGsat;
                return true;
            }

            if (So < 0.0 && Sg > 0.0 && FluidSystem::enableVaporizedOil()) {
                // we switch to the oil mole fraction in the gas phase if some gas phase
                // is present and vaporized oil is enabled. TODO (?): use oil instead of
                // gas pressure!
                Scalar Rvsat = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temperature_, pg);
                Scalar XgOsat = FluidSystem::convertRvToXgO(Rvsat, pvtRegionIdx_);
                Scalar xgOsat = FluidSystem::convertXgOToxgO(XgOsat, pvtRegionIdx_);

                setSwitchingVarMeaning(OilMoleFractionInGas);
                (*this)[Indices::switchIdx] = xgOsat;
                return true;
            }
            return false;
        }
        else if (switchingVarMeaning() == OilMoleFractionInGas) {
            // only gas. oil appears as soon as more oil is present as can be dissolved
            Scalar Sw = (*this)[Indices::waterSaturationIdx];
            Scalar xgO = (*this)[Indices::switchIdx];
            Scalar Rvsat = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, temperature_, pg);
            Scalar XgOsat = FluidSystem::convertRvToXgO(Rvsat, pvtRegionIdx_);
            Scalar xgOsat = FluidSystem::convertXgOToxgO(XgOsat, pvtRegionIdx_);

            if (xgO > xgOsat) {
                setSwitchingVarMeaning(GasSaturation);
                (*this)[Indices::switchIdx] = 1 - Sw;
                return true;
            }

            return false;
        }
        else {
            assert(switchingVarMeaning() == GasMoleFractionInOil);

            // only oil. gas appears as soon as more gas is present as can be dissolved
            Scalar xoG = (*this)[Indices::switchIdx];
            Scalar Rssat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, temperature_, pg);
            Scalar XoGsat = FluidSystem::convertRsToXoG(Rssat, pvtRegionIdx_);
            Scalar xoGsat = FluidSystem::convertXoGToxoG(XoGsat, pvtRegionIdx_);

            if (xoG > xoGsat) {
                setSwitchingVarMeaning(GasSaturation);
                (*this)[Indices::switchIdx] = 0.0;
                return true;
            }

            return false;
        }

        assert(false);
        return false;
    }

    BlackOilPrimaryVariables& operator=(const BlackOilPrimaryVariables& other)
    {
        ParentType::operator=(other);
        temperature_ = other.temperature_;
        switchingVarMeaning_ = other.switchingVarMeaning_;
        pvtRegionIdx_ = other.pvtRegionIdx_;
        return *this;
    }

    BlackOilPrimaryVariables& operator=(Scalar value)
    {
        for (int i = 0; i < numEq; ++i)
            (*this)[i] = value;

        return *this;
    }

private:
    Scalar temperature_; // so far this is just a pseudo-primary variable
    SwitchingVarMeaning switchingVarMeaning_;
    unsigned char pvtRegionIdx_;
};


} // namespace Ewoms

#endif
