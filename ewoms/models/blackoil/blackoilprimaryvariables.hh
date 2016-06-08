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
 * \copydoc Ewoms::BlackOilPrimaryVariables
 */
#ifndef EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH
#define EWOMS_BLACK_OIL_PRIMARY_VARIABLES_HH

#include "blackoilproperties.hh"

#include <ewoms/disc/common/fvbaseprimaryvariables.hh>

#include <dune/common/fvector.hh>

#include <opm/material/constraintsolvers/NcpFlash.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

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
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // number of equations
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    // primary variable indices
    enum { waterSaturationIdx = Indices::waterSaturationIdx };
    enum { pressureSwitchIdx = Indices::pressureSwitchIdx };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

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
    enum PrimaryVarsMeaning {
        Sw_po_Sg, // threephase case
        Sw_po_Rs, // water + oil case
        Sw_po_Rv, // water + gas case
    };

    BlackOilPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetUndefined(*this);
        pvtRegionIdx_ = 0;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    BlackOilPrimaryVariables(Scalar value)
        : ParentType(value)
    {
        Valgrind::SetUndefined(primaryVarsMeaning_);
        pvtRegionIdx_ = 0;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables &)
     */
    BlackOilPrimaryVariables(const BlackOilPrimaryVariables &value) = default;

    void setPvtRegionIndex(int value)
    { pvtRegionIdx_ = value; }

    unsigned pvtRegionIndex() const
    { return pvtRegionIdx_; }

    void setOilPressure( const Scalar& po )
    {
        (*this)[ pressureSwitchIdx ] = po;
    }

    void setWaterSaturation( const Scalar& sw )
    {
        (*this)[ waterSaturationIdx ] = sw;
    }

    void setSwitchingVariable( const Scalar& x )
    {
        (*this)[ compositionSwitchIdx ] = x;
    }

    PrimaryVarsMeaning primaryVarsMeaning() const
    { return primaryVarsMeaning_; }

    void setPrimaryVarsMeaning(PrimaryVarsMeaning newMeaning)
    { primaryVarsMeaning_ = newMeaning; }

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
        typename FluidSystem::template ParameterCache<Scalar> paramCache;
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

        paramCache.updateAll(fsFlash);
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
        NcpFlash::template solve<MaterialLaw>(fsFlash, matParams, paramCache, globalMolarities);

        // use the result to assign the primary variables
        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        typedef typename std::remove_reference<typename FluidState::Scalar>::type ConstEvaluation;
        typedef typename std::remove_const<ConstEvaluation>::type FsEvaluation;
        typedef typename Opm::MathToolbox<FsEvaluation> FsToolbox;

        bool gasPresent = (fluidState.saturation(gasPhaseIdx) > 0.0);
        bool oilPresent = (fluidState.saturation(oilPhaseIdx) > 0.0);

        // determine the meaning of the primary variables
        if ((gasPresent && oilPresent) || (!gasPresent && !oilPresent))
            // gas and oil: both hydrocarbon phases are in equilibrium (i.e., saturated
            // with the "protagonist" component of the other phase.)
            primaryVarsMeaning_ = Sw_po_Sg;
        else if (oilPresent) {
            // only oil: if dissolved gas is enabled, we need to consider the oil phase
            // composition, if it is disabled, the gas component must stick to its phase
            if (FluidSystem::enableDissolvedGas())
                primaryVarsMeaning_ = Sw_po_Rs;
            else
                primaryVarsMeaning_ = Sw_po_Sg;
        }
        else {
            assert(gasPresent);
            // only gas: if vaporized oil is enabled, we need to consider the gas phase
            // composition, if it is disabled, the oil component must stick to its phase
            if (FluidSystem::enableVaporizedOil())
                primaryVarsMeaning_ = Sw_po_Rv;
            else
                primaryVarsMeaning_ = Sw_po_Sg;
        }

        // assign the actual primary variables
        if (primaryVarsMeaning() == Sw_po_Sg) {
            (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(oilPhaseIdx));
            (*this)[compositionSwitchIdx] = FsToolbox::value(fluidState.saturation(gasPhaseIdx));
        }
        else if (primaryVarsMeaning() == Sw_po_Rs) {
            const auto& Rs = Opm::BlackOil::getRs_<FluidSystem, Scalar, FluidState>(fluidState, pvtRegionIdx_);

            (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(oilPhaseIdx));
            (*this)[compositionSwitchIdx] = Rs;
        }
        else {
            assert(primaryVarsMeaning() == Sw_po_Rv);

            const auto& Rv = Opm::BlackOil::getRv_<FluidSystem, Scalar, FluidState>(fluidState, pvtRegionIdx_);
            (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(gasPhaseIdx));
            (*this)[compositionSwitchIdx] = Rv;
        }
    }

    /*!
     * \brief Adapt the interpretation of the switching variables to be physically
     *        meaningful.
     *
     * \return true Iff the interpretation of one of the switching variables was changed
     */
    bool adaptPrimaryVariables()
    {
        // this function accesses some low level functions directly for better
        // performance (instead of going the canonical way through the
        // IntensiveQuantities). The reason is that most intensive quantities are not
        // required to be able to decide if the primary variables needs to be switched or
        // not, so it would be a waste to compute them.

        Scalar Sw = (*this)[Indices::waterSaturationIdx];
        if (primaryVarsMeaning() == Sw_po_Sg) {
            // both hydrocarbon phases are present.
            Scalar Sg = (*this)[Indices::compositionSwitchIdx];
            Scalar So = 1.0 - Sw - Sg;

            Scalar So2 = 1.0 - Sw;
            if (Sg < 0.0 && So2 > 0.0 && FluidSystem::enableDissolvedGas()) {
                // the gas phase disappeared, i.e., switch the primary variables to { Sw,
                // po, xoG }.
                //
                // by a lucky coincidence the pressure switching variable already
                // represents the oil phase pressure, so we do not need to change
                // this. For the gas mole fraction, we use the low level blackoil PVT
                // objects to calculate the mole fraction of gas saturated oil.
                Scalar po = (*this)[Indices::pressureSwitchIdx];
                Scalar T = asImp_().temperature_();
                Scalar RsSat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, T, po);

                setPrimaryVarsMeaning(Sw_po_Rs);
                (*this)[Indices::pressureSwitchIdx] = po;
                (*this)[Indices::compositionSwitchIdx] = RsSat;

                return true;
            }

            Scalar Sg2 = 1.0 - Sw;
            if (So < 0.0 && Sg2 > 0.0 && FluidSystem::enableVaporizedOil()) {
                // the oil phase disappeared, i.e., switch the primary variables to { Sw,
                // pg, xgO }.
                //
                // here we're less lucky than above because we only have the oil phas
                // pressure but need the gas phase pressure. Since the required
                // quantities for calculating the capillary pressures are not easily
                // available here, we just assume that the capillary pressures are 0,
                // i.e., that gas phase pressure equals the oil phase pressure.
                Scalar pg = (*this)[Indices::pressureSwitchIdx]; // TODO: capillary pressure
                Scalar T = asImp_().temperature_();
                Scalar RvSat = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, T, pg);

                setPrimaryVarsMeaning(Sw_po_Rv);
                (*this)[Indices::pressureSwitchIdx] = pg;
                (*this)[Indices::compositionSwitchIdx] = RvSat;

                return true;
            }

            return false;
        }
        else if (primaryVarsMeaning() == Sw_po_Rs) {
            // only the oil and the water phases are present. the gas phase appears as
            // soon as more of the gas component is present in the oil phase than what
            // the saturated phase contains. Note that we use the blackoil specific
            // low-level PVT objects here for performance reasons.
            Scalar T = asImp_().temperature_();
            Scalar po = (*this)[Indices::pressureSwitchIdx];
            Scalar So = 1 - Sw;

            if (So <= 0.0) {
                // switch back to phase equilibrium mode if the oil phase vanishes (i.e.,
                // the water-only case)
                setPrimaryVarsMeaning(Sw_po_Sg);
                (*this)[Indices::waterSaturationIdx] = 1.0; // water saturation
                (*this)[Indices::pressureSwitchIdx] = po;
                (*this)[Indices::compositionSwitchIdx] = 0.0; // gas saturation

                return true;
            }

            Scalar RsSat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_, T, po);

            Scalar Rs = (*this)[Indices::compositionSwitchIdx];
            if (Rs > RsSat) {
                // the gas phase appears, i.e., switch the primary variables to { Sw, po,
                // Sg }.
                setPrimaryVarsMeaning(Sw_po_Sg);
                (*this)[Indices::pressureSwitchIdx] = po; // TODO (?): capillary pressure
                (*this)[Indices::compositionSwitchIdx] = 0.0; // gas saturation

                return true;
            }

            return false;
        }
        else {
            assert(primaryVarsMeaning() == Sw_po_Rv);

            // only the gas and the water phases are present. the oil phase appears as
            // soon as more of the oil component is present in the gas phase than what
            // the saturated phase contains. Note that we use the blackoil specific
            // low-level PVT objects here for performance reasons.
            Scalar T = asImp_().temperature_();
            Scalar pg = (*this)[Indices::pressureSwitchIdx];
            Scalar Sg = 1 - Sw;

            if (Sg <= 0.0) {
                // switch back to phase equilibrium mode if the gas phase also vanishes
                Scalar po = pg;

                setPrimaryVarsMeaning(Sw_po_Sg);
                (*this)[Indices::waterSaturationIdx] = 1.0; // water saturation
                (*this)[Indices::pressureSwitchIdx] = po;
                (*this)[Indices::compositionSwitchIdx] = 0.0; // gas saturation

                return true;
            }

            Scalar RvSat = FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_, T, pg);

            Scalar Rv = (*this)[Indices::compositionSwitchIdx];
            if (Rv > RvSat) {
                // the oil phase appears, i.e., switch the primary variables to { Sw,
                // po, Sg }.
                //
                // For this, we should calculate the oil phase pressure, but since this
                // it not easily possible at this point, let's just assume a capillary
                // pressure of 0 here. (like above, but for the switch into the opposite
                // direction.)
                Scalar po = pg;

                setPrimaryVarsMeaning(Sw_po_Sg);
                (*this)[Indices::pressureSwitchIdx] = po;
                (*this)[Indices::compositionSwitchIdx] = 1.0 - Sw; // gas saturation, i.e., So = 0

                return true;
            }

            return false;
        }

        assert(false);
        return false;
    }

    BlackOilPrimaryVariables& operator=(const BlackOilPrimaryVariables& other) = default;
    BlackOilPrimaryVariables& operator=(Scalar value)
    {
        for (int i = 0; i < numEq; ++i)
            (*this)[i] = value;

        return *this;
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the standard blackoil model is isothermal
    Scalar temperature_() const
    { return FluidSystem::surfaceTemperature; }

    PrimaryVarsMeaning primaryVarsMeaning_;
    unsigned char pvtRegionIdx_;
};


} // namespace Ewoms

#endif
