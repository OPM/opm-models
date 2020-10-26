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
 * \copydoc Opm::BlackOilIndices
 */
#ifndef EWOMS_BLACK_OIL_DYN_INDICES_HH
#define EWOMS_BLACK_OIL_DYN_INDICES_HH

#include <opm/parser/eclipse/EclipseState/Runspec.hpp>

namespace Opm {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the black-oil model.
 */
template <unsigned numEqV, unsigned numWellEqV>
class BlackOilDynIndices
{
public:
    //! The number of equations
    static const int numEq = numEqV;
    static const int numWellEq = numWellEqV;
    static const int numPhases = 3;
    static const int maxPhase = 3;
    static const int maxEq = 4;

    //using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    //enum { gasCompIdx = FluidSystem::gasCompIdx };
    //enum { waterCompIdx = FluidSystem::waterCompIdx };
    //enum { oilCompIdx = FluidSystem::oilCompIdx };

    ////////
    // Component indices
    ////////

    //Inherit from BlackoilDefaultIndexTraits?
    // assumes canonical oil = 0, water = 1, gas = 2;
    static const int gasCompIdx = 2;
    static const int waterCompIdx = 1;
    static const int oilCompIdx = 0;

    ////////
    // Primary variable indices
    ////////

    //! The index of the water saturation
    static const int waterSaturationIdx = 0;

    //! Index of the oil pressure in a vector of primary variables
    static const int pressureSwitchIdx = 1;

    /*!
     * \brief Index of the switching variable which determines the composition of the
     *        hydrocarbon phases.
     *
     * Depending on the phases present, this variable is either interpreted as the
     * saturation of the gas phase, as the mole fraction of the gas component in the oil
     * phase or as the mole fraction of the oil component in the gas phase.
     */
    static const int compositionSwitchIdx = 2;

    //! Index of the primary variable for the first solvent
    static const int solventSaturationIdx =  3;

    //! Index of the primary variable for the first polymer
    static const int polymerConcentrationIdx = -1;

    //! Index of the primary variable for the second polymer primary variable (molecular weight)
    static const int polymerMoleWeightIdx = -1;

    //! Index of the primary variable for the foam
    static const int foamConcentrationIdx = -1;

    //! Index of the primary variable for the brine
    static const int saltConcentrationIdx = -1;

    //! Index of the primary variable for temperature
    static const int temperatureIdx  = -1;


    ////////
    // Equation indices
    ////////

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = 0;
    // two continuity equations follow

    //! Index of the continuity equation for the first solvent component
    static const int contiSolventEqIdx = 3;

    //! Index of the continuity equation for the first polymer component
    static const int contiPolymerEqIdx = -1;

    //! Index of the continuity equation for the second polymer component (molecular weight)
    static const int contiPolymerMWEqIdx = -1;

    //! Index of the continuity equation for the foam component
    static const int contiFoamEqIdx = -1;

    //! Index of the continuity equation for the salt water component
    static const int contiBrineEqIdx = -1;

    //! Index of the continuity equation for energy
    static const int contiEnergyEqIdx = -1;


    BlackOilDynIndices()
    {
    }

    static void init(const Phases& phases)
    {
        // The equation index map
        int numActiveComponents = 0;
        if (phases.active(Phase::OIL)){
            activeComponentMap_[oilCompIdx] = numActiveComponents;
            canonicalComponentMap_[numActiveComponents] = oilCompIdx;
            activeEquationVariableMap_[conti0EqIdx + oilCompIdx] = numActiveComponents;
            canonicalEquationVariableMap_[numActiveComponents] = oilCompIdx;
            numActiveComponents++;
        } else {
            activeComponentMap_[oilCompIdx] = -1;
            canonicalEquationVariableMap_[conti0EqIdx + oilCompIdx] = -1;
        }
        if (phases.active(Phase::WATER)){
            activeComponentMap_[waterCompIdx] = numActiveComponents;
            canonicalComponentMap_[numActiveComponents] = waterCompIdx;
            activeEquationVariableMap_[conti0EqIdx + waterCompIdx] = numActiveComponents;
            canonicalEquationVariableMap_[numActiveComponents] = waterCompIdx;
            numActiveComponents++;
        } else {
            activeComponentMap_[waterCompIdx] = -1;
            activeEquationVariableMap_[conti0EqIdx + waterCompIdx] = -1;
        }
        if (phases.active(Phase::GAS)){
            activeComponentMap_[gasCompIdx] = numActiveComponents;
            canonicalComponentMap_[numActiveComponents] = gasCompIdx;
            activeEquationVariableMap_[conti0EqIdx + gasCompIdx] = numActiveComponents;
            canonicalEquationVariableMap_[numActiveComponents] = gasCompIdx;
            numActiveComponents++;
        } else {
            activeComponentMap_[gasCompIdx] = -1;
            activeEquationVariableMap_[conti0EqIdx + gasCompIdx] = -1;
        }
        numActivePhase_ = numActiveComponents;
        for (int i = 0; i < activeComponentMap_.size(); ++i) {
            std::cout << activeComponentMap_[i] << std::endl;
        }
        for (int i = 0; i < canonicalEquationVariableMap_.size(); ++i) {
            std::cout << canonicalEquationVariableMap_[i] << std::endl;
        }

        // The equation indices
        // we always have want the pressure equation.
        int numActiveEquations = 0;

        if (phases.active(Phase::WATER) ) {
            activePrimaryVariableMap_[waterSaturationIdx] = numActiveEquations;
            canonicalPrimaryVariableMap_[numActiveEquations] = waterSaturationIdx;
            numActiveEquations++;
        } else {
            activeEquationVariableMap_[waterSaturationIdx] = -1;
        }

        activePrimaryVariableMap_[pressureSwitchIdx] = numActiveEquations;
        canonicalPrimaryVariableMap_[numActiveEquations] = pressureSwitchIdx;
        numActiveEquations++;

        if (phases.active(Phase::GAS) && phases.active(Phase::OIL)) {
            activePrimaryVariableMap_[compositionSwitchIdx] = numActiveEquations;
            canonicalPrimaryVariableMap_[numActiveEquations] = compositionSwitchIdx;
            numActiveEquations++;
        } else {
            activeEquationVariableMap_[compositionSwitchIdx] = -1;
        }

        assert(numActiveComponents == numActiveEquations);

        //The extra equtions
        if (phases.active(Phase::SOLVENT)){
            activePrimaryVariableMap_[solventSaturationIdx] = numActiveComponents;
            activeEquationVariableMap_[contiSolventEqIdx] = numActiveComponents;
            canonicalEquationVariableMap_[numActiveComponents] = contiSolventEqIdx;
            canonicalPrimaryVariableMap_[numActiveComponents] = solventSaturationIdx;
            numActiveComponents++;
        } else {
            activePrimaryVariableMap_[solventSaturationIdx] = -1;
            activeEquationVariableMap_[contiSolventEqIdx] = -1;
        }
        assert(numActiveComponents == numEq);

        for (int i = 0; i < activePrimaryVariableMap_.size(); ++i) {
            std::cout << "PV " << activePrimaryVariableMap_[i] << std::endl;
        }
        for (int i = 0; i < canonicalPrimaryVariableMap_.size(); ++i) {
            std::cout << "cpV " << canonicalPrimaryVariableMap_[i] << std::endl;
        }
#warning Add the rest

    }

    //! \brief returns the index of "active" component
    static int canonicalToActiveComponentIndex(unsigned compIdx)
    {
        assert(activeComponentMap_[compIdx] != -1);
        return activeComponentMap_[compIdx];
    }

    //! \brief returns the index of "canonical" component
    static int activeToCanonicalComponentIndex(unsigned activeCompIdx)
    {
        assert(activeCompIdx < numEq);
        return canonicalComponentMap_[activeCompIdx];
    }

    //! \brief returns the index of "active" primary variable
    static int canonicalToActivePrimaryVariableIndex(unsigned pvIdx)
    {
        assert(activePrimaryVariableMap_[pvIdx] != -1);
        return activePrimaryVariableMap_[pvIdx];
    }

    //! \brief returns the index of "canonical" primary variable
    static int activeToCanonicalPrimaryVariableIndex(unsigned activePvIdx)
    {
        assert(activePvIdx < numEq);
        return canonicalPrimaryVariableMap_[activePvIdx];
    }

    static int activePressureSwitchIdx()
    {
        return canonicalToActivePrimaryVariableIndex(pressureSwitchIdx);
    }

    static int activeWaterSaturationIdx()
    {
        return canonicalToActivePrimaryVariableIndex(waterSaturationIdx);
    }

    static int activeCompositionSwitchIdx()
    {
        return canonicalToActivePrimaryVariableIndex(compositionSwitchIdx);
    }

    static int activeSolventSaturationIdx()
    {
        return canonicalToActivePrimaryVariableIndex(solventSaturationIdx);
    }

    //! \brief returns the index of "active" component
    static int canonicalToActiveEquationVariableIndex(unsigned eqIdx)
    {
        assert(activeEquationVariableMap_[eqIdx] != -1);
        return activeEquationVariableMap_[eqIdx];
    }

    //! \brief returns the index of "canonical" component
    static int activeToCanonicalEquationVariableIndex(unsigned activeEqIdx)
    {
        assert(activeEqIdx < numEq);
        return canonicalEquationVariableMap_[activeEqIdx];
    }

    static int activeContiSolventEqIdx()
    {
        return canonicalToActiveEquationVariableIndex(contiSolventEqIdx);
    }

    static bool waterIsActive()
    {
        return activeComponentMap_[waterCompIdx] != -1;
    }

    static bool oilIsActive()
    {
        return activeComponentMap_[oilCompIdx] != -1;
    }

    static bool gasIsActive()
    {
        return activeComponentMap_[gasCompIdx] != -1;
    }

    static bool solventIsActive()
    {
        return activePrimaryVariableMap_[solventSaturationIdx] != -1;
    }

    static int numActivePhase()
    {
        return numActivePhase_;
    }

private:
    static int numActivePhase_;
    static std::array<int, maxPhase> activeComponentMap_;
    static std::array<int, numPhases> canonicalComponentMap_;
    static std::array<int, maxEq> activePrimaryVariableMap_;
    static std::array<int, numEq> canonicalPrimaryVariableMap_;
    static std::array<int, maxEq> activeEquationVariableMap_;
    static std::array<int, numEq> canonicalEquationVariableMap_;
};
    template <unsigned numEq, unsigned numWellEq>
    int BlackOilDynIndices<numEq, numWellEq>::numActivePhase_;

    template <unsigned numEq, unsigned numWellEq>
    std::array<int, BlackOilDynIndices<numEq, numWellEq>::maxPhase> BlackOilDynIndices<numEq, numWellEq>::activeComponentMap_;

    template <unsigned numEq, unsigned numWellEq>
    std::array<int, BlackOilDynIndices<numEq, numWellEq>::numPhases> BlackOilDynIndices<numEq, numWellEq>::canonicalComponentMap_;

    template <unsigned numEq, unsigned numWellEq>
    std::array<int, BlackOilDynIndices<numEq, numWellEq>::maxEq> BlackOilDynIndices<numEq, numWellEq>::activePrimaryVariableMap_;

    template <unsigned numEq, unsigned numWellEq>
    std::array<int, BlackOilDynIndices<numEq, numWellEq>::numEq> BlackOilDynIndices<numEq, numWellEq>::canonicalPrimaryVariableMap_;

    template <unsigned numEq, unsigned numWellEq>
    std::array<int, BlackOilDynIndices<numEq, numWellEq>::maxEq> BlackOilDynIndices<numEq, numWellEq>::activeEquationVariableMap_;

    template <unsigned numEq, unsigned numWellEq>
    std::array<int, BlackOilDynIndices<numEq, numWellEq>::numEq> BlackOilDynIndices<numEq, numWellEq>::canonicalEquationVariableMap_;


} // namespace Opm

#endif
