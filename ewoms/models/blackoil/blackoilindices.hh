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
 * \copydoc Ewoms::BlackOilIndices
 */
#ifndef EWOMS_BLACK_OIL_INDICES_HH
#define EWOMS_BLACK_OIL_INDICES_HH

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the black-oil model.
 */
template <bool enableSolventV, bool enablePolymerV, bool enableEnergyV, unsigned PVOffset>
struct BlackOilIndices
{
    //! Number of phases active at all times
    static const int numPhases = 3;

    //! All phases are enabled
    static const bool oilEnabled = true;
    static const bool waterEnabled = true;
    static const bool gasEnabled = true;

    //! Are solvents involved?
    static const bool enableSolvent = enableSolventV;

    //! Are polymers involved?
    static const bool enablePolymer = enablePolymerV;

    //! Shall energy be conserved?
    static const bool enableEnergy = enableEnergyV;

private:
    //! Number of solvent components to be considered
    static const int numSolvents_ = enableSolvent ? 1 : 0;

    //! Number of polymer components to be considered
    static const int numPolymers_ = enablePolymer ? 1 : 0;

    //! Number of energy equations to be considered
    static const int numEnergy_ = enableEnergy ? 1 : 0;

public:
    //! The number of equations
    static const int numEq = numPhases + numSolvents_ + numPolymers_ + numEnergy_;

    //! \brief returns the index of "active" component
    static constexpr unsigned canonicalToActiveComponentIndex(unsigned compIdx)
    { return compIdx; }

    static constexpr unsigned activeToCanonicalComponentIndex(unsigned compIdx)
    { return compIdx; }

    ////////
    // Primary variable indices
    ////////

    //! The index of the water saturation
    static const int waterSaturationIdx = PVOffset + 0;

    //! Index of the oil pressure in a vector of primary variables
    static const int pressureSwitchIdx = PVOffset + 1;

    /*!
     * \brief Index of the switching variable which determines the composition of the
     *        hydrocarbon phases.
     *
     * Depending on the phases present, this variable is either interpreted as the
     * saturation of the gas phase, as the mole fraction of the gas component in the oil
     * phase or as the mole fraction of the oil component in the gas phase.
     */
    static const int compositionSwitchIdx = PVOffset + 2;

    //! Index of the primary variable for the first solvent
    static const int solventSaturationIdx =
        enableSolvent ? PVOffset + (numPhases-1) + numSolvents_ : -1000;

    //! Index of the primary variable for the first polymer
    static const int polymerConcentrationIdx =
        enablePolymer ? PVOffset + (numPhases-1) + numSolvents_ + numPolymers_ : -1000;

    //! Index of the primary variable for temperature
    static const int temperatureIdx  =
        enableEnergy ? PVOffset + (numPhases-1) + numSolvents_ + numPolymers_ + numEnergy_ : - 1000;

    ////////
    // Equation indices
    ////////

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = PVOffset + 0;
    // two continuity equations follow

    //! Index of the continuity equation for the first solvent component
    static const int contiSolventEqIdx =
        enableSolvent ? PVOffset + (numPhases-1) + numSolvents_ : -1000;

    //! Index of the continuity equation for the first polymer component
    static const int contiPolymerEqIdx =
        enablePolymer > 0 ? PVOffset + (numPhases-1) + numSolvents_ + numPolymers_ : -1000;

    //! Index of the continuity equation for energy
    static const int contiEnergyEqIdx =
        enableEnergy ? PVOffset + (numPhases-1) + numSolvents_ + numPolymers_ + numEnergy_ : -1000;
};

} // namespace Ewoms

#endif
