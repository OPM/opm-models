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
 * \copydoc Ewoms::BlackOilTwoPhaseIndices
 */
#ifndef EWOMS_BLACK_OIL_TWO_PHASE_INDICES_HH
#define EWOMS_BLACK_OIL_TWO_PHASE_INDICES_HH

#include <cassert>

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 *
 * \brief The primary variable and equation indices for the black-oil model.
 */
template <unsigned numSolventsV, unsigned numPolymersV, unsigned PVOffset, unsigned disabledCanonicalCompIdx>
struct BlackOilTwoPhaseIndices
{

    //! Is phase enabled or not
    static const bool oilEnabled = disabledCanonicalCompIdx==0? false:true;
    static const bool waterEnabled = disabledCanonicalCompIdx==1? false:true;
    static const bool gasEnabled = disabledCanonicalCompIdx==2? false:true;

    //! Number of phases active at all times
    static const int numPhases = 2;

    //! Number of solvent components considered
    static const int numSolvents = numSolventsV;

    //! Number of polymer components considered
    static const int numPolymers = numPolymersV;

    //! The number of equations
    static const int numEq = numPhases + numSolvents + numPolymers;

    //////////////////////////////
    // Primary variable indices
    //////////////////////////////

    //! The index of the water saturation. For two-phase oil gas models this is disabled.
    static const int waterSaturationIdx  = waterEnabled ? PVOffset + 0: -10000;

    //! Index of the oil pressure in a vector of primary variables
    static const int pressureSwitchIdx  = PVOffset + 1;

    /*!
     * \brief Index of the switching variable which determines the composition of the
     *        hydrocarbon phases.
     *
     * \note For two-phase water oil models this is disabled.
     */
    static const int compositionSwitchIdx = gasEnabled ? PVOffset + 0: -10000;

    //! Index of the primary variable for the first solvent
    static const int solventSaturationIdx  = PVOffset + numPhases;

    //! Index of the primary variable for the first polymer
    static const int polymerConcentrationIdx  = solventSaturationIdx + numPolymers;

    // numSolvents-1 primary variables follow


    //////////////////////
    // Equation indices
    //////////////////////
    //! \brief returns the index of "active" component
    static unsigned canonicalToActiveComponentIndex(unsigned compIdx)
    {
        // assumes canonical oil = 0, water = 1, gas = 2;
        if(!gasEnabled) {
            assert(compIdx != 2);
            // oil = 0, water = 1
            return compIdx;
        } else if (!waterEnabled) {
            assert(compIdx != 1);
            // oil = 0, gas = 1
            return compIdx / 2;
        } else {
            assert(!oilEnabled);
            assert(compIdx != 0);
        }
        // water = 0, gas = 1;
        return compIdx-1;
    }

    static unsigned activeToCanonicalComponentIndex(unsigned compIdx)
    {
        // assumes canonical oil = 0, water = 1, gas = 2;
        assert(compIdx < 2);
        if(!gasEnabled) {
            // oil = 0, water = 1
            return compIdx;
        } else if (!waterEnabled) {
            // oil = 0, gas = 1
            return compIdx * 2;
        } else {
            assert(!oilEnabled);
        }
        // water = 0, gas = 1;
        return compIdx+1;
    }

    //! Index of the continuity equation of the first phase
    static const int conti0EqIdx = PVOffset + 0;
    // two continuity equations follow

    //! Index of the continuity equation for the first solvent component
    static const int contiSolventEqIdx = PVOffset + numPhases - 1 + numSolvents;

    //! Index of the continuity equation for the first polymer component
    static const int contiPolymerEqIdx = contiSolventEqIdx + numPolymers;
    // numSolvents-1 continuity equations follow

};

} // namespace Ewoms

#endif
