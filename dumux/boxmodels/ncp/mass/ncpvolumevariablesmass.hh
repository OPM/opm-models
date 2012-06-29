// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Contains the mass conservation part of the volume variables
 */
#ifndef DUMUX_NCP_VOLUME_VARIABLES_MASS_HH
#define DUMUX_NCP_VOLUME_VARIABLES_MASS_HH

#include <dumux/boxmodels/ncp/ncpproperties.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include <dune/common/fvector.hh>

#include <sstream>
#include <string>

namespace Dumux
{
/*!
 * \brief The compositional part of the volume variables if chemical
 *        equilibrium is assumed
 */
template <class TypeTag>
class NcpVolumeVariablesMass
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, CompositionFromFugacitiesSolver) CompositionFromFugacitiesSolver;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { fugacityOverPressure0Idx = Indices::fugacityOverPressure0Idx };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef typename FluidSystem::ParameterCache ParameterCache;

public:
    /*!
     * \brief The fluid state which is used by the volume variables to
     *        store the thermodynamic state.
     *
     * If chemical equilibrium is assumed, we use the fluid state
     * which saves some memory.
     */
    typedef CompositionalFluidState<Scalar, FluidSystem, /*enableEnthalpy=*/GET_PROP_VALUE(TypeTag, EnableEnergy)> FluidState;

    /*!
     * \brief Update composition of all phases in the mutable
     *        parameters from the primary variables.
     */
    void update(FluidState &fs,
                ParameterCache &paramCache,
                const ElementContext &elemCtx,
                int scvIdx,
                int timeIdx)
    {
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        const auto *hint = elemCtx.hint(scvIdx, timeIdx);


        // calculate phase compositions
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // initial guess
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar xIJ = 1.0/numComponents;
                if (hint)
                    // use the hint for the initial mole fraction!
                    xIJ = hint->fluidState().moleFraction(phaseIdx, compIdx);

                // set initial guess of the component's mole fraction
                fs.setMoleFraction(phaseIdx, compIdx, xIJ);

            }

            ComponentVector fug;
            // retrieve component fugacities
            Scalar pAlpha = fs.pressure(phaseIdx);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fug[compIdx] = priVars[fugacityOverPressure0Idx + compIdx] * pAlpha;

            // calculate the phase composition from the component
            // fugacities
            if (!hint)
                // if no hint was given, we ask the constraint solver
                // to make an initial guess
                CompositionFromFugacitiesSolver::guessInitial(fs, paramCache, phaseIdx, fug);
            CompositionFromFugacitiesSolver::solve(fs, paramCache, phaseIdx, fug);

            /*
              std::cout << "final composition: " << FluidSystem::phaseName(phaseIdx) << "="
              << fs.moleFraction(phaseIdx, 0) << " "
              << fs.moleFraction(phaseIdx, 1) << " "
              << fs.moleFraction(phaseIdx, 2) << " "
              << fs.moleFraction(phaseIdx, 3) << " "
              << fs.moleFraction(phaseIdx, 4) << " "
              << fs.moleFraction(phaseIdx, 5) << " "
              << fs.moleFraction(phaseIdx, 6) << "\n";
            */

        }
    }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
    }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    static std::string primaryVarName(int pvIdx)
    {
        std::ostringstream oss;
        if (Indices::fugacityOverPressure0Idx <= pvIdx && pvIdx < Indices::fugacityOverPressure0Idx + numComponents)
            oss << "f^" << FluidSystem::componentName(pvIdx - Indices::fugacityOverPressure0Idx) << "/p";
        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    static std::string eqName(int eqIdx)
    {
        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "continuity^" << FluidSystem::componentName(eqIdx - Indices::conti0EqIdx);
        return oss.str();
    }
};

} // end namepace

#endif
