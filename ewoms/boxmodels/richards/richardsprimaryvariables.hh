// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::RichardsPrimaryVariables
 */
#ifndef EWOMS_RICHARDS_PRIMARY_VARIABLES_HH
#define EWOMS_RICHARDS_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <ewoms/material/constraintsolvers/immiscibleflash.hh>
#include <ewoms/material/fluidstates/immisciblefluidstate.hh>

#include "richardsproperties.hh"

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 *
 * \brief Represents the primary variables used in the Richards box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class RichardsPrimaryVariables
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    // primary variable indices
    enum { pwIdx = Indices::pwIdx };

    enum {
        wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex),
        nPhaseIdx = 1 - wPhaseIdx
    };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) EnergyModule;

    typedef Ewoms::ImmiscibleFlash<Scalar, FluidSystem> ImmiscibleFlash;

public:
    RichardsPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetUndefined(*this);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    RichardsPrimaryVariables(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables &)
     */
    RichardsPrimaryVariables(const RichardsPrimaryVariables &value)
        : ParentType(value)
    { }

    /*!
     * \brief Set the primary variables with the wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pw The pressure of the wetting phase [Pa]
     * \param Sw The saturation of the wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    void assignImmiscibleFromWetting(Scalar T,
                                     Scalar pw,
                                     Scalar Sw,
                                     const MaterialLawParams &matParams)
    {
        ImmiscibleFluidState<Scalar, FluidSystem> fs;

        fs.setTemperature(T);
        fs.setSaturation(wPhaseIdx, Sw);
        fs.setSaturation(nPhaseIdx, 1 - Sw);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(wPhaseIdx, pw);
        fs.setPressure(nPhaseIdx, pw + (pC[nPhaseIdx] - pC[wPhaseIdx]));

        assignNaive(fs);
    }

    /*!
     * \brief Set the primary variables with the non-wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pn The pressure of the non-wetting phase [Pa]
     * \param Sn The saturation of the non-wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    void assignImmiscibleFromNonWetting(Scalar T,
                                        Scalar pn,
                                        Scalar Sn,
                                        const MaterialLawParams &matParams)
    {
        ImmiscibleFluidState<Scalar, FluidSystem> fs;

        fs.setTemperature(T);
        fs.setSaturation(wPhaseIdx, 1 - Sn);
        fs.setSaturation(nPhaseIdx, Sn);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(nPhaseIdx, pn);
        fs.setPressure(nPhaseIdx, pn + (pC[wPhaseIdx] - pC[nPhaseIdx]));

        assignNaive(fs);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState &fluidState,
                                const MaterialLawParams &matParams,
                                bool isInEquilibrium = false)
    {
        ComponentVector globalMolarities(0.0);
        for (int compIdx = 0; compIdx < numComponents; ++ compIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                globalMolarities[compIdx] +=
                    fluidState.molarity(phaseIdx, compIdx)
                    * fluidState.saturation(phaseIdx);
            }
        }

        ImmiscibleFluidState<Scalar, FluidSystem> fsFlash;
        fsFlash.assign(fluidState);
        typename FluidSystem::ParameterCache paramCache;
        ImmiscibleFlash::template solve<MaterialLaw>(fsFlash, paramCache, matParams, globalMolarities);

        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        (*this)[pwIdx] = fluidState.pressure(wPhaseIdx);
    }
};

} // end namepace

#endif
