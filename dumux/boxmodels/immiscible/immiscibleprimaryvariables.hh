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
 * \copydoc Dumux::ImmisciblePrimaryVariables
 */
#ifndef DUMUX_IMMISCIBLE_PRIMARY_VARIABLES_HH
#define DUMUX_IMMISCIBLE_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>
#include <dumux/material/constraintsolvers/immiscibleflash.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

#include "immiscibleproperties.hh"

namespace Dumux {

/*!
 * \ingroup ImmiscibleModel
 *
 * \brief Represents the primary variables used by the immiscible
 *        multi-phase, model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class ImmisciblePrimaryVariables
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
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Dumux::ImmiscibleFlash<Scalar, FluidSystem> ImmiscibleFlash;
    typedef BoxMultiPhaseEnergyModule<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergy)> EnergyModule;

public:
    /*!
     * \brief Default constructor
     */
    ImmisciblePrimaryVariables()
        : ParentType()
    {
        Valgrind::SetUndefined(*this);
    }

    /*!
     * \brief Constructor with assignment from scalar
     *
     * \param value The scalar value to which all entries of the vector will be set.
     */
    ImmisciblePrimaryVariables(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \brief Copy constructor
     *
     * \param value The primary variables that will be duplicated.
     */
    ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables &value)
        : ParentType(value)
    { }

    /*!
     * \brief Set the primary variables from an arbitrary fluid state
     *        in a mass conservative way.
     *
     * If an energy equation is included, the fluid temperatures are
     * the same as the one given in the fluid state, *not* the
     * enthalpy.
     *
     * \param fluidState The fluid state which should be represented
     *                   by the primary variables. The temperatures,
     *                   pressures, compositions and densities of all
     *                   phases must be defined.
     * \param matParams The capillary pressure law parameters
     * \param isInEquilibrium If true, the fluid state expresses
     *                        thermodynamic equilibrium assuming the
     *                        relations expressed by the fluid
     *                        system. This implies that in addition to
     *                        the quantities mentioned above, the
     *                        fugacities are also defined.
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
     * \brief Directly retrieve the primary variables from an
     *        arbitrary fluid state.
     *
     * This method retrieves all primary variables from an abitrary
     * fluid state without careing whether the state which is
     * represented by the resulting primary variables features the
     * equivalent mass as the given fluid state. This method is
     * massively cheaper and simpler than assignMassConservative() but
     * it should be used with care!
     *
     * \param fluidState The fluid state which should be represented
     *                   by the primary variables. The temperatures,
     *                   pressures, compositions and densities of all
     *                   phases must be defined.
     */
    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        (*this)[pressure0Idx] = fluidState.pressure(/*phaseIdx=*/0);
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++ phaseIdx)
            (*this)[saturation0Idx + phaseIdx] = fluidState.saturation(phaseIdx);
    }
};

} // end namepace

#endif
