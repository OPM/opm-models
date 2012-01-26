// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \brief Represents the primary variables used in the single phase box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
#ifndef DUMUX_1P_PRIMARY_VARIABLES_HH
#define DUMUX_1P_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

namespace Dumux
{
/*!
 * \ingroup 1PModel
 *
 * \brief Represents the primary variables used in the M-phase,
 *        N-component box model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class OnePPrimaryVariables 
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;


    typedef typename GET_PROP_TYPE(TypeTag, OnePIndices) Indices;

    // primary variable indices
    enum { pressureIdx = Indices::pressureIdx };



    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) EnergyModule;

public:
    /*!
     * \brief Default constructor
     */
    OnePPrimaryVariables()
        : ParentType()
    {
        Valgrind::SetUndefined(*this);
    };

    /*!
     * \brief Constructor with assignment from scalar
     */
    OnePPrimaryVariables(Scalar value)
        : ParentType(value)
    { };

    /*!
     * \brief Copy constructor
     */
    OnePPrimaryVariables(const OnePPrimaryVariables &value)
        : ParentType(value)
    { };

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
                                bool isInEquilibrium = false)
    {
        assignNaive(fluidState);
    }

    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);
        
        (*this)[pressureIdx] = fluidState.pressure(/*phaseIdx=*/0);
    }
};

} // end namepace

#endif
