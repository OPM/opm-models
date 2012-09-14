// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \brief Represents the primary variables used in the fully-implicit
 *        black-oil model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
#ifndef DUMUX_BLACK_OIL_PRIMARY_VARIABLES_HH
#define DUMUX_BLACK_OIL_PRIMARY_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/material/constraintsolvers/immiscibleflash.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

#include "blackoilproperties.hh"

namespace Dumux
{
/*!
 * \ingroup BlackOilModel
 *
 * \brief Represents the primary variables used by the black-oil model.
 */
template <class TypeTag>
class BlackOilPrimaryVariables 
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq) >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldVector<Scalar, numEq> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };


    static_assert(numPhases == 3, "The black-oil model has three phases!");
    static_assert(numComponents == 3, "The black-oil model has three components!");

public:
    /*!
     * \brief Default constructor
     */
    BlackOilPrimaryVariables()
        : ParentType()
    { Valgrind::SetUndefined(*this); }

    /*!
     * \brief Constructor with assignment from scalar
     */
    BlackOilPrimaryVariables(Scalar value)
        : ParentType(value)
    { }

    /*!
     * \brief Copy constructor
     */
    BlackOilPrimaryVariables(const BlackOilPrimaryVariables &value)
        : ParentType(value)
    { }

    template <class FluidState>
    void assignNaive(const FluidState &fluidState)
    {
        (*this)[pressure0Idx] = fluidState.pressure(/*phaseIdx=*/0);
        (*this)[saturation0Idx] = fluidState.saturation(/*phaseIdx=*/0);
        (*this)[saturation0Idx + 1] = fluidState.saturation(/*phaseIdx=*/1);
    }
};

} // end namepace

#endif
