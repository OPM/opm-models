/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Ewoms::RichardsRateVector
 */
#ifndef EWOMS_RICHARDS_RATE_VECTOR_HH
#define EWOMS_RICHARDS_RATE_VECTOR_HH

#include <dune/common/fvector.hh>

#include <opm/material/Valgrind.hpp>
#include <opm/material/constraintsolvers/NcpFlash.hpp>

#include "richardsintensivequantities.hh"

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 *
 * \brief Implements a vector representing mass, molar or volumetric rates.
 *
 * This class is basically a Dune::FieldVector which can be set using either mass, molar
 * or volumetric rates.
 */
template <class TypeTag>
class RichardsRateVector
    : public Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                               GET_PROP_VALUE(TypeTag, NumEq)>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) EnergyModule;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { contiEqIdx = Indices::contiEqIdx };
    enum { liquidCompIdx = GET_PROP_VALUE(TypeTag, LiquidComponentIndex) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dune::FieldVector<Scalar, numEq> ParentType;

public:
    RichardsRateVector() : ParentType()
    { Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(Scalar)
     */
    RichardsRateVector(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::ImmiscibleRateVector(const
     * ImmiscibleRateVector &)
     */
    RichardsRateVector(const RichardsRateVector &value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmiscibleRateVector::setMassRate
     */
    void setMassRate(const ParentType &value)
    { ParentType::operator=(value); }

    /*!
     * \copydoc ImmiscibleRateVector::setMolarRate
     */
    void setMolarRate(const ParentType &value)
    {
        // convert to mass rates
        ParentType massRate(value);
        massRate[contiEqIdx] *= FluidSystem::molarMass(liquidCompIdx);

        // set the mass rate
        setMassRate(massRate);
    }

    /*!
     * \copydoc ImmiscibleRateVector::setEnthalpyRate
     */
    void setEnthalpyRate(Scalar rate)
    { EnergyModule::setEnthalpyRate(*this, rate); }

    /*!
     * \copydoc ImmiscibleRateVector::setVolumetricRate
     */
    template <class FluidState>
    void setVolumetricRate(const FluidState &fluidState, int phaseIdx,
                           Scalar volume)
    {
       (*this)[contiEqIdx] =
            fluidState.density(phaseIdx)
            * fluidState.massFraction(phaseIdx, liquidCompIdx)
            * volume;

        EnergyModule::setEnthalpyRate(*this, fluidState, phaseIdx, volume);
    }
};

} // namespace Ewoms

#endif
