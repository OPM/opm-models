// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \copydoc Ewoms::NcpNewtonMethod
 */
#ifndef EWOMS_NCP_NEWTON_METHOD_HH
#define EWOMS_NCP_NEWTON_METHOD_HH

#include "ncpproperties.hh"

#include <algorithm>

namespace Ewoms {

/*!
 * \ingroup NcpModel
 *
 * \brief A Newton solver specific to the NCP model.
 */
template <class TypeTag>
class NcpNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };

public:
    /*!
     * \copydoc FvBaseNewtonMethod::FvBaseNewtonMethod(Problem &)
     */
    NcpNewtonMethod(Simulator &simulator) : ParentType(simulator)
    {
        choppedIterations_ = EWOMS_GET_PARAM(TypeTag, int, NcpNewtonNumChoppedIterations);
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);
    }

    /*!
     * \copydoc NewtonMethod::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, int, NcpNewtonNumChoppedIterations,
                             "The number of Newton iterations for which the "
                             "update gets limited");
    }

    // HACK: this is necessary since GCC 4.4 does not support
    // befriending typedefs...
/*
private:
    friend class NewtonMethod<TypeTag>;
    friend class ParentType;
*/
    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(int globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        // normal Newton-Raphson update
        nextValue = currentValue;
        nextValue -= update;

        // put crash barriers along the update path at the
        // beginning...
        if (this->numIterations_ < choppedIterations_) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                saturationChop_(nextValue[saturation0Idx + phaseIdx],
                                currentValue[saturation0Idx + phaseIdx]);
            pressureChop_(nextValue[pressure0Idx],
                          currentValue[pressure0Idx]);

            // fugacities
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar &val = nextValue[fugacity0Idx + compIdx];
                Scalar oldVal = currentValue[fugacity0Idx + compIdx];

                // allow the mole fraction of the component to change
                // at most 70% (assuming composition independent
                // fugacity coefficients)
                Scalar minPhi = this->problem().model().minActivityCoeff(globalDofIdx, compIdx);
                Scalar maxDelta = 0.7 * minPhi;

                clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);

                // do not allow mole fractions lager than 101% or
                // smaller than -1%
                val = std::max<Scalar>(-0.01 * minPhi, val);
                val = std::min<Scalar>(1.01 * minPhi, val);
            }
        }
    }

private:
    void clampValue_(Scalar &val, Scalar minVal, Scalar maxVal) const
    { val = std::max(minVal, std::min(val, maxVal)); }

    void pressureChop_(Scalar &val, Scalar oldVal) const
    {
        // limit pressure updates to 20% per iteration
        clampValue_(val, oldVal * 0.8, oldVal * 1.2);
    }

    void saturationChop_(Scalar &val, Scalar oldVal) const
    {
        // limit saturation updates to 20% per iteration
        const Scalar maxDelta = 0.20;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
    }

    int choppedIterations_;
};
} // namespace Ewoms

#endif
