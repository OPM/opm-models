/*
  Copyright (C) 2008-2014 by Andreas Lauser

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
 * \copydoc Ewoms::BlackOilNewtonMethod
 */
#ifndef EWOMS_BLACK_OIL_NEWTON_METHOD_HH
#define EWOMS_BLACK_OIL_NEWTON_METHOD_HH

#include "blackoilproperties.hh"

namespace Ewoms {

/*!
 * \ingroup Newton
 * \ingroup BlackOilModel
 *
 * \brief A newton solver which is specific to the black oil model.
 */
template <class TypeTag>
class BlackOilNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;

    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    BlackOilNewtonMethod(Simulator &simulator) : ParentType(simulator)
    { numChoppedIterations_ = EWOMS_GET_PARAM(TypeTag, int, BlackoilNumChoppedIterations); }

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, int, BlackoilNumChoppedIterations,
                             "Number of Newton-Raphson iterations for which the update gets"
                             " limited");
    }


    // HACK which is necessary because GCC 4.4 does not support
    // being a friend of typedefs
/*
protected:
    friend class NewtonMethod<TypeTag>;
    friend class ParentType;
*/

    /*!
     * \copydoc FvBaseNewtonMethod::endIteration_
     */
    void endIteration_(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ParentType::endIteration_(uCurrentIter, uLastIter);
        this->problem().model().switchPrimaryVars_();
    }

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(int globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            // calculate the update of the current primary variable. For the
            // black-oil model we limit the pressure and saturation updates, but do
            // we not clamp anything after the specified number of iterations was
            // reached
            Scalar delta = update[eqIdx];
            if (this->numIterations_ < numChoppedIterations_) {
                // limit changes in pressure to 20% of the pressure value at the
                // beginning of the current iteration
                if (eqIdx == Indices::gasPressureIdx
                    && std::abs(delta/currentValue[eqIdx]) > 0.2)
                {
                    delta /= std::abs(delta/(0.2*currentValue[eqIdx]));
                }
                // limit changes in saturation to 20%
                else if ((eqIdx == Indices::waterSaturationIdx ||
                          (eqIdx == Indices::switchIdx
                           && currentValue.switchingVarMeaning() == PrimaryVariables::GasSaturation))
                         && std::abs(delta) > 0.2)
                {
                    delta /= std::abs(delta/0.2);
                }

            }

            // do the actual update
            nextValue[eqIdx] = currentValue[eqIdx] - delta;
        }
    }

private:
    int numChoppedIterations_;
};
} // namespace Ewoms

#endif
