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
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

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
                             "Number of Newton-Raphson iterations for which the update ought"
                             " to be limited");
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
     * \copydoc FvBaseNewtonMethod::update_
     */
    void update_(SolutionVector &currentSolution,
                 const SolutionVector &previousSolution,
                 const GlobalEqVector &solutionUpdate,
                 const GlobalEqVector &previousResidual)
    {
        auto& assembler = this->model().jacobianAssembler();

        // first, write out the current solution to make convergence
        // analysis possible
        this->writeConvergence_(previousSolution, solutionUpdate);

        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(solutionUpdate.two_norm2()))
            OPM_THROW(Opm::NumericalProblem, "Non-finite update!");

        // compute the degree of freedom and element colors for
        // partial relinearization
        if (this->enablePartialRelinearization_()) {
            // rationale: the change of the derivatives of the
            // residual are relatively small if the solution is
            // largely unchanged and a solution is largly unchanged if
            // the right hand side is close to zero. This argument may
            // not be bullet proof, but it is a heuristic that usually
            // works.
            Scalar linearTol =
                this->error_ * EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
            Scalar newtonTol = this->tolerance();

            Scalar relinearizationTol = 0.01*linearTol;
            if (relinearizationTol < newtonTol/10)
                relinearizationTol = newtonTol/10;

            assembler.updateDiscrepancy(previousResidual);
            assembler.computeColors(relinearizationTol);
        }

        // update the solution
        for (unsigned i = 0; i < previousSolution.size(); ++i) {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                // calculate the update of the curren primary variable. Here we limit the
                // pressure update, but do we not clamp anything after the specified
                // number of iterations was reached
                Scalar delta = solutionUpdate[i][eqIdx];
                if (this->numIterations_ < numChoppedIterations_) {
                    // limit changes in pressure to 10% of the value of the previous iteration
                    if (eqIdx == Indices::gasPressureIdx
                        && std::abs(delta/previousSolution[i][eqIdx]) > 0.1)
                    {
                        delta /= std::abs(delta/(0.1*previousSolution[i][eqIdx]));
                    }
                    // limit changes in saturation to 20%
                    else if ((eqIdx == Indices::waterSaturationIdx ||
                              (eqIdx == Indices::switchIdx
                               && currentSolution[i].switchingVariableIsGasSaturation()))
                             && std::abs(delta) > 0.2)
                    {
                        delta /= std::abs(delta/0.2);
                    }

                }

                // do the actual update
                currentSolution[i][eqIdx] = previousSolution[i][eqIdx] - delta;
            }

            this->model_().invalidateIntensiveQuantitiesCacheEntry(i, /*timeIdx=*/0);
        }
    }

private:
    int numChoppedIterations_;
};
} // namespace Ewoms

#endif
