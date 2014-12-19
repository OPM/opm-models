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
 * \ingroup Newton
 * \brief A Newton solver specific to the NCP model.
 */
template <class TypeTag>
class NcpNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
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
     * \copydoc FvBaseNewtonMethod::update_
     */
    void update_(SolutionVector &uCurrentIter,
                 const SolutionVector &uLastIter,
                 const GlobalEqVector &deltaU,
                 const GlobalEqVector &previousResidual)
    {
        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(deltaU.two_norm2()))
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
                this->error_*EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
            Scalar newtonTol = this->tolerance();

            Scalar relinearizationTol = 0.01*linearTol;
            if (relinearizationTol < newtonTol/10)
                relinearizationTol = newtonTol/10;

            this->model_().linearizer().updateDiscrepancy(previousResidual);
            this->model_().linearizer().computeColors(relinearizationTol);
        }

        // normal Newton-Raphson update
        for (size_t i = 0; i < uLastIter.size(); ++i) {
            for (int j = 0; j < numEq; ++j) {
                uCurrentIter[i][j] = uLastIter[i][j] - deltaU[i][j];
            }
            this->model_().invalidateIntensiveQuantitiesCacheEntry(i, /*timeIdx=*/0);
        }

        // put crash barriers along the update path at the
        // beginning...
        if (this->numIterations_ < choppedIterations_) {
            chopUpdate_(uCurrentIter, uLastIter);
        }
    }

private:
    void chopUpdate_(SolutionVector &uCurrentIter,
                     const SolutionVector &uLastIter)
    {
        for (size_t i = 0; i < uLastIter.size(); ++i) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                saturationChop_(uCurrentIter[i][saturation0Idx + phaseIdx],
                                uLastIter[i][saturation0Idx + phaseIdx]);
            pressureChop_(uCurrentIter[i][pressure0Idx],
                          uLastIter[i][pressure0Idx]);

            // fugacities
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar &val = uCurrentIter[i][fugacity0Idx + compIdx];
                Scalar oldVal = uLastIter[i][fugacity0Idx + compIdx];

                // allow the mole fraction of the component to change
                // at most 70% (assuming composition independent
                // fugacity coefficients)
                Scalar minPhi = this->problem().model().minActivityCoeff(i, compIdx);
                Scalar maxDelta = 0.7 * minPhi;

                clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);

                // do not allow mole fractions lager than 101% or
                // smaller than -1%
                val = std::max(-0.01 * minPhi, val);
                val = std::min(1.01 * minPhi, val);
            }
        }
    }

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
