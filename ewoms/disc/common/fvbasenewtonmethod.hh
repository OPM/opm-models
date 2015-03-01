/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \copydoc Ewoms::FvBaseNewtonMethod
 */
#ifndef EWOMS_FV_BASE_NEWTON_METHOD_HH
#define EWOMS_FV_BASE_NEWTON_METHOD_HH

#include "fvbasenewtonconvergencewriter.hh"

#include <ewoms/nonlinear/newtonmethod.hh>
#include <opm/core/utility/PropertySystem.hpp>

namespace Ewoms {

template <class TypeTag>
class FvBaseNewtonMethod;

template <class TypeTag>
class FvBaseNewtonConvergenceWriter;
} // namespace Ewoms

namespace Opm {
namespace Properties {
//! create a type tag for the Newton method of the finite-volume discretization
NEW_TYPE_TAG(FvBaseNewtonMethod, INHERITS_FROM(NewtonMethod));

//! The class dealing with the balance equations
NEW_PROP_TAG(Model);

//! The class storing primary variables plus pseudo primary variables
NEW_PROP_TAG(PrimaryVariables);

//! The number of balance equations.
NEW_PROP_TAG(NumEq);

//! The discretization specific part of he implementing the Newton algorithm
NEW_PROP_TAG(DiscNewtonMethod);

//! The class implementing the Newton algorithm
NEW_PROP_TAG(NewtonMethod);

//! Specifies whether the linearization should only be relinearized if
//! the current solution deviates too much from the evaluation point
NEW_PROP_TAG(EnablePartialRelinearization);

//! some user adjustable knob to increase or decrease the reliniearization tolerance.
NEW_PROP_TAG(RelinearizationToleranceFactor);

//! Enable linearization recycling?
NEW_PROP_TAG(EnableLinearizationRecycling);

/*!
 * \brief The desired residual reduction of the linear solver.
 *
 * In the context of this class, this is needed to implement partial
 * relinearization of the linearized system of equations.
 */
NEW_PROP_TAG(LinearSolverTolerance);

// set default values
SET_TYPE_PROP(FvBaseNewtonMethod, DiscNewtonMethod,
              Ewoms::FvBaseNewtonMethod<TypeTag>);
SET_TYPE_PROP(FvBaseNewtonMethod, NewtonMethod,
              typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod));
SET_TYPE_PROP(FvBaseNewtonMethod, NewtonConvergenceWriter,
              Ewoms::FvBaseNewtonConvergenceWriter<TypeTag>);
}} // namespace Properties, Opm

namespace Ewoms {
/*!
 * \ingroup Discretization
 * \ingroup Newton
 *
 * \brief A Newton method for models using a finite volume discretization.
 *
 * This class is sufficient for most models which use an Element or a
 * Vertex Centered Finite Volume discretization.
 */
template <class TypeTag>
class FvBaseNewtonMethod : public NewtonMethod<TypeTag>
{
    typedef Ewoms::NewtonMethod<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;


public:
    FvBaseNewtonMethod(Simulator &simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters of the Newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, RelinearizationToleranceFactor,
                             "The scaling factor from the default tolerance used for partial "
                             "relinearization");
    }

protected:
    friend class Ewoms::NewtonMethod<TypeTag>;

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the converged() and
     * proceed() methods should be updated inside this method.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param nextSolution The solution vector at the end of the current iteration
     * \param currentSolution The solution vector at the beginning of the current iteration
     * \param solutionUpdate The delta as calculated by solving the linear system of
     *                       equations. This parameter also stores the updated solution.
     * \param currentResidual The residual (i.e., right-hand-side) of the current solution.
     */
    void update_(SolutionVector &nextSolution,
                 const SolutionVector &currentSolution,
                 const GlobalEqVector &solutionUpdate,
                 const GlobalEqVector &currentResidual)
    {
        auto& model = this->model();
        auto& linearizer = model.linearizer();

        // first, write out the next solution to make convergence
        // analysis possible
        this->writeConvergence_(currentSolution, solutionUpdate);

        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(solutionUpdate.two_norm2()))
            OPM_THROW(Opm::NumericalProblem, "Non-finite update!");

        // compute the degree of freedom and element colors for
        // partial relinearization
        if (enablePartialRelinearization_()) {
            // rationale: the change of the derivatives of the
            // residual are relatively small if the solution is
            // largely unchanged and a solution is largly unchanged if
            // the right hand side is close to zero. This argument may
            // not be bullet proof, but it is a heuristic that usually
            // works.
            Scalar linearTol =
                this->error_ * EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
            Scalar newtonTol = this->tolerance();

            Scalar relinearizationTol = linearTol/500;
            if (relinearizationTol < newtonTol/100)
                relinearizationTol = newtonTol/100;

            relinearizationTol *= EWOMS_GET_PARAM(TypeTag, Scalar, RelinearizationToleranceFactor);

            linearizer.updateDiscrepancy(currentResidual);
            linearizer.computeColors(relinearizationTol);
        }

        // update the solution for the grid DOFs
        for (unsigned dofIdx = 0; dofIdx < model.numGridDof(); ++dofIdx) {
            asImp_().updatePrimaryVariables_(dofIdx,
                                             nextSolution[dofIdx],
                                             currentSolution[dofIdx],
                                             solutionUpdate[dofIdx],
                                             currentResidual[dofIdx]);
            model.invalidateIntensiveQuantitiesCacheEntry(dofIdx, /*timeIdx=*/0);
        }

        // update the DOFs of the auxiliary equations
        int nTotalDof = model.numTotalDof();
        for (int dofIdx = model.numGridDof(); dofIdx < nTotalDof; ++dofIdx) {
            nextSolution[dofIdx] = currentSolution[dofIdx];
            nextSolution[dofIdx] -= solutionUpdate[dofIdx];
        }
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void beginIteration_()
    {
        model_().syncOverlap();

        ParentType::beginIteration_();
    }

    /*!
     * \brief Called if the Newton method broke down.
     */
    void failed_()
    {
        ParentType::failed_();

        model_().linearizer().relinearizeAll();
    }

    /*!
     * \brief Called when the Newton method was successful.
     */
    void succeeded_()
    {
        ParentType::succeeded_();

        if (enableLinearizationRecycling_())
            model_().linearizer().setLinearizationReusable(true);
        else
            model_().linearizer().relinearizeAll();
    }

    /*!
     * \brief Returns a reference to the model.
     */
    Model &model_()
    { return ParentType::model(); }

    /*!
     * \brief Returns a reference to the model.
     */
    const Model &model_() const
    { return ParentType::model(); }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    /*!
     * \brief Returns true iff the linearizer uses partial
     *        relinearization.
     */
    bool enablePartialRelinearization_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, EnablePartialRelinearization); }

    /*!
     * \brief Returns true iff the linearizer recycles the
     *        linearization if possible.
     */
    bool enableLinearizationRecycling_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableLinearizationRecycling); }
};
} // namespace Ewoms

#endif
