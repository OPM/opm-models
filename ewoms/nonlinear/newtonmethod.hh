// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2013 by Andreas Lauser                               *
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
 * \copydoc Ewoms::NewtonMethod
 */
#ifndef EWOMS_NEWTON_METHOD_HH
#define EWOMS_NEWTON_METHOD_HH

#include "nullconvergencewriter.hh"

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/PropertySystem.hpp>
#include <opm/core/utility/ClassName.hpp>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/parallel/mpihelper.hh>

#include <dune/common/timer.hh>
#include <dune/istl/istlexception.hh>

#include <iostream>
#include <sstream>

namespace Ewoms {
// forward declaration of classes
template <class TypeTag>
class NewtonMethod;
}

namespace Opm {
// forward declaration of property tags
namespace Properties {
//! The type tag on which the default properties for the Newton method
//! are attached
NEW_TYPE_TAG(NewtonMethod);

//! The physical model which we would like to solve
NEW_PROP_TAG(Problem);

//! The model describing the PDEs for the conservation quantities
NEW_PROP_TAG(Model);

//! The type of scalar values
NEW_PROP_TAG(Scalar);

//! Specifies the type of the actual Newton method
NEW_PROP_TAG(NewtonMethod);

//! Specifies the type of a solution
NEW_PROP_TAG(SolutionVector);

//! Vector containing a quantity of for equation on the whole grid
NEW_PROP_TAG(GlobalEqVector);

//! Specifies the class of the physical problem
NEW_PROP_TAG(Problem);

//! The class which linearizes the non-linear system of equations
NEW_PROP_TAG(JacobianAssembler);

//! Specifies the type of a global Jacobian matrix
NEW_PROP_TAG(JacobianMatrix);

//! Specifies the type of the linear solver to be used
NEW_PROP_TAG(LinearSolverBackend);

//! Specifies whether the Newton method should print messages or not
NEW_PROP_TAG(NewtonVerbose);

//! Specifies the type of the class which writes out the Newton convergence
NEW_PROP_TAG(NewtonConvergenceWriter);

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
NEW_PROP_TAG(NewtonWriteConvergence);

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
NEW_PROP_TAG(ConvergenceWriter);

//! Indicate whether the relative error should be used
NEW_PROP_TAG(NewtonEnableRelativeCriterion);

//! The value for the relative error below which convergence is
//! declared
NEW_PROP_TAG(NewtonRelativeTolerance);

//! The maximum relative error which may occur in a simulation before
//! the Newton method is aborted
NEW_PROP_TAG(NewtonMaxRelativeError);

//! Indicate whether the absolute error should be used
NEW_PROP_TAG(NewtonEnableAbsoluteCriterion);

//! The value for the absolute error reduction below which convergence
//! is declared
NEW_PROP_TAG(NewtonAbsoluteTolerance);

//! Indicate whether both of the criteria should be satisfied to
//! declare convergence
NEW_PROP_TAG(NewtonSatisfyAbsoluteAndRelative);

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time-step size. The heuristic used
 * is to scale the last time-step size by the deviation of the
 * number of iterations used from the target steps.
 */
NEW_PROP_TAG(NewtonTargetIterations);

//! Number of maximum iterations for the Newton method.
NEW_PROP_TAG(NewtonMaxIterations);

// set default values for the properties
SET_TYPE_PROP(NewtonMethod, NewtonMethod, Ewoms::NewtonMethod<TypeTag>);
SET_TYPE_PROP(NewtonMethod, NewtonConvergenceWriter, Ewoms::NullConvergenceWriter<TypeTag>);
SET_BOOL_PROP(NewtonMethod, NewtonWriteConvergence, false);
SET_BOOL_PROP(NewtonMethod, NewtonEnableRelativeCriterion, true);
SET_BOOL_PROP(NewtonMethod, NewtonEnableAbsoluteCriterion, false);
SET_BOOL_PROP(NewtonMethod, NewtonSatisfyAbsoluteAndRelative, false);
SET_BOOL_PROP(NewtonMethod, NewtonVerbose, true);
SET_SCALAR_PROP(NewtonMethod, NewtonRelativeTolerance, 1e-8);
SET_SCALAR_PROP(NewtonMethod, NewtonAbsoluteTolerance, 1e-5);
SET_SCALAR_PROP(NewtonMethod, NewtonMaxRelativeError, 1e100); // effectively disabled if not overwritten at run-time
SET_INT_PROP(NewtonMethod, NewtonTargetIterations, 10);
SET_INT_PROP(NewtonMethod, NewtonMaxIterations, 18);
}
}

namespace Ewoms {
/*!
 * \ingroup Newton
 * \brief The multi-dimensional Newton method.
 *
 * This class uses static polymorphism to allow implementations to
 * implement different update/convergence strategies.
 */
template <class TypeTag>
class NewtonMethod
{
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) LinearSolverBackend;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonConvergenceWriter) ConvergenceWriter;

    typedef typename Dune::MPIHelper::MPICommunicator Communicator;
    typedef Dune::CollectiveCommunication<Communicator> CollectiveCommunication;

public:
    NewtonMethod(Problem &problem)
        : problem_(problem)
        , endIterMsgStream_(std::ostringstream::out)
        , linearSolver_(problem)
        , comm_(Dune::MPIHelper::getCommunicator())
        , convergenceWriter_(asImp_())
    {
        if (!enableRelativeCriterion_() && !enableAbsoluteCriterion_())
        {
            OPM_THROW(std::logic_error,
                       "At least one of NewtonEnableRelativeCriterion or "
                       "NewtonEnableAbsoluteCriterion has to be set to true");
        }
/*
        else if (satisfyAbsAndRel_() &&
                 (!enableRelativeCriterion_() || !enableAbsoluteCriterion_()))
        {
            OPM_THROW(std::logic_error,
                       "If you set NewtonSatisfyAbsoluteAndRelative to true, you also must set "
                       "NewtonEnableRelativeCriterion and NewtonEnableAbsoluteCriterion");

        }
*/

        lastRelError_ = 1e100;
        lastAbsError_ = 1e100;

        relError_ = 1e100;
        absError_ = 1e100;

        numIterations_ = 0;
    }

    /*!
     * \brief Register all run-time parameters for the Newton method.
     */
    static void registerParameters()
    {
        LinearSolverBackend::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonVerbose, "Specify whether the Newton method should inform the user about its progress or not");
        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonWriteConvergence, "Write the convergence behaviour of the Newton method to a VTK file");
        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonEnableRelativeCriterion, "Make the Newton method consider the relative convergence criterion");
        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonEnableAbsoluteCriterion, "Make the Newton method consider the absolute convergence criterion");
        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonSatisfyAbsoluteAndRelative, "Let the Newton method only consider a solution to be converged if the relative _and_ the absolute criterion are fulfilled");
        EWOMS_REGISTER_PARAM(TypeTag, int, NewtonTargetIterations, "The 'optimimum' number of Newton iterations per time step");
        EWOMS_REGISTER_PARAM(TypeTag, int, NewtonMaxIterations, "The maximum number of Newton iterations per time step");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonRelativeTolerance, "The maximum relative error between two iterations tolerated by the Newton method");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonMaxRelativeError, "The maximum relative error for which the next iteration of the Newton method is executed");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonAbsoluteTolerance, "The maximum residual tolerated by the Newton method");
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool converged() const
    {
        if (enableRelativeCriterion_() && !enableAbsoluteCriterion_())
            // only look at the relative error
            return relError_ <= relTolerance_();
        else if (!enableRelativeCriterion_() && enableAbsoluteCriterion_())
            // only look at the absolute error
            return absError_ <= absTolerance_();
        else if (satisfyAbsAndRel_())
            // both, the absolute and the relative tolerances must be
            // attained
            return
                relError_ <= relTolerance_()
                && absError_ <= absTolerance_();

        // we're done as soon as either the absolute or the
        // relative tolerance is achieved.
        return relError_ <= relTolerance_()
            || absError_ <= absTolerance_();
    }

    /*!
     * \brief Returns a reference to the current numeric problem.
     */
    Problem &problem()
    { return problem_; }

    /*!
     * \brief Returns a reference to the current numeric problem.
     */
    const Problem &problem() const
    { return problem_; }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    Model &model()
    { return problem().model(); }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    const Model &model() const
    { return problem().model(); }

    /*!
     * \brief Returns the number of iterations done since the Newton method
     *        was invoked.
     */
    int numIterations() const
    { return numIterations_; }

    /*!
     * \brief Run the Newton method.
     *
     * The implementation is responsible
     * for all the strategic decisions.
     */
    bool apply()
    {
        SolutionVector &uCurrentIter = model().solution(/*historyIdx=*/0);
        SolutionVector uLastIter(uCurrentIter);
        GlobalEqVector deltaU(uCurrentIter.size());

        JacobianAssembler &jacobianAsm = model().jacobianAssembler();

        // tell the implementation that we begin solving
        asImp_().begin_(uCurrentIter);

        assembleTimer_.reset();
        solveTimer_.reset();
        updateTimer_.reset();

        // execute the method as long as the implementation thinks
        // that we should do another iteration
        while (asImp_().proceed_())
        {
            // notify the implementation that we're about to start
            // a new iteration
            asImp_().beginIteration_();

            // make the current solution to the old one
            uLastIter = uCurrentIter;

            if (asImp_().verbose_()) {
                std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad r";
                std::cout.flush();
            }

            ///////////////
            // assemble
            ///////////////

            // linearize the problem at the current solution
            try {
                assembleTimer_.start();
                asImp_().linearize_();
                assembleTimer_.stop();
            }
            catch (const Dune::Exception &e) {
                if (asImp_().verbose_())
                    std::cout << "Newton: Caught exception during linearization: \"" << e.what() << "\"\n";
                asImp_().failed_();
                return false;
            };

            ///////////////
            // linear solve
            ///////////////

            // Clear the current line using an ansi escape
            // sequence.  for an explanation see
            // http://en.wikipedia.org/wiki/ANSI_escape_code
            const char clearRemainingLine[] = { 0x1b, '[', 'K', 0 };

            if (asImp_().verbose_()) {
                std::cout << "\rSolve: M deltax^k = r";
                std::cout << clearRemainingLine;
                std::cout.flush();
            }

            // solve the resulting linear equation system
            solveTimer_.start();

            // set the delta vector to zero before solving the linear system!
            deltaU = 0;
            // ask the implementation to solve the linearized system
            if (!asImp_().solveLinear_(jacobianAsm.matrix(),
                                       deltaU,
                                       jacobianAsm.residual()))
            {
                if (asImp_().verbose_())
                    std::cout << "Newton: Linear solver did not converge\n";
                solveTimer_.stop();
                asImp_().failed_();
                return false;
            }
            solveTimer_.stop();

            ///////////////
            // update
            ///////////////
            if (asImp_().verbose_()) {
                std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
                std::cout << clearRemainingLine;
                std::cout.flush();
            }

            try {
                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                updateTimer_.start();
                asImp_().updateErrors_(uCurrentIter, uLastIter, deltaU);
                asImp_().update_(uCurrentIter, uLastIter, deltaU);
                updateTimer_.stop();

                // tell the implementation that we're done with this
                // iteration
                asImp_().endIteration_(uCurrentIter, uLastIter);
            }
            catch (const Dune::Exception &e) {
                updateTimer_.stop();
                if (asImp_().verbose_())
                    std::cout << "Newton: Caught exception during update: \"" << e.what() << "\"\n";
                asImp_().failed_();
                return false;
            };

        }

        // tell the implementation that we're done
        asImp_().end_();

        if (asImp_().verbose_()) {
            Scalar elapsedTot = assembleTimer_.elapsed() + solveTimer_.elapsed() + updateTimer_.elapsed();
            std::cout << "Assemble/solve/update time: "
                      <<  assembleTimer_.elapsed() << "(" << 100*assembleTimer_.elapsed()/elapsedTot << "%)/"
                      <<  solveTimer_.elapsed() << "(" << 100*solveTimer_.elapsed()/elapsedTot << "%)/"
                      <<  updateTimer_.elapsed() << "(" << 100*updateTimer_.elapsed()/elapsedTot << "%)"
                      << "\n";
        }

        if (!asImp_().converged()) {
            asImp_().failed_();
            return false;
        }

        asImp_().succeeded_();
        return true;
    }

    /*!
     * \brief Returns the wall time spend so far for linearizing the
     *        non-linear system for all iterations of the current time
     *        step.
     */
    Scalar assembleTime() const
    { return assembleTimer_.elapsed(); }

    /*!
     * \brief Returns the wall time spend so far for solving the
     *        linear systems for all iterations of the current time
     *        step.
     */
    Scalar solveTime() const
    { return solveTimer_.elapsed(); }

    /*!
     * \brief Returns the wall time spend so far for updating the
     *        iterative solutions of the non-linear system for all
     *        iterations of the current time step.
     */
    Scalar updateTime() const
    { return updateTimer_.elapsed();}

    /*!
     * \brief Suggest a new time-step size based on the old time-step
     *        size.
     *
     * The default behavior is to suggest the old time-step size
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time-step.
     */
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // be aggressive reducing the time-step size but
        // conservative when increasing it. the rationale is
        // that we want to avoid failing in the next time
        // integration which would be quite expensive
        if (numIterations_ > targetIterations_()) {
            Scalar percent = Scalar(numIterations_ - targetIterations_())/targetIterations_();
            return oldTimeStep/(1.0 + percent);
        }

        Scalar percent = Scalar(targetIterations_() - numIterations_)/targetIterations_();
        return oldTimeStep*(1.0 + percent/1.2);
    }

    /*!
     * \brief Message that should be printed for the user after the
     *        end of an iteration.
     */
    std::ostringstream &endIterMsg()
    { return endIterMsgStream_; }

protected:
    /*!
     * \brief Returns true if the Newton method ought to be chatty.
     */
    bool verbose_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, NewtonVerbose) && comm_.rank() == 0; }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    void begin_(const SolutionVector &u)
    {
        numIterations_ = 0;

        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence))
            convergenceWriter_.beginTimeStep();
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void beginIteration_()
    {
        lastRelError_ = relError_;
        lastAbsError_ = absError_;
    }

    /*!
     * \brief Assemble the global linear system of equations.
     */
    void linearize_()
    { model().jacobianAssembler().assemble(); }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Opm::NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    bool solveLinear_(const JacobianMatrix &A,
                      GlobalEqVector &x,
                      const GlobalEqVector &b)
    { return linearSolver_.solve(A, x, b);  }

    /*!
     * \brief Update the relative and absolute errors given delta vector.
     *
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void updateErrors_(const SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter,
                       const GlobalEqVector &deltaU)
    {
        lastRelError_ = relError_;
        lastAbsError_ = absError_;
        asImp_().updateRelError_(uCurrentIter, uLastIter, deltaU);
        asImp_().updateAbsError_(uCurrentIter, uLastIter, deltaU);
        writeConvergence_(uLastIter, deltaU);
    }

    /*!
     * \brief Update the relative error of the solution compared to
     *        the previous iteration.
     *
     * The relative error can be seen as a norm of the difference
     * between the current and the next iteration. By default, this is
     * the maxiumum of the difference.
     *
     * \param uCurrentIter The current iterative solution
     * \param uLastIter The last iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    void updateRelError_(const SolutionVector &uCurrentIter,
                         const SolutionVector &uLastIter,
                         const GlobalEqVector &deltaU)
    {
        if (!this->enableRelativeCriterion_)
            return;

        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        relError_ = 0;
        for (unsigned i = 0; i < uLastIter.size(); ++i) {
            const auto &delta = deltaU[i];
            for (unsigned j = 0; j < delta.size(); ++j)
                relError_ = std::max(std::abs(deltaU[i][j]), relError_);
        }

        // take the other processes into account
        relError_ = comm_.max(relError_);

        if (relError_ > EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxRelativeError))
            OPM_THROW(Opm::NumericalProblem,
                      "Newton: Relative error " << relError_
                      << " is larger than maximum allowed error of "
                      << EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxRelativeError));
    }

    /*!
     * \brief Update  the absolute error  of the solution  compared to
     *        the previous iteration.
     *
     * Since this method might be very computationally expensive, the
     * generic Newton method does not really implement it and only
     * throws an exception.
     */
    void updateAbsError_(const SolutionVector &uCurrentIter,
                         const SolutionVector &uLastIter,
                         const GlobalEqVector &deltaU)
    {
        OPM_THROW(std::logic_error,
                  "Not implemented: The Newton controller (" << Opm::className<Implementation>() << ") "
                  "does not provide a newtonUpdateAbsError() method!");
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * Different update strategies, such as line search and chopped
     * updates can be implemented. The default behavior is just to
     * subtract deltaU from uLastIter, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution vector after the current iteration
     * \param uLastIter The solution vector after the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void update_(SolutionVector &uCurrentIter,
                 const SolutionVector &uLastIter,
                 const GlobalEqVector &deltaU)
    {
        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(deltaU.two_norm2()))
            OPM_THROW(Opm::NumericalProblem, "Non-finite update!");

        for (unsigned i = 0; i < uLastIter.size(); ++i) {
            uCurrentIter[i] = uLastIter[i];
            uCurrentIter[i] -= deltaU[i];
        }
    }

    /*!
     * \brief Write the convergence behaviour of the newton method to
     *        disk.
     *
     * This method is called as part of the update proceedure.
     */
    void writeConvergence_(const SolutionVector &uLastIter,
                           const GlobalEqVector &deltaU)
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence)) {
            convergenceWriter_.beginIteration();
            convergenceWriter_.writeFields(uLastIter, deltaU);
            convergenceWriter_.endIteration();
        }
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void endIteration_(const SolutionVector &uCurrentIter,
                  const SolutionVector &uLastIter)
    {
        ++numIterations_;

        if (asImp_().verbose_())
        {
            std::cout << "\rNewton iteration " << numIterations_ << " done";
            if (enableRelativeCriterion_())
                std::cout << ", relative error = " << relError_;
            if (enableAbsoluteCriterion_())
                std::cout << ", absolute error = " << absError_;
            std::cout << endIterMsg().str() << "\n";
        }

        endIterMsgStream_.str("");
    }

    /*!
     * \brief Returns true iff another Newton iteration should be done.
     */
    bool proceed_() const
    {
        if (asImp_().numIterations() < 1)
            return true; // we always do at least one iteration
        else if (asImp_().converged()) {
            // we are below the specified tolerance, so we don't have to
            // do more iterations
            return false;
        }
        else if (asImp_().numIterations() >= asImp_().maxIterations_()) {
            // we have exceeded the allowed number of steps.  If the
            // error was reduced by a factor of at least 4,
            // in the last iterations we proceed even if we are above
            // the maximum number of steps
            if (enableRelativeCriterion_())
                return relError_*4.0 < lastRelError_;
            else
                return absError_*4.0 < lastAbsError_;
        }

        return true;
    }

    /*!
     * \brief Indicates that we're done solving the non-linear system
     *        of equations.
     */
    void end_()
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence))
            convergenceWriter_.endTimeStep();
    }

    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ end_()
     */
    void failed_()
    { numIterations_ = targetIterations_()*2; }

    /*!
     * \brief Called if the Newton method was successful.
     *
     * This method is called _after_ end_()
     */
    void succeeded_()
    { }

    Problem &problem_;

    Dune::Timer assembleTimer_;
    Dune::Timer solveTimer_;
    Dune::Timer updateTimer_;

    std::ostringstream endIterMsgStream_;

    // relative errors and tolerance
    Scalar relError_;
    Scalar lastRelError_;
    Scalar relTolerance_() const
    { return EWOMS_GET_PARAM(TypeTag, Scalar, NewtonRelativeTolerance); }

    // absolute errors and tolerance
    Scalar absError_;
    Scalar lastAbsError_;
    Scalar absTolerance_() const
    { return EWOMS_GET_PARAM(TypeTag, Scalar, NewtonAbsoluteTolerance); }

    // which criteria do we have to satisfy in order to be converged?
    bool enableRelativeCriterion_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, NewtonEnableRelativeCriterion); }
    bool enableAbsoluteCriterion_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, NewtonEnableAbsoluteCriterion); }
    bool satisfyAbsAndRel_() const
    { return EWOMS_GET_PARAM(TypeTag, bool, NewtonSatisfyAbsoluteAndRelative); }

    // optimal number of iterations we want to achieve
    int targetIterations_() const
    { return EWOMS_GET_PARAM(TypeTag, int, NewtonTargetIterations); }
    // maximum number of iterations we do before giving up
    int maxIterations_() const
    { return EWOMS_GET_PARAM(TypeTag, int, NewtonMaxIterations); }
    // actual number of iterations done so far
    int numIterations_;

    // the linear solver
    LinearSolverBackend linearSolver_;

    // the collective communication used by the simulation (i.e. fake
    // or MPI)
    CollectiveCommunication comm_;

    // the object which writes the convergence behaviour of the Newton
    // method to disk
    ConvergenceWriter convergenceWriter_;

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
