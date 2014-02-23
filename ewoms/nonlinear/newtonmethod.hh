/*
  Copyright (C) 2008-2013 by Andreas Lauser
  Copyright (C) 2009-2011 by Bernd Flemisch

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

//! The value for the error below which convergence is declared
NEW_PROP_TAG(NewtonTolerance);

//! The maximum error which may occur in a simulation before the
//! Newton method for the time step is aborted
NEW_PROP_TAG(NewtonMaxError);

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
SET_TYPE_PROP(NewtonMethod, NewtonConvergenceWriter,
              Ewoms::NullConvergenceWriter<TypeTag>);
SET_BOOL_PROP(NewtonMethod, NewtonWriteConvergence, false);
SET_BOOL_PROP(NewtonMethod, NewtonVerbose, true);
SET_SCALAR_PROP(NewtonMethod, NewtonTolerance, 1e-8);
// set the abortion tolerace to some very large value. if not
// overwritten at run-time this basically disables abortions
SET_SCALAR_PROP(NewtonMethod, NewtonMaxError, 1e100);
SET_INT_PROP(NewtonMethod, NewtonTargetIterations, 10);
SET_INT_PROP(NewtonMethod, NewtonMaxIterations, 18);
} // namespace Properties
} // namespace Opm

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
        : problem_(problem), endIterMsgStream_(std::ostringstream::out),
          linearSolver_(problem), comm_(Dune::MPIHelper::getCommunicator()),
          convergenceWriter_(asImp_())
    {
        lastError_ = 1e100;
        error_ = 1e100;

        numIterations_ = 0;
    }

    /*!
     * \brief Register all run-time parameters for the Newton method.
     */
    static void registerParameters()
    {
        LinearSolverBackend::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonVerbose,
                             "Specify whether the Newton method should inform "
                             "the user about its progress or not");
        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonWriteConvergence,
                             "Write the convergence behaviour of the Newton "
                             "method to a VTK file");
        EWOMS_REGISTER_PARAM(TypeTag, int, NewtonTargetIterations,
                             "The 'optimimum' number of Newton iterations per "
                             "time step");
        EWOMS_REGISTER_PARAM(TypeTag, int, NewtonMaxIterations,
                             "The maximum number of Newton iterations per time "
                             "step");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonTolerance,
                             "The maximum error tolerated by the Newton "
                             "method for considering a solution to be "
                             "converged");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonMaxError,
                             "The maximum error tolerated by the Newton "
                             "method to which does not cause an abort");
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool converged() const
    { return error_ <= tolerance_(); }

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
     * The actual implementation can influence all the strategic
     * decisions via callbacks using static polymorphism.
     */
    bool apply()
    {
        SolutionVector &currentSolution = model().solution(/*historyIdx=*/0);
        SolutionVector previousSolution(currentSolution);
        GlobalEqVector solutionUpdate(currentSolution.size());

        JacobianAssembler &jacobianAsm = model().jacobianAssembler();

        // tell the implementation that we begin solving
        asImp_().begin_(currentSolution);

        assembleTimer_.reset();
        solveTimer_.reset();
        updateTimer_.reset();

        // execute the method as long as the implementation thinks
        // that we should do another iteration
        while (asImp_().proceed_()) {
            // notify the implementation that we're about to start
            // a new iteration
            asImp_().beginIteration_();

            // make the current solution to the old one
            previousSolution = currentSolution;

            if (asImp_().verbose_()) {
                std::cout << "Assemble: r(x^k) = dS/dt + div F - q;   M = grad "
                             "r";
                std::cout.flush();
            }

            // linearize the problem at the current solution
            try
            {
                assembleTimer_.start();
                asImp_().linearize_();
                assembleTimer_.stop();
            }
            catch (const Dune::Exception &e)
            {
                if (asImp_().verbose_())
                    std::cout << "Newton: Caught Dune exception during linearization: \""
                              << e.what() << "\"\n";
                asImp_().failed_();
                return false;
            }
            catch (const Opm::NumericalProblem &e)
            {
                if (asImp_().verbose_())
                    std::cout << "Newton: Caught Opm exception during linearization: \""
                              << e.what() << "\"\n";
                asImp_().failed_();
                return false;
            };

            // Clear the current line using an ansi escape
            // sequence.  For an explanation see
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
            solutionUpdate = 0;
            // ask the implementation to solve the linearized system
            if (!asImp_().solveLinear_(jacobianAsm.matrix(), solutionUpdate,
                                       jacobianAsm.residual())) {
                if (asImp_().verbose_())
                    std::cout << "Newton: Linear solver did not converge\n";
                solveTimer_.stop();
                asImp_().failed_();
                return false;
            }
            solveTimer_.stop();

            // update the solution
            if (asImp_().verbose_()) {
                std::cout << "\rUpdate: x^(k+1) = x^k - deltax^k";
                std::cout << clearRemainingLine;
                std::cout.flush();
            }

            try
            {
                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                updateTimer_.start();
                asImp_().updateError_(currentSolution, previousSolution,
                                      jacobianAsm.residual(), solutionUpdate);
                asImp_().update_(currentSolution, previousSolution, solutionUpdate);
                updateTimer_.stop();

                // tell the implementation that we're done with this
                // iteration
                asImp_().endIteration_(currentSolution, previousSolution);
            }
            catch (const Dune::Exception &e)
            {
                updateTimer_.stop();
                if (asImp_().verbose_())
                    std::cout << "Newton: Caught exception during update: \""
                              << e.what() << "\"\n";
                asImp_().failed_();
                return false;
            }
            catch (const Opm::NumericalProblem &e)
            {
                if (asImp_().verbose_())
                    std::cout << "Newton: Caught Opm exception during linearization: \""
                              << e.what() << "\"\n";
                asImp_().failed_();
                return false;
            };
        }

        // tell the implementation that we're done
        asImp_().end_();

        if (asImp_().verbose_()) {
            Scalar elapsedTot =
                assembleTimer_.elapsed() + solveTimer_.elapsed() + updateTimer_.elapsed();
            std::cout << "Assemble/solve/update time: "
                      << assembleTimer_.elapsed() << "("
                      << 100 * assembleTimer_.elapsed() / elapsedTot << "%)/"
                      << solveTimer_.elapsed() << "("
                      << 100 * solveTimer_.elapsed() / elapsedTot << "%)/"
                      << updateTimer_.elapsed() << "("
                      << 100 * updateTimer_.elapsed() / elapsedTot << "%)"
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
    { return updateTimer_.elapsed(); }

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
            Scalar percent = Scalar(numIterations_ - targetIterations_())
                             / targetIterations_();
            return oldTimeStep / (1.0 + percent);
        }

        Scalar percent = Scalar(targetIterations_() - numIterations_)
                         / targetIterations_();
        return oldTimeStep * (1.0 + percent / 1.2);
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
    {
        return EWOMS_GET_PARAM(TypeTag, bool, NewtonVerbose) && comm_.rank()
                                                                == 0;
    }

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
    { lastError_ = error_; }

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
    bool solveLinear_(const JacobianMatrix &A, GlobalEqVector &x,
                      const GlobalEqVector &b)
    { return linearSolver_.solve(A, x, b); }

    /*!
     * \brief Update the error of the solution given the previous
     *        iteration.
     *
     * For our purposes, the error of a solution is defined as the
     * maximum of the weighted residual of a given solution.
     *
     * \param currentSolution The current iterative solution
     * \param previousSolution The last iterative solution
     * \param previousResidual The residual (i.e., right-hand-side) of
     *                         the previous iteration's solution.
     * \param solutionUpdate The difference between the current and the next solution
     */
    void updateError_(const SolutionVector &currentSolution,
                      const SolutionVector &previousSolution,
                      const GlobalEqVector &previousResidual,
                      const GlobalEqVector &solutionUpdate)
    {
        lastError_ = error_;

        // calculate the error as the maximum weighted tolerance of
        // the solution's residual
        error_ = 0;
        for (unsigned i = 0; i < previousResidual.size(); ++i) {
            if (model().dofTotalVolume(i) <= 0)
                continue;

            const auto &r = previousResidual[i];
            for (unsigned j = 0; j < r.size(); ++j)
                error_ = std::max(std::abs(r[j] * model().eqWeight(i, j)), error_);
        }

        // take the other processes into account
        error_ = comm_.max(error_);

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (error_ > EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxError))
            OPM_THROW(Opm::NumericalProblem,
                      "Newton: Error "
                      << error_
                      << " is larger than maximum allowed error of "
                      << EWOMS_GET_PARAM(TypeTag, Scalar,
                                         NewtonMaxError));
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * Different update strategies, such as chopped updates can be
     * implemented by overriding this method. The default behavior is
     * use the standard Newton-Raphson update strategy, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param currentSolution The solution vector after the current iteration
     * \param previousSolution The solution vector after the last iteration
     * \param solutionUpdate The delta vector as calculated by solving
     *                       the linear system of equations
     */
    void update_(SolutionVector &currentSolution,
                 const SolutionVector &previousSolution,
                 const GlobalEqVector &solutionUpdate)
    {
        // first, write out the current solution to make convergence
        // analysis possible
        writeConvergence_(previousSolution, solutionUpdate);

        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(solutionUpdate.two_norm2()))
            OPM_THROW(Opm::NumericalProblem, "Non-finite update!");

        for (unsigned i = 0; i < previousSolution.size(); ++i) {
            currentSolution[i] = previousSolution[i];
            currentSolution[i] -= solutionUpdate[i];
        }
    }

    /*!
     * \brief Write the convergence behaviour of the newton method to
     *        disk.
     *
     * This method is called as part of the update proceedure.
     */
    void writeConvergence_(const SolutionVector &previousSolution,
                           const GlobalEqVector &solutionUpdate)
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence)) {
            convergenceWriter_.beginIteration();
            convergenceWriter_.writeFields(previousSolution, solutionUpdate);
            convergenceWriter_.endIteration();
        }
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param currentSolution The solution after the current Newton iteration
     * \param previousSolution The solution at the beginning of the current Newton iteration
     */
    void endIteration_(const SolutionVector &currentSolution,
                       const SolutionVector &previousSolution)
    {
        ++numIterations_;

        if (asImp_().verbose_()) {
            std::cout << "\rNewton iteration " << numIterations_ << "";
            std::cout << " error: " << error_;
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
            return error_ * 4.0 < lastError_;
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
    { numIterations_ = targetIterations_() * 2; }

    /*!
     * \brief Called if the Newton method was successful.
     *
     * This method is called _after_ end_()
     */
    void succeeded_()
    {}

    Problem &problem_;

    Dune::Timer assembleTimer_;
    Dune::Timer solveTimer_;
    Dune::Timer updateTimer_;

    std::ostringstream endIterMsgStream_;

    // errors and tolerance
    Scalar error_;
    Scalar lastError_;
    Scalar tolerance_() const
    { return EWOMS_GET_PARAM(TypeTag, Scalar, NewtonTolerance); }

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
    { return *static_cast<Implementation *>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Ewoms

#endif
