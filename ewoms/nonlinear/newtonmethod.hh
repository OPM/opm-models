// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/ErrorMacros.hpp>
#include <ewoms/common/propertysystem.hh>
#include <opm/material/common/ClassName.hpp>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/timer.hh>

#include <dune/istl/istlexception.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <iostream>
#include <sstream>

#include <unistd.h>

namespace Ewoms {
// forward declaration of classes
template <class TypeTag>
class NewtonMethod;
}

namespace Ewoms {
// forward declaration of property tags
namespace Properties {
//! The type tag on which the default properties for the Newton method
//! are attached
NEW_TYPE_TAG(NewtonMethod);

//! The simulation management class of the simulation
NEW_PROP_TAG(Simulator);

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

//! Specifies the type of a solution for a single degee of freedom
NEW_PROP_TAG(PrimaryVariables);

//! Vector containing a quantity of for equation on the whole grid
NEW_PROP_TAG(GlobalEqVector);

//! Vector containing a quantity of for equation for a single degee of freedom
NEW_PROP_TAG(EqVector);

//! The class which linearizes the non-linear system of equations
NEW_PROP_TAG(Linearizer);

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

/*!
 * \brief The value for the error below which convergence is declared
 *
 * This value can (and for the porous media models will) be changed to account for grid
 * scaling and other effects.
 */
NEW_PROP_TAG(NewtonRawTolerance);

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
SET_TYPE_PROP(NewtonMethod, NewtonConvergenceWriter, Ewoms::NullConvergenceWriter<TypeTag>);
SET_BOOL_PROP(NewtonMethod, NewtonWriteConvergence, false);
SET_BOOL_PROP(NewtonMethod, NewtonVerbose, true);
SET_SCALAR_PROP(NewtonMethod, NewtonRawTolerance, 1e-8);
// set the abortion tolerace to some very large value. if not
// overwritten at run-time this basically disables abortions
SET_SCALAR_PROP(NewtonMethod, NewtonMaxError, 1e100);
SET_INT_PROP(NewtonMethod, NewtonTargetIterations, 10);
SET_INT_PROP(NewtonMethod, NewtonMaxIterations, 18);
} // namespace Properties
} // namespace Ewoms

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
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) LinearSolverBackend;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonConvergenceWriter) ConvergenceWriter;

    typedef typename Dune::MPIHelper::MPICommunicator Communicator;
    typedef Dune::CollectiveCommunication<Communicator> CollectiveCommunication;

public:
    NewtonMethod(Simulator &simulator)
        : simulator_(simulator), endIterMsgStream_(std::ostringstream::out),
          linearSolver_(simulator), comm_(Dune::MPIHelper::getCommunicator()),
          convergenceWriter_(asImp_())
    {
        lastError_ = 1e100;
        error_ = 1e100;
        tolerance_ = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonRawTolerance);

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
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonRawTolerance,
                             "The maximum raw error tolerated by the Newton"
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
    { return error_ <= tolerance(); }

    /*!
     * \brief Returns a reference to the object describing the current physical problem.
     */
    Problem &problem()
    { return simulator_.problem(); }

    /*!
     * \brief Returns a reference to the object describing the current physical problem.
     */
    const Problem &problem() const
    { return simulator_.problem(); }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    Model &model()
    { return simulator_.model(); }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    const Model &model() const
    { return simulator_.model(); }

    /*!
     * \brief Returns the number of iterations done since the Newton method
     *        was invoked.
     */
    int numIterations() const
    { return numIterations_; }

    /*!
     * \brief Return the current tolerance at which the Newton method considers itself to
     *        be converged.
     */
    Scalar tolerance() const
    { return tolerance_; }

    /*!
     * \brief Set the current tolerance at which the Newton method considers itself to
     *        be converged.
     */
    void setTolerance(Scalar value)
    { tolerance_ = value; }

    /*!
     * \brief Run the Newton method.
     *
     * The actual implementation can influence all the strategic
     * decisions via callbacks using static polymorphism.
     */
    bool apply()
    {
        // Clear the current line using an ansi escape
        // sequence.  For an explanation see
        // http://en.wikipedia.org/wiki/ANSI_escape_code
        const char *clearRemainingLine = "\n";
        if (isatty(fileno(stdout))) {
            static const char blubb[] = { 0x1b, '[', 'K', '\r', 0 };
            clearRemainingLine = blubb;
        }

        SolutionVector &nextSolution = model().solution(/*historyIdx=*/0);
        SolutionVector currentSolution(nextSolution);
        GlobalEqVector solutionUpdate(nextSolution.size());

        Linearizer &linearizer = model().linearizer();

        // tell the implementation that we begin solving
        asImp_().begin_(nextSolution);

        linearizeTime_ = 0.0;
        solveTime_ = 0.0;
        updateTime_ = 0.0;

        Timer prePostProcessTimer;

        try {
            // execute the method as long as the implementation thinks
            // that we should do another iteration
            while (asImp_().proceed_()) {
                // linearize the problem at the current solution

                // notify the implementation that we're about to start
                // a new iteration
                prePostProcessTimer.start();
                asImp_().beginIteration_();
                prePostProcessTimer.stop();
                simulator_.addPrePostProcessTime(prePostProcessTimer.realTimeElapsed());


                // make the current solution to the old one
                currentSolution = nextSolution;

                if (asImp_().verbose_()) {
                    std::cout << "Linearize: r(x^k) = dS/dt + div F - q;   M = grad r"
                              << clearRemainingLine
                              << std::flush;
                }

                linearizeTimer_.start();
                asImp_().linearize_();
                linearizeTimer_.stop();
                linearizeTime_ += linearizeTimer_.realTimeElapsed();
                linearizeTimer_.halt();

                if (asImp_().verbose_()) {
                    std::cout << "Solve: M deltax^k = r"
                              << clearRemainingLine
                              << std::flush;
                }

                // solve the resulting linear equation system
                solveTimer_.start();

                // set the delta vector to zero before solving the linear system!
                solutionUpdate = 0;
                // ask the implementation to solve the linearized system
                auto b = linearizer.residual();
                bool converged = asImp_().solveLinear_(linearizer.matrix(), solutionUpdate, b);
                solveTimer_.stop();
                solveTime_ += solveTimer_.realTimeElapsed();;
                solveTimer_.halt();

                if (!converged) {
                    if (asImp_().verbose_())
                        std::cout << "Newton: Newton solver did not converge\n" << std::flush;
                    solveTimer_.stop();

                    prePostProcessTimer.start();
                    asImp_().failed_();
                    prePostProcessTimer.stop();
                    simulator_.addPrePostProcessTime(prePostProcessTimer.realTimeElapsed());

                    return false;
                }

                // update the solution
                if (asImp_().verbose_()) {
                    std::cout << "Update: x^(k+1) = x^k - deltax^k"
                              << clearRemainingLine
                              << std::flush;
                }

                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                updateTimer_.start();
                asImp_().updateError_(nextSolution,
                                      currentSolution,
                                      b,
                                      solutionUpdate);
                asImp_().update_(nextSolution, currentSolution, solutionUpdate, b);
                updateTimer_.stop();
                updateTime_ += updateTimer_.realTimeElapsed();
                updateTimer_.halt();

                // tell the implementation that we're done with this iteration
                prePostProcessTimer.start();
                asImp_().endIteration_(nextSolution, currentSolution);
                prePostProcessTimer.stop();
                simulator_.addPrePostProcessTime(prePostProcessTimer.realTimeElapsed());
            }
        }
        catch (const Dune::Exception &e)
        {
            linearizeTime_ += linearizeTimer_.realTimeElapsed();
            solveTime_ += solveTimer_.realTimeElapsed();
            updateTime_ += updateTimer_.realTimeElapsed();
            linearizeTimer_.halt();
            solveTimer_.halt();
            updateTimer_.halt();
            if (asImp_().verbose_())
                std::cout << "Newton method caught exception: \""
                          << e.what() << "\"\n" << std::flush;
            asImp_().failed_();
            return false;
        }
        catch (const Opm::NumericalIssue &e)
        {
            linearizeTime_ += linearizeTimer_.realTimeElapsed();
            solveTime_ += solveTimer_.realTimeElapsed();
            updateTime_ += updateTimer_.realTimeElapsed();
            linearizeTimer_.halt();
            solveTimer_.halt();
            updateTimer_.halt();

            if (asImp_().verbose_())
                std::cout << "Newton method caught exception: \""
                          << e.what() << "\"\n" << std::flush;
            asImp_().failed_();
            return false;
        }

        // clear current line on terminal
        if (asImp_().verbose_() && isatty(fileno(stdout)))
            std::cout << clearRemainingLine
                      << std::flush;

        // tell the implementation that we're done
        asImp_().end_();

        // print the timing summary of the time step
        if (asImp_().verbose_()) {
            Scalar elapsedTot =
                linearizeTime_
                + solveTime_
                + updateTime_;
            std::cout << "Linearization/solve/update time: "
                      << linearizeTime_ << "("
                      << 100 * linearizeTime_/elapsedTot << "%)/"
                      << solveTime_ << "("
                      << 100 * solveTime_/elapsedTot << "%)/"
                      << updateTime_ << "("
                      << 100 * updateTime_/elapsedTot << "%)"
                      << "\n" << std::flush;
        }


        // if we're not converged, tell the implementation that we've failed, else tell
        // it that we succeeded
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
    Scalar linearizeTime() const
    { return linearizeTime_; }

    /*!
     * \brief Returns the wall time spend so far for solving the
     *        linear systems for all iterations of the current time
     *        step.
     */
    Scalar solveTime() const
    { return solveTime_; }

    /*!
     * \brief Returns the wall time spend so far for updating the
     *        iterative solutions of the non-linear system for all
     *        iterations of the current time step.
     */
    Scalar updateTime() const
    { return updateTime_; }

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
        return EWOMS_GET_PARAM(TypeTag, bool, NewtonVerbose) && (comm_.rank() == 0);
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
    {
        problem().beginIteration();
        lastError_ = error_;
    }

    /*!
     * \brief Linearize the global non-linear system of equations.
     */
    void linearize_()
    { model().linearizer().linearize(); }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Opm::NumericalIssue if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    bool solveLinear_(const JacobianMatrix &A,
                      GlobalEqVector &x,
                      GlobalEqVector &b)
    { return linearSolver_.solve(A, x, b); }

    /*!
     * \brief Update the error of the solution given the previous
     *        iteration.
     *
     * For our purposes, the error of a solution is defined as the
     * maximum of the weighted residual of a given solution.
     *
     * \param nextSolution The solution after the current iteration
     * \param currentSolution The solution at the beginning the current iteration
     * \param currentResidual The residual (i.e., right-hand-side) of the current
     *                        iteration's solution.
     * \param solutionUpdate The difference between the current and the next solution
     */
    void updateError_(const SolutionVector &nextSolution,
                      const SolutionVector &currentSolution,
                      const GlobalEqVector &currentResidual,
                      const GlobalEqVector &solutionUpdate)
    {
        lastError_ = error_;

        // calculate the error as the maximum weighted tolerance of
        // the solution's residual
        error_ = 0;
        for (unsigned i = 0; i < currentResidual.size(); ++i) {
            if (i >= model().numGridDof() || model().dofTotalVolume(i) <= 0)
                continue;

            const auto &r = currentResidual[i];
            for (unsigned j = 0; j < r.size(); ++j)
                error_ = std::max(std::abs(r[j] * model().eqWeight(i, j)), error_);
        }

        // take the other processes into account
        error_ = comm_.max(error_);

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (error_ > EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxError))
            OPM_THROW(Opm::NumericalIssue,
                      "Newton: Error " << error_
                      << " is larger than maximum allowed error of "
                      << EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxError));
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * Different update strategies, such as chopped updates can be
     * implemented by overriding this method. The default behavior is
     * use the standard Newton-Raphson update strategy, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param nextSolution The solution vector after the current iteration
     * \param currentSolution The solution vector after the last iteration
     * \param solutionUpdate The delta vector as calculated by solving the linear system
     *                       of equations
     * \param currentResidual The residual vector of the current Newton-Raphson iteraton
     */
    void update_(SolutionVector &nextSolution,
                 const SolutionVector &currentSolution,
                 const GlobalEqVector &solutionUpdate,
                 const GlobalEqVector &currentResidual)
    {
        // first, write out the current solution to make convergence
        // analysis possible
        asImp_().writeConvergence_(currentSolution, solutionUpdate);

        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(solutionUpdate.two_norm2()))
            OPM_THROW(Opm::NumericalIssue, "Non-finite update!");

        for (unsigned dofIdx = 0; dofIdx < currentSolution.size(); ++dofIdx) {
            asImp_().updatePrimaryVariables_(dofIdx,
                                             nextSolution[dofIdx],
                                             currentSolution[dofIdx],
                                             solutionUpdate[dofIdx],
                                             currentResidual[dofIdx]);
        }
    }

    /*!
     * \brief Update a single primary variables object.
     */
    void updatePrimaryVariables_(int globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        nextValue = currentValue;
        nextValue -= update;
    }

    /*!
     * \brief Write the convergence behaviour of the newton method to
     *        disk.
     *
     * This method is called as part of the update proceedure.
     */
    void writeConvergence_(const SolutionVector &currentSolution,
                           const GlobalEqVector &solutionUpdate)
    {
        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence)) {
            convergenceWriter_.beginIteration();
            convergenceWriter_.writeFields(currentSolution, solutionUpdate);
            convergenceWriter_.endIteration();
        }
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param nextSolution The solution after the current Newton iteration
     * \param currentSolution The solution at the beginning of the current Newton iteration
     */
    void endIteration_(const SolutionVector &nextSolution,
                       const SolutionVector &currentSolution)
    {
        ++numIterations_;
        problem().endIteration();

        if (asImp_().verbose_()) {
            std::cout << "Newton iteration " << numIterations_ << ""
                      << " error: " << error_
                      << endIterMsg().str() << "\n" << std::flush;
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

    // optimal number of iterations we want to achieve
    int targetIterations_() const
    { return EWOMS_GET_PARAM(TypeTag, int, NewtonTargetIterations); }
    // maximum number of iterations we do before giving up
    int maxIterations_() const
    { return EWOMS_GET_PARAM(TypeTag, int, NewtonMaxIterations); }

    Simulator &simulator_;

    Ewoms::Timer linearizeTimer_;
    Ewoms::Timer solveTimer_;
    Ewoms::Timer updateTimer_;

    Scalar linearizeTime_;
    Scalar solveTime_;
    Scalar updateTime_;

    std::ostringstream endIterMsgStream_;

    Scalar error_;
    Scalar lastError_;
    Scalar tolerance_;

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
