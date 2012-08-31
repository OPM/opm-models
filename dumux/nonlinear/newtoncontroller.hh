// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2011 by Klaus Mosthaf                                     *
 *   Copyright (C) 2012 by Christoph Grueninger                              *
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
 * \brief Base controller class for the Newton solver.
 *
 * The actual discretizations are supposed to provide a more
 * specialized class which handles stuff like actually updating the
 * solution, etc.
 */
#ifndef DUMUX_NEWTON_CONTROLLER_HH
#define DUMUX_NEWTON_CONTROLLER_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/math.hh>

#include <dune/istl/istlexception.hh>
#include <dune/istl/ilu.hh> // required for Dune::MatrixBlockError. WTF?
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <sstream>
#include <string>

namespace Dumux {
namespace Properties {

//! create a new type tag on which the properties for the Newton method can be attached
NEW_TYPE_TAG(NewtonMethod);

//! The type of scalar values
NEW_PROP_TAG(Scalar);

//! Specifies the implementation of the Newton controller
NEW_PROP_TAG(NewtonController);

//! Specifies the type of the actual Newton method
NEW_PROP_TAG(NewtonMethod);

//! Specifies the type of a solution
NEW_PROP_TAG(SolutionVector);

//! Vector containing a quantity of for equation on the whole grid
NEW_PROP_TAG(GlobalEqVector);

//! Specifies the class of the physical problem
NEW_PROP_TAG(Problem);

//! Specifies the type of a global Jacobian matrix
NEW_PROP_TAG(JacobianMatrix);

//! Specifies the type of the linear solver to be used
NEW_PROP_TAG(LinearSolver);

//! Specifies the type of the class which writes out the Newton convergence
NEW_PROP_TAG(NewtonConvergenceWriter);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration (default is false)
NEW_PROP_TAG(NewtonWriteConvergence);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration (default is false)
NEW_PROP_TAG(ConvergenceWriter);

//! indicate whether the relative error should be used
NEW_PROP_TAG(NewtonEnableRelativeCriterion);

//! the value for the relative error below which convergence is declared
NEW_PROP_TAG(NewtonRelTolerance);

//! the maximum relative error which may occur in a simulation before
//! the newton method is aborted
NEW_PROP_TAG(NewtonMaxRelError);

//! indicate whether the absolute error should be used
NEW_PROP_TAG(NewtonEnableAbsoluteCriterion);

//! the value for the absolute error reduction below which convergence is declared
NEW_PROP_TAG(NewtonAbsTolerance);

//! indicate whether both of the criteria should be satisfied to declare convergence
NEW_PROP_TAG(NewtonSatisfyAbsAndRel);

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time-step size. The heuristic used
 * is to scale the last time-step size by the deviation of the
 * number of iterations used from the target steps.
 */
NEW_PROP_TAG(NewtonTargetSteps);

//! Number of maximum iterations for the Newton method.
NEW_PROP_TAG(NewtonMaxSteps);

// set default values
SET_BOOL_PROP(NewtonMethod, NewtonWriteConvergence, false);
SET_BOOL_PROP(NewtonMethod, NewtonEnableRelativeCriterion, true);
SET_BOOL_PROP(NewtonMethod, NewtonEnableAbsoluteCriterion, false);
SET_BOOL_PROP(NewtonMethod, NewtonSatisfyAbsAndRel, false);
SET_SCALAR_PROP(NewtonMethod, NewtonRelTolerance, 1e-8);
SET_SCALAR_PROP(NewtonMethod, NewtonAbsTolerance, 1e-5);
SET_SCALAR_PROP(NewtonMethod, NewtonMaxRelError, 1e100); // effectively disabled if not overwritten at run-time
SET_INT_PROP(NewtonMethod, NewtonTargetSteps, 10);
SET_INT_PROP(NewtonMethod, NewtonMaxSteps, 18);
}

/*!
 * \ingroup Newton
 * \brief A reference implementation of a Newton controller specific
 *        for the box scheme.
 *
 * If you want to specialize only some methods but are happy with the
 * defaults of the reference controller, derive your controller from
 * this class and simply overload the required methods.
 */
template <class TypeTag>
class NewtonController
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonConvergenceWriter) NewtonConvergenceWriterNewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonConvergenceWriter) NewtonConvergenceWriter;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;

    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;

public:
    NewtonController(Problem &problem)
        : endIterMsgStream_(std::ostringstream::out)
        , problem_(problem)
        , convergenceWriter_(asImp_())
        , linearSolver_(problem)
        , comm_(Dune::MPIHelper::getCommunicator())
    {
        enableRelativeCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableRelativeCriterion);
        enableAbsoluteCriterion_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, EnableAbsoluteCriterion);
        satisfyAbsAndRel_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, SatisfyAbsAndRel);
        if (!enableRelativeCriterion_ && !enableAbsoluteCriterion_)
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "At least one of Newton.EnableRelativeCriterion or "
                       "Newton.EnableAbsoluteCriterion has to be set to true");
        }
/*
        else if (satisfyAbsAndRel &&
                 (!enableRelativeCriterion_ || !enableAbsoluteCriterion_))
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "If you set Newton.SatisfyAbsAndRel to true, you also must set "
                       "Newton.EnableRelativeCriterion and Newton.EnableAbsoluteCriterion");

        }
*/

        setRelTolerance(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, RelTolerance));
        setAbsTolerance(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Newton, AbsTolerance));
        setTargetSteps(GET_PARAM_FROM_GROUP(TypeTag, int, Newton, TargetSteps));
        setMaxSteps(GET_PARAM_FROM_GROUP(TypeTag, int, Newton, MaxSteps));

        lastError_ = 1e100;
        lastAbsoluteError_ = 1e100;

        error_ = 1e100;
        absoluteError_ = 1e100;

        verbose_ = true;
        numSteps_ = 0;
    }

    /*!
     * \brief Returns a reference to the object which describes the
     *        physical problem.
     */
    const Problem &problem() const
    { return problem_; }
    Problem &problem()
    { return problem_; }

    /*!
     * \brief Set the maximum acceptable difference for convergence of
     *        any primary variable between two iterations.
     *
     * \param tolerance The maximum relative error between two Newton
     *                  iterations at which the scheme is considered
     *                  finished
     */
    void setRelTolerance(Scalar tolerance)
    { tolerance_ = tolerance; }

    /*!
     * \brief Set the maximum acceptable residual norm reduction.
     *
     * \param tolerance The maximum reduction of the residual norm
     *                  at which the scheme is considered finished
     */
    void setAbsTolerance(Scalar tolerance)
    { absoluteTolerance_ = tolerance; }

    /*!
     * \brief Set the number of iterations at which the Newton method
     *        should aim at.
     *
     * This is used to control the time-step size. The heuristic used
     * is to scale the last time-step size by the deviation of the
     * number of iterations used from the target steps.
     *
     * \param targetSteps Number of iterations which are considered "optimal"
     */
    void setTargetSteps(int targetSteps)
    { targetSteps_ = targetSteps; }

    /*!
     * \brief Set the number of iterations after which the Newton
     *        method gives up.
     *
     * \param maxSteps Number of iterations after we give up
     */
    void setMaxSteps(int maxSteps)
    { maxSteps_ = maxSteps; }

    /*!
     * \brief Returns true if another iteration should be done.
     *
     * \param uCurrentIter The solution of the current Newton iteration
     */
    bool newtonProceed(const SolutionVector &uCurrentIter)
    {
        if (numSteps_ < 2)
            return true; // we always do at least two iterations
        else if (asImp_().newtonConverged()) {
            return false; // we are below the desired tolerance
        }
        else if (numSteps_ >= maxSteps_) {
            // we have exceeded the allowed number of steps.  if the
            // relative error was reduced by a factor of at least 4,
            // we proceed even if we are above the maximum number of
            // steps
            if (enableRelativeCriterion_)
                return error_*4.0 < lastError_;
            else
                return absoluteError_*4.0 < lastAbsoluteError_;
        }

        return true;
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool newtonConverged() const
    {
        if (enableRelativeCriterion_ && !enableAbsoluteCriterion_)
            // only look at the relative error
            return error_ <= tolerance_;
        else if (!enableRelativeCriterion_ && enableAbsoluteCriterion_)
            // only look at the absolute error
            return absoluteError_ <= absoluteTolerance_;
        else if (satisfyAbsAndRel_)
            // both, the absolute and the relative tolerances must be attained
            return 
                error_ <= tolerance_
                && absoluteError_ <= absoluteTolerance_;
        else
            // we're done as soon as either the absolute or the
            // relative tolerance is achieved.
            return error_ <= tolerance_
                || absoluteError_ <= absoluteTolerance_;

        return false;
    }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param method The object where the NewtonMethod is executed
     * \param u The initial solution
     */
    void newtonBegin(NewtonMethod &method, const SolutionVector &u)
    {
        method_ = &method;
        numSteps_ = 0;

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence))
            convergenceWriter_.beginTimestep();
    }

    /*!
     * \brief Assemble the global linear system of equations.
     */
    void newtonAssemble()
    { problem().model().jacobianAssembler().assemble(); }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void newtonBeginStep()
    {
        lastError_ = error_;
        lastAbsoluteError_ = absoluteError_;
    }

    /*!
     * \brief Returns the number of steps done since newtonBegin() was
     *        called.
     */
    int newtonNumSteps()
    { return numSteps_; }

    /*!
     * \brief Update the relative error of the solution compared to
     *        the previous iteration.
     *
     * The relative error can be seen as a norm of the difference
     * between the current and the next iteration.
     *
     * \param uLastIter The current iterative solution
     * \param deltaU The difference between the current and the next solution
     */
    void newtonUpdateRelError(const SolutionVector &uCurrentIter,
                              const SolutionVector &uLastIter,
                              const GlobalEqVector &deltaU)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The Newton controller (" << Dune::className<Implementation>() << ") "
                   "does not provide a newtonUpdateRelError() method!");
    }

    /*!
     * \brief Update the absolute error of the solution compared to
     *        the previous iteration.
     */
    void newtonUpdateAbsError(const SolutionVector &uCurrentIter,
                              const SolutionVector &uLastIter,
                              const GlobalEqVector &deltaU)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "The Newton controller (" << Dune::className<Implementation>() << ") "
                   "does not provide a newtonUpdateAbsError() method!");
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     *
     * \param A The matrix of the linear system of equations
     * \param x The vector which solves the linear system
     * \param b The right hand side of the linear system
     */
    void newtonSolveLinear(const JacobianMatrix &A,
                           GlobalEqVector &x,
                           const GlobalEqVector &b)
    {
        try {
            if (numSteps_ == 0)
            {
                Scalar norm2 = b.two_norm2();
                norm2 = comm_.sum(norm2);
            }
            
            int converged = linearSolver_.solve(A, x, b);

            // make sure all processes converged
            int convergedRemote = converged;
            convergedRemote = comm_.min(converged);

            if (!converged) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge");
            }
            else if (!convergedRemote) {
                DUNE_THROW(NumericalProblem,
                           "Linear solver did not converge on a remote process");
            }
        }
        catch (Dune::MatrixBlockError e) {
            // make sure all processes converged
            int converged = 0;
            converged = comm_.min(converged);

            Dumux::NumericalProblem p;
            std::string msg;
            std::ostringstream ms(msg);
            ms << e.what() << "M=" << A[e.r][e.c];
            p.message(ms.str());
            throw p;
        }
        catch (const Dune::Exception &e) {
            // make sure all processes converged
            int converged = 0;
            converged = comm_.min(converged);

            Dumux::NumericalProblem p;
            p.message(e.what());
            throw p;
        }
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * The error estimates required for the newtonConverged() and
     * newtonProceed() methods should be updated inside this method.
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
    void newtonUpdateErrors(const SolutionVector &uCurrentIter,
                            const SolutionVector &uLastIter,
                            const GlobalEqVector &deltaU)
    {
        lastError_ = error_;
        lastAbsoluteError_ = absoluteError_;
        asImp_().newtonUpdateRelError(uCurrentIter, uLastIter, deltaU);
        asImp_().newtonUpdateAbsError(uCurrentIter, uLastIter, deltaU);       
        writeConvergence_(uLastIter, deltaU);
    }

    /*!
     * \brief Indicates that one Newton iteration was finished.
     *
     * \param uCurrentIter The solution after the current Newton iteration
     * \param uLastIter The solution at the beginning of the current Newton iteration
     */
    void newtonEndStep(const SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ++numSteps_;

        if (verbose())
        {
            std::cout << "\rNewton iteration " << numSteps_ << " done";
            if (enableRelativeCriterion_)
                std::cout << ", relative error = " << error_;
            if (enableAbsoluteCriterion_)
                std::cout << ", absolute error = " << absoluteError_;
            std::cout << endIterMsg().str() << "\n";
        }

        endIterMsgStream_.str("");
    }

    /*!
     * \brief Indicates that we're done solving the non-linear system
     *        of equations.
     */
    void newtonEnd()
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence))
            convergenceWriter_.endTimestep();
    }

    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonFail()
    { numSteps_ = targetSteps_*2; }

    /*!
     * \brief Called when the Newton method was successful.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed()
    { }

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
        // that we want to avoid failing in the next Newton
        // iteration which would require another linearization
        // of the problem.
        if (numSteps_ > targetSteps_) {
            Scalar percent = Scalar(numSteps_ - targetSteps_)/targetSteps_;
            return oldTimeStep/(1.0 + percent);
        }
        
        Scalar percent = Scalar(targetSteps_ - numSteps_)/targetSteps_;
        return oldTimeStep*(1.0 + percent/1.2);
    }

    std::ostringstream &endIterMsg()
    { return endIterMsgStream_; }

    /*!
     * \brief Specifies if the Newton method ought to be chatty.
     */
    void setVerbose(bool val)
    { verbose_ = val; }

    /*!
     * \brief Returns true if the Newton method ought to be chatty.
     */
    bool verbose() const
    { return verbose_ && comm_.rank() == 0; }

protected:
    std::ostringstream endIterMsgStream_;

    NewtonMethod *method_;
    Problem &problem_;

    NewtonConvergenceWriterNewtonMethod convergenceWriter_;

    // relative errors and tolerance
    Scalar error_;
    Scalar lastError_;
    Scalar tolerance_;

    // absolute errors and tolerance
    Scalar absoluteError_;
    Scalar lastAbsoluteError_;
    Scalar absoluteTolerance_;
    
    // which criteria do we have to satisfy in order to be converged?
    bool enableRelativeCriterion_;
    bool enableAbsoluteCriterion_;
    bool satisfyAbsAndRel_;

    // optimal number of iterations we want to achieve
    int targetSteps_;
    // maximum number of iterations we do before giving up
    int maxSteps_;
    // actual number of steps done so far
    int numSteps_;

    // are we supposed to be chatty or not?
    bool verbose_;

    // the linear solver
    LinearSolver linearSolver_;

    // the collective communication used by the simulation (i.e. fake
    // or MPI)
    CollectiveCommunication comm_;

private:
    void writeConvergence_(const SolutionVector &uLastIter,
                           const GlobalEqVector &deltaU)
    {
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Newton, WriteConvergence)) {
            convergenceWriter_.beginIteration();
            convergenceWriter_.writeFields(uLastIter, deltaU);
            convergenceWriter_.endIteration();
        }
    }

    // returns the actual implementation for the controller we do
    // it this way in order to allow "poor man's virtual methods",
    // i.e. methods of subclasses which can be called by the base
    // class.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

};
} // namespace Dumux

#endif
