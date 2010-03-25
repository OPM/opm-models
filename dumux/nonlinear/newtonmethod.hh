// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
#ifndef DUMUX_NEWTON_METHOD_HH
#define DUMUX_NEWTON_METHOD_HH

#include "newtonconvergencewriter.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/preconditionerpdelab.hh>

// HACK!!!! this is actually not specifc to DuMuX and should be put
// into PDELab after some cleanup!
#include <dumux/common/preconditionerpdelab.hh>

namespace Dumux
{

namespace Properties
{
//! specifies the verbosity of the linear solver (by default it is 0,
//! i.e. it doesn't print anything)
NEW_PROP_TAG(NewtonLinearSolverVerbosity);

//! specifies whether the convergence rate and the global defect gets
//! written out to disk for every newton iteration (default is false)
NEW_PROP_TAG(NewtonWriteConvergence);

//! specifies whether the update should be done using the line search
//! method instead of the "raw" newton method. Whether this property
//! has any effect depends on wether the line search method is
//! implemented for the actual model's newton method. By default we do
//! not use line search.
NEW_PROP_TAG(NewtonUseLineSearch);

SET_PROP_DEFAULT(NewtonLinearSolverVerbosity)
{ static const int value = 0; };

SET_PROP_DEFAULT(NewtonWriteConvergence)
{ static const bool value = false; };

SET_PROP_DEFAULT(NewtonUseLineSearch)
{ static const bool value = false; };
};

/*!
 * \brief The algorithmic part of the multi dimensional newton method.
 *
 * In order to use the method you need a NewtonController.
 */
template <class TypeTag>
class NewtonMethod
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridSolution)) GridSolution;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridFunctionSpace)) GridFunctionSpace;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridJacobian)) GridJacobian;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridDefect))   GridDefect;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVarVector)) PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridBlockMatrixContainer)) GridBlockMatrixContainer;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridBlockVectorContainer)) GridBlockVectorContainer;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridBlockVector)) GridBlockVector;
    
    enum { writeConvergence = GET_PROP_VALUE(TypeTag, PTAG(NewtonWriteConvergence)) };
    enum { useLineSearch = GET_PROP_VALUE(TypeTag, PTAG(NewtonUseLineSearch)) };
    enum { dim = GridView::dimensionworld };

    typedef Dune::NewtonConvergenceWriter<TypeTag, writeConvergence> ConvergenceWriter;


public:
    NewtonMethod(const Box::GlobalContext<TypeTag> &ctx)
        : endIterMsgStream_(std::ostringstream::out),
          update_(ctx.gridFunctionSpace())
    { 
        setTolerance(1e-7);
        setTargetIterations(8);
        setMaxIterations(15);
    }

    ~NewtonMethod()
    { }

    /*!
     * \brief Set the acceptable change for a primary variable between
     *        two iterations.
     */
    void setTolerance(Scalar tol)
    { tolerance_ = tol; };

    /*!
     * \brief Returns the acceptable change for a primary variable
     *        between two iterations.
     */
    Scalar tolerance() const
    { return tolerance_; }

    /*!
     * \brief Set the number of iterations which are considered
     *        optimal by the newton method.
     */
    void setTargetIterations(int value)
    { targetIterations_ = value; }

    /*!
     * \brief Returns the number of iterations which are considered
     *        optimal by the newton method.
     */
    int targetIterations()
    { return targetIterations_; }

    /*!
     * \brief Set the maximum number of iterations before the newton
     *        method is aborted.
     */
    void setMaxIterations(int value)
    { maxIterations_ = value; }

    /*!
     * \brief Returns the maximum number of iterations before the
     *        newton method is aborted.
     */
    void maxIterations(int value)
    { return maxIterations_; }
    
	/*!
     * \brief Returns the number of iterations which where done since
     *        the invokation of the newton method.
     */
	int numIterations() const
	{ return numIterations_; }

    /*!
     * \brief Returns true iff the newton method should be verbose
     *        about what it is going on;
     */
    template <class GlobalContext>
    bool verbose(const GlobalContext &context) const
    { return context.gridView().comm().rank() == 0; }

    /*!
     * \brief Returns the ostream for the message which is printed at
     *        the end of each iteration.
     *
     * eoi means "End Of Iteration"
     */
	std::ostringstream &eoiMsg() const
	{ return endIterMsgStream_; }

    /*!
     * \brief Run the newton method.
     */
    template <class GlobalContext>
    bool nonLinearSolve(GlobalContext &context)
    {
        try {
            convergenceWriter_.beginNonLinearSolve();
            bool result = nonLinearSolve_(context);
            convergenceWriter_.endNonLinearSolve();
            return result;
        }
        catch (const Dumux::NumericalProblem &e) {
            if (verbose(context))
                std::cout << "Newton: Caught non-critical exception: \"" << e.what() << "\"\n";
            asImp_().newtonEnd_(context);
            asImp_().newtonFailed_(context);
            // tell the convergence writer that we're finished
            convergenceWriter_.endNonLinearSolve(); 
            return false;
        };
    }

	/*!
     * \brief Suggest a size size the next time step [sec] based on
     *        the current time step size.
     *
     * The default behaviour is to suggest the old time step size
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time step.
     */
    template <class GlobalContext>
	Scalar suggestTimeStepSize(GlobalContext &context) const
	{
        Scalar dt = context.timeManager().timeStepSize();

		// be agressive when reducing the time step size but
		// conservative when increasing it. the rationale is that we
		// want to avoid failing for the next invocation of the newton
		// method which would be probably quite expensive.
		if (numIterations_ > targetIterations_) {
			Scalar percent = 
                (Scalar(numIterations_) - targetIterations_)/targetIterations_;
			return dt/(1.0 + percent);
		}
		else {
			Scalar percent = 
                (Scalar(targetIterations_) - numIterations_)/targetIterations_;
			return dt*(1.0 + percent/1.2);
		}
	}


protected:
	/*!
	* \brief Returns true iff another iteration should be done.
	*/
    template <class GlobalContext>
	bool proceed_(const GlobalContext &context) const
	{
		if (numIterations_ < 2)
			return true; // we always do at least two iterations
		else if (numIterations_ <= maxIterations_) 
            // we have not yet reached the maximum number of
            // iterations, so we do another iteration if we have not
            // yet converged!
            return !asImp_().isConverged_(context);

        // we have exceeded the allowed number of steps.  if the
        // relative error was reduced by a factor of at least 4, we
        // proceed even in this case
        if (error_*4.0 < lastError_)
            return !asImp_().isConverged_(context);

        // we have exceeded the maximum number of iterations
        return false;
	}

	/*!
     * \brief Returns true iff the error of the current iterative
     *        solution is below the desired tolerance.
     */
    template <class GlobalContext>
	bool isConverged_(const GlobalContext &context) const
	{ return error_ <= tolerance_; }

    /*!
     * \brief Called before the newton method is applied to an
     *        non-linear system of equations.
     */
    template <class GlobalContext>
	void newtonBegin_(const GlobalContext &context)
	{
		numIterations_ = 0;
        error_ = 1e100;
        lastError_ = 0.0;
    }

	/*!
     * \brief Indidicates the beginning of a newton iteration.
     */
    template <class GlobalContext>
	void beginIteration_(const GlobalContext &context)
	{ lastError_ = error_; }

	/*!
     * \brief Update the error of the solution compared to the
     *        previous iteration.
     */
    template <class GlobalContext>
	void updateRelError_(const GlobalContext &context,
                         const GridBlockVector &deltaU)
	{
        const GridSolution &sol = context.nextGridSolution();

		// calculate the relative error as the maximum relative
		// deflection in any degree of freedom.
		error_ = 0;
		for (int i = 0; i < deltaU.size(); ++i) {
			for (int j = 0; j < PrimaryVarVector::size; ++j) {
				Scalar tmp =
                    std::abs((*deltaU)[i][j])
                    / 
                    std::max(std::abs(sol[i][j]), Scalar(1e-4));
				error_ = std::max(error_, tmp);
			}
		};

		error_ = context.gridView().comm().max(error_);
	}

    template <class GlobalContext>
    void addUpdate_(GlobalContext &context, 
                    GridBlockVector &update)
    {
        if (!useLineSearch) {
            // let the grid defect do the update of the next solution,
            // so that the global defect can be calculated in one go
            context.gridDefect().addUpdate(context, update, error_);
            return;
        }
        
        // save the solution of the current iteration
        GridSolution tmpSol(context.nextGridSol());

        Scalar oldDefect = context.gridDefect().norm();
        const int maxDivisions = 3;
        int nDiv;
        for (nDiv = 0; nDiv <= maxDivisions; ++nDiv) {
            // let the grid defect do the update of the next solution,
            // so that the global defect can be calculated in one go
            if (nDiv == 0)
                context.gridDefect().addUpdate(context, update, error_);
            else {
                Scalar dummy;
                context.gridDefect().addUpdate(context, update, dummy);
            }
            
            Scalar newDefect = context.gridDefect().norm();
            if (oldDefect > newDefect)
                break;

            // the old defect was smaller than the one of the current
            // iterative solution. We restore the old solution and try
            // with half the update length.
            context.setNextGridSol(tmpSol);
            update /= 2;
        }
        
        eoiMsg() << ", defect="<<context.gridDefect().norm()<<", divisions="<<nDiv;
    };

	/*!
     * \brief Indidicates the end of a newton iteration.
     */
    template <class GlobalContext>
	void endIteration_(const GlobalContext &context)
	{
        ++numIterations_;
        lastError_ = error_;

        if (verbose(context))
            std::cout << "\rNewton iteration " << numIterations_ << " done: "
                      << "error=" << error_ << eoiMsg().str() << "\n";
        endIterMsgStream_.str("");
    }

    /*!
     * \brief Called after the newton method was applied to an
     *        non-linear system of equations.
     */
    template <class GlobalContext>
	void newtonEnd_(const GlobalContext &context)
	{ }

    /*!
     * \brief Called after newtonEnd() if the newton method
     *        converged.
     */
    template <class GlobalContext>
	void newtonSucceeded_(const GlobalContext &context)
	{ }

    /*!
     * \brief Called after newtonEnd() if the newton method
     *        did not converge.
     */
    template <class GlobalContext>
	void newtonFailed_(const GlobalContext &context)
	{ }

    /*! 
     * \brief The actual newton method.
     */
    template <class GlobalContext>
    bool nonLinearSolve_(GlobalContext &context)
    {
        // tell the implementation that we begin solving the
        // non-linear system of equations
        asImp_().newtonBegin_(context);
       
        // execute the method as long as the controller thinks
        // that we should do another iteration
        while (asImp_().proceed_(context))
        {
            // notify the controller that we're about to start
            // a new timestep
            asImp_().beginIteration_(context);

            if (verbose(context)) {
                std::cout << "Assembling global jacobian";
                std::cout.flush();
            }

            // linearize the problem at the current solution
            context.gridJacobian().update(context);
            context.gridDefect().update(context, 
                                        true /* applyDirichlet? yes! */);
            
            // solve the resulting linear equation system
            if (verbose(context)) {
                std::cout << "\rSolve Mx = r               ";
                std::cout.flush();
            }

            std::cout << "HALLO:" << context.gridDefect().blockVector() << "\n";
            
            // let the implementation solve the linear system of
            // equations
            asImp_().linearSolve_(context,
                                  context.gridJacobian().blockMatrix(),
                                  update_,
                                  context.gridDefect().blockVector());

            // we want to add the delta vector to the current
            // solution, so we multiply everything with -1
            update_ *= -1.0;

            // make the implementation add the update vector. This
            // method also updates the global defect.
            asImp_().addUpdate_(context, update_);
            
            // write the convergence behaviour into a VTK file (only
            // if enabled)
            convergenceWriter_.writeFields(context, update_);

            // tell the controller that we're done with this iteration
            asImp_().endIteration_(context);
        }


        // tell the implementation that we're done
        asImp_().newtonEnd_(context);
        
        if (!asImp_().isConverged_(context)) {
            asImp_().newtonFailed_(context);
            return false;
        }
        
        asImp_().newtonSucceeded_(context);
        return true;
    };

	/*!
     * \brief Solve a linear system of equations \f$ \mathbf{A}x - b =
     *        0\f$.
     *
     * Throws Dumux::NumericalProblem if the linear solver didn't
     * converge.
     */
	void linearSolve_(Box::GlobalContext<TypeTag> &context,
                      GridBlockMatrixContainer&A,
                      GridBlockVectorContainer &x,
                      GridBlockVectorContainer &b)
	{
		// if the deflection of the newton method is large, we do not
		// need to solve the linear approximation accurately. Assuming
		// that the initial value for the delta vector u is quite
		// close to the final value, a reduction of 6 orders of
		// magnitute in the defect should be sufficient...
		Scalar residReduction = 1e-6;

        int verbosity = 0;
        if (verbose(context))
            verbosity =
                GET_PROP_VALUE(TypeTag,
                               PTAG(NewtonLinearSolverVerbosity));
        
        /*
        // Parallel Loop solver with Pardiso as sequential solver
        typedef Dune::PDELab::ISTLBackend_NoOverlap_Loop_Pardiso<TypeTag> Solver;
        Solver solver(context.gridFunctionSpace(),
                      context.constraintsTrafo(), 5000, verbosity);
        */
        
        /*
        // Parallel AMG with no overlap and SSOR as preconditioner
        typedef Dune::PDELab::ISTLBackend_NoOverlap_AMG_SSOR<GridFunctionSpace> Solver;
        Solver solver(context.gridFunctionSpace(), verbosity);
        */
        
        /*
        // Parallel stabilized BiCG with no overlap without
        // preconditioner
        typedef Dune::PDELab::ISTLBackend_NoOverlap_BCGS_NOPREC<GridFunctionSpace> Solver;
        Solver solver(context.gridFunctionSpace(), 5000, verbosity);
        */

        // Parallel stabilized BiCG with no overlap and ILU as
        // preconditioner
        typedef Dumux::PDELab::ISTLBackend_NoOverlap_BCGS_ILU<TypeTag> Solver;
        Solver solver(context, 5000, verbosity);

		try {
            // set the initial guess for x to the zero vector
            x = 0.0;
            // apply the solver
            solver.apply(A, x, b, residReduction);
		}
		catch (Dune::MatrixBlockError e) {
            // the solver could not invert a matrix block, that's
            // usually because it is singular (or almost singular).
            // We throw a NumericalProblem and try again with a
            // smaller time step size if the models wants this...
            Dumux::NumericalProblem p;
			p.message(e.what());
			throw p;
		}
        
        if (!solver.result().converged)
            DUNE_THROW(Dumux::NumericalProblem,
                       "Solving the linear system of equations did not converge.");
    }
    
private:
    Implementation &asImp_() 
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
    
    Scalar   lastError_; // relative error of the last iteration
    Scalar   error_; // relative error of the current iteration
    Scalar   tolerance_; // the maximal relative error which is tolerated

    int numIterations_;
    int targetIterations_;
    int maxIterations_;

    mutable std::ostringstream endIterMsgStream_;
    
    GridBlockVectorContainer update_;
    ConvergenceWriter convergenceWriter_;
};

} // namespace Dumux

#endif
