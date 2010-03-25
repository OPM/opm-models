// $Id$
/*****************************************************************************
*   Copyright (C) 2008 by Bernd Flemisch, Andreas Lauser                    *
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
* \brief Reference implementation of a newton controller solver.
*
* Usually for most cases this controller should be sufficient.
*/
#ifndef DUMUX_NEWTON_CONTROLLER_HH
#define DUMUX_NEWTON_CONTROLLER_HH

#include <dumux/common/exceptions.hh>

#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include "dune/istl/owneroverlapcopy.hh"

#include <dune/istl/io.hh>

#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

#include <dumux/pardiso/pardiso.hh>

#if HAVE_DUNE_PDELAB
#include <dumux/boxmodels/pdelab/preconditionerpdelab.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#endif


namespace Dumux
{
namespace Properties
{
//! specifies the verbosity of the linear solver (by default it is 0,
//! i.e. it doesn't print anything)
NEW_PROP_TAG(NewtonLinearSolverVerbosity);

//! specifies whether the convergence rate and the global residual
//! gets written out to disk for every newton iteration (default is false)
NEW_PROP_TAG(NewtonWriteConvergence);

//! specifies whether the update should be done using the line search
//! method instead of the "raw" newton method. whether this property
//! has any effect depends on wether the line search method is
//! implemented for the actual model's newton controller's update()
//! method. By default we do not use line search.
NEW_PROP_TAG(NewtonUseLineSearch);

SET_PROP_DEFAULT(NewtonLinearSolverVerbosity)
{public:
	static const int value = 0;
};

SET_PROP_DEFAULT(NewtonWriteConvergence)
{public:
	static const bool value = false;
};

SET_PROP_DEFAULT(NewtonUseLineSearch)
{public:
	static const bool value = false;
};
};


template <class TypeTag, bool enable>
struct NewtonConvergenceWriter
{
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

	typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
	typedef typename SolutionTypes::SolutionFunction SolutionFunction;
	typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

	NewtonConvergenceWriter(NewtonController &ctl)
	: ctl_(ctl)
	{
		timeStepNum_ = 0;
		iteration_ = 0;
		vtkMultiWriter_ = new VtkMultiWriter("convergence");
	}

	~NewtonConvergenceWriter()
	{ delete vtkMultiWriter_; };

	void beginTimestep()
	{
		++timeStepNum_;
		iteration_ = 0;
	};

	void beginIteration(const GridView &gv)
	{
		++ iteration_;
		vtkMultiWriter_->beginTimestep(timeStepNum_ + iteration_ / 100.0,
				gv);
	};

	void writeFields(const SolutionFunction &uOld,
			const SolutionFunction &deltaU)
	{
		ctl_.model().localJacobian().addConvergenceVtkFields(*vtkMultiWriter_, uOld, deltaU);
	};

	void endIteration()
	{
		vtkMultiWriter_->endTimestep();
	};

	void endTimestep()
	{
		++timeStepNum_;
		iteration_ = 0;
	};

private:
	int timeStepNum_;
	int iteration_;
	VtkMultiWriter *vtkMultiWriter_;
	NewtonController &ctl_;
};

template <class TypeTag>
struct NewtonConvergenceWriter<TypeTag, false>
{ 
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) NewtonController;

	typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
	typedef typename SolutionTypes::SolutionFunction SolutionFunction;
	typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

	NewtonConvergenceWriter(NewtonController &ctl)
	{};

	void beginTimestep()
	{ };

	void beginIteration(const GridView &gv)
	{ };

	void writeFields(const SolutionFunction &uOld,
			const SolutionFunction &deltaU)
	{ };

	void endIteration()
	{ };

	void endTimestep()
	{ };
};

/*!
* \brief The reference implementation of a newton controller.
*
* If you want to specialize only some methods but are happy with
* the defaults of the reference controller, derive your
* controller from this class and simply overload the required
* methods.
*/
template <class TypeTag>
class NewtonController
{
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController))  Implementation;
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))    GridView;

	typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))        Model;
	typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;

	typedef typename GET_PROP(TypeTag, PTAG(PDELabTypes))  PDELabTypes;
	typedef typename PDELabTypes::GridFunctionSpace GridFunctionSpace;
	typedef typename PDELabTypes::ConstraintsTrafo ConstraintsTrafo;

	typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))  SolutionTypes;
	typedef typename SolutionTypes::SolutionFunction         SolutionFunction;

	typedef NewtonConvergenceWriter<TypeTag, GET_PROP_VALUE(TypeTag, PTAG(NewtonWriteConvergence))>  ConvergenceWriter;

public:
	NewtonController(Scalar tolerance  = 1e-7, // maximum tolerated deflection between two iterations
			int targetSteps   = 8,
			int maxSteps      = 12)
	: endIterMsgStream_(std::ostringstream::out),
	  convergenceWriter_(asImp_())
	{
		assert(maxSteps > targetSteps + 3);
		numSteps_ = 0;
		tolerance_ = tolerance;
		targetSteps_ = targetSteps;
		maxSteps_ = maxSteps;

		curPhysicalness_ = 0;
		maxPhysicalness_ = 0;
	};

	/*!
	* \brief Returns true iff another iteration should be done.
	*/
	bool newtonProceed(SolutionFunction &u)
	{
		if (numSteps_ < 2)
			return true; // we always do at least two iterations
		else if (numSteps_ >= maxSteps_) {
			// we have exceeded the allowed number of steps.  if the
			// relative error was reduced by a factor of at least 5,
			// we proceed even if we are above the maximum number of
			// steps
			return error_*4.0 < lastError_ && !asImp_().newtonConverged();
		}
		else if (asImp_().newtonConverged())
			return false; // we are below the desired tolerance

		Scalar tmp = asImp_().physicalness_(u);
		curPhysicalness_ = model().gridView().comm().min(tmp);
		curPhysicalness_ = std::min(curPhysicalness_, Scalar(1.0));


		// check for the physicalness of the solution
		if (curPhysicalness_ <= 0)
			// not physical enough even for a temporary
			// solution
			return false;
		else if (curPhysicalness_ < ((Scalar) numSteps_)/(maxSteps_ - 1)) {
			// we require that the solution gets more physical
			// with every step and at the last step the
			// solution must be completely physical.
			return false;
		}
		else if (curPhysicalness_ < maxPhysicalness_)
		{
			if (probationCount_ > 1) {
				// an iterative solution was more physical
				// than the current solution and at least 2
				// others.
				return false;
			}
			else {
				// we are physical enough, but some earlier
				// solution was more physical, so we let the
				// solver continue on probation.
				++probationCount_;
				return true;
			}
		}
		else {
			// everything's fine: the solution is physical
			// enough for the number of iterations we did and
			// it is the most physical so far.
			maxPhysicalness_ = curPhysicalness_;
			probationCount_ = std::max(0, probationCount_ - 1);

			return true; // do another round
		};
	}

	/*!
	* \brief Returns true iff the error of the solution is below the
	*        tolerance.
	*/
	bool newtonConverged()
	{
		return (error_ <= tolerance_) && (curPhysicalness_ >= 1.0);
	}

	/*!
	* \brief Called before the newton method is applied to an
	*        non-linear system of equations.
	*/
	void newtonBegin(NewtonMethod *method, SolutionFunction &u)
	{
		method_ = method;
		numSteps_ = 0;
		probationCount_ = 0;
		maxPhysicalness_ = 0;
		curPhysicalness_ = 0;

		convergenceWriter_.beginTimestep();
	}

	/*!
	* \brief Indidicates the beginning of a newton iteration.
	*/
	void newtonBeginStep()
	{ lastError_ = error_; }

	/*!
	* \brief Returns the number of steps done since newtonBegin() was
	*        called.
	*/
	int newtonNumSteps()
	{ return numSteps_; }


	/*!
	* \brief Update the error of the solution compared to the
	*        previous iteration.
	*/
	void newtonUpdateRelError(const SolutionFunction &uOld,
			const SolutionFunction &deltaU)
	{
		// calculate the relative error as the maximum relative
		// deflection in any degree of freedom.
		typedef typename SolutionFunction::BlockType FV;
		error_ = 0;
		for (int i = 0; i < int((*uOld).size()); ++i) {
			for (int j = 0; j < FV::size; ++j) {
				Scalar tmp
				=
						std::abs((*deltaU)[i][j])
				/ std::max(std::abs((*uOld)[i][j]), Scalar(1e-4));
				error_ = std::max(error_, tmp);
			}
		};
		error_ = model().gridView().comm().max(error_);
	}


	/*!
	* \brief Solve the linear system of equations \f$ \mathbf{A}x - b
	*        = 0\f$.
	*
	* Throws Dumux::NumericalProblem if the linear solver didn't
	* converge.
	*/
	template <class Matrix, class Vector>
	void newtonSolveLinear(Matrix &A,
			SolutionFunction &u,
			Vector &b)
	{
		// if the deflection of the newton method is large, we do not
		// need to solve the linear approximation accurately. Assuming
		// that the initial value for the delta vector u is quite
		// close to the final value, a reduction of 6 orders of
		// magnitute in the defect should be sufficient...
		Scalar residReduction = 1e-6;

		try {
#if HAVE_MPI
			solveParallel_(A, u, b, residReduction);
#else
			solveSequential_(A, *u, b, residReduction);
#endif
		}
		catch (Dune::MatrixBlockError e) {
			Dumux::NumericalProblem p;
			p.message(e.what());
			throw p;
		}
	}

	/*!
	* \brief Update the current solution function with a delta vector.
	*
	* The error estimates required for the newtonConverged() and
	* newtonProceed() methods should be updated here.
	*
	* Different update strategies, such as line search and chopped
	* updates can be implemented. The default behaviour is just to
	* subtract deltaU from uOld.
	*
	* \param deltaU The delta as calculated from solving the linear
	*               system of equations. This parameter also stores
	*               the updated solution.
	* \param uOld   The solution of the last iteration
	*/
	void newtonUpdate(SolutionFunction &deltaU, const SolutionFunction &uOld)
	{
		writeConvergence_(uOld, deltaU);

		newtonUpdateRelError(uOld, deltaU);

		*deltaU *= -1;
		*deltaU += *uOld;
	}

	/*!
	* \brief Indicates that one newton iteration was finished.
	*/
	void newtonEndStep(SolutionFunction &u, SolutionFunction &uOld)
	{
		++numSteps_;

		curPhysicalness_ = asImp_().physicalness_(u);
		if (this->method().verbose())
			std::cout << "\rNewton iteration " << numSteps_ << " done: "
			<< "error=" << error_ << endIterMsg().str() << "\n";
		endIterMsgStream_.str("");
	}

	/*!
	* \brief Indicates that we're done solving the non-linear system of equations.
	*/
	void newtonEnd()
	{
		convergenceWriter_.endTimestep();
	}

	/*!
	* \brief Called if the newton method broke down.
	*
	* This method is called _after_ newtonEnd()
	*/
	void newtonFail()
	{
		numSteps_ = targetSteps_*2;
	}

	/*!
	* \brief Called when the newton method was sucessful.
	*
	* This method is called _after_ newtonEnd()
	*/
	void newtonSucceed()
	{ }

	/*!
	* \brief Suggest a new time stepsize based on the old time step size.
	*
	* The default behaviour is to suggest the old time step size
	* scaled by the ratio between the target iterations and the
	* iterations required to actually solve the last time step.
	*/
	Scalar suggestTimeStepSize(Scalar oldTimeStep) const
	{
		// be agressive reducing the timestep size but
		// conservative when increasing it. the rationale is
		// that we want to avoid failing in the next newton
		// iteration which would require another linearization
		// of the problem.
		if (numSteps_ > targetSteps_) {
			Scalar percent = ((Scalar) numSteps_ - targetSteps_)/targetSteps_;
			return oldTimeStep/(1.0 + percent);
		}
		else {
			/*Scalar percent = (Scalar(1))/targetSteps_;
              return oldTimeStep*(1 + percent);
			 */
			Scalar percent = ((Scalar) targetSteps_ - numSteps_)/targetSteps_;
			return oldTimeStep*(1.0 + percent/1.2);
		}
	}

	/*!
	* \brief Returns a reference to the current newton method
	*        which is controlled by this controller.
	*/
	NewtonMethod &method()
	{ return *method_; }

	/*!
	* \brief Returns a reference to the current newton method
	*        which is controlled by this controller.
	*/
	const NewtonMethod &method() const
	{ return *method_; }

	/*!
	* \brief Returns a reference to the current numeric model.
	*/
	Model &model()
	{ return method_->model(); }

	/*!
	* \brief Returns a reference to the current numeric model.
	*/
	const Model &model() const
	{ return method_->model(); }

	std::ostringstream &endIterMsg()
	{ return endIterMsgStream_; }

protected:
	// returns the actual implementation for the cotroller we do
	// it this way in order to allow "poor man's virtual methods",
	// i.e. methods of subclasses which can be called by the base
	// class.
	Implementation &asImp_()
	{ return *static_cast<Implementation*>(this); }
	const Implementation &asImp_() const
	{ return *static_cast<const Implementation*>(this); }

	void writeConvergence_(const SolutionFunction &uOld,
			const SolutionFunction &deltaU)
	{
		convergenceWriter_.beginIteration(this->model().gridView());
		convergenceWriter_.writeFields(uOld, deltaU);
		convergenceWriter_.endIteration();
	};


#if HAVE_MPI
template <class Matrix, class Vector>
void solveParallel_(Matrix &A,
		SolutionFunction &u,
		Vector &b,
		Scalar residReduction)
{
	Vector &x = *u;

	int verbosity = GET_PROP_VALUE(TypeTag,
			PTAG(NewtonLinearSolverVerbosity));
	if (model().gridView().grid().comm().rank() != 0)
		verbosity = 0;

#ifdef HAVE_DUNE_PDELAB

#if HAVE_PARDISO
	typedef  Dune::PDELab::ISTLBackend_NoOverlap_Loop_Pardiso<TypeTag> Solver;
	Solver solver(model().jacobianAssembler().gridFunctionSpace(),
			model().jacobianAssembler().constraintsTrafo(), 5000, verbosity);
#else // !HAVE_PARDISO
//	typedef  Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GridFunctionSpace> Solver;
//	Solver solver(model().jacobianAssembler().gridFunctionSpace(), verbosity);
//	typedef  Dune::PDELab::ISTLBackend_NOVLP_BCGS_NOPREC<GridFunctionSpace> Solver;
//	Solver solver(model().jacobianAssembler().gridFunctionSpace(), 5000, verbosity);
	typedef  Dune::PDELab::ISTLBackend_NoOverlap_BCGS_ILU<TypeTag> Solver;
	Solver solver(model(), 5000, verbosity);
#endif // HAVE_PARDISO
	solver.apply(A, x, b, residReduction);

	if (!solver.result().converged)
		DUNE_THROW(Dumux::NumericalProblem,
				"Solving the linear system of equations did not converge.");

#else // !HAVE_DUNE_PDELAB

#if HAVE_PARDISO
	typedef Dune::SeqPardiso<Matrix,Vector,Vector> SeqPreCond;
	typedef Dune::LoopSolver<Vector> Solver;
	SeqPreCond seqPreCond(A);
#else // !HAVE_PARDISO
	typedef Dune::SeqILU0<Matrix,Vector,Vector> SeqPreCond;
	typedef Dune::BiCGSTABSolver<Vector> Solver;
	SeqPreCond seqPreCond(A, 0.9);
#endif

	typedef typename Grid::Traits::GlobalIdSet::IdType GlobalId;
	typedef Dune::OwnerOverlapCopyCommunication<GlobalId,int> Communication;
	Dune::IndexInfoFromGrid<GlobalId,int> indexinfo;
	u.fillIndexInfoFromGrid(indexinfo);
	Communication comm(indexinfo, MPIHelper::getCommunicator());
	Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,Communication> parallelOperator(A, comm);
	Dune::OverlappingSchwarzScalarProduct<Vector,Communication> scalarProduct(comm);
	Dune::BlockPreconditioner<Vector,Vector,Communication> parPreCond(seqPreCond, comm);
	Solver solver(parallelOperator,
			scalarProduct,
			parPreCond,
			residReduction,
			300,
			verbosity);

	Dune::InverseOperatorResult result;

	solver.apply(x, b, result);

	if (!solver.result().converged)
		DUNE_THROW(Dumux::NumericalProblem,
				"Solving the linear system of equations did not converge.");

#endif // HAVE_DUNE_PDELAB
}


#else // !HAVE_MPI
template <class Matrix, class Vector>
void solveSequential_(Matrix &A,
		Vector &x,
		Vector &b,
		Scalar residReduction)
{
	typedef typename SolutionTypes::JacobianAssembler::RepresentationType AsmRep;
	typedef typename SolutionFunction::RepresentationType FnRep;
	typedef Dune::MatrixAdapter<AsmRep, FnRep, FnRep>  Dune::MatrixAdapter;
	Dune::MatrixAdapter opA(A);

#ifdef HAVE_PARDISO
	Dune::SeqPardiso<Matrix,Vector,Vector> pardiso;
	pardiso.factorize(A);
	Dune::LoopSolver<Vector> solver(opA, pardiso, residReduction, 100, 0);         // an inverse operator
#else // HAVE_PARDISO
	// initialize the preconditioner
	Dune::SeqILU0<Matrix,Vector,Vector> precond(A, 1.0);
	//                Dune::SeqSSOR<OpAsmRep,FnRep,FnRep> precond(*opAsm, 3, 1.0);
	//                SeqIdentity<OpAsmRep,FnRep,FnRep> precond(*opAsm);
	// invert the linear equation system
	Dune::BiCGSTABSolver<Vector> solver(opA, precond, residReduction, 1000, 0);
#endif // HAVE_PARDISO

	Dune::InverseOperatorResult result;
	solver.apply(x, b, result);

	if (!result.converged)
		DUNE_THROW(Dumux::NumericalProblem,
				"Solving the linear system of equations did not converge.");
}
#endif // HAVE_MPI

//! this function is an indication of how "physically
//! meaningful" a temporary solution is. 0 means it isn't
//! meaningful at all (maybe because it has highly negative
//! pressures, etc) and the newton method can be stopped
//! immediately. Conversly 1 means that the solution is
//! perfectly physically meaningful (although it doesn't need
//! to be the solution in any way) and the method can to run.
//! Values inbetween mean that the funtion is not meaningfull,
//! but can be tolerated as temporary solution at some
//! iteration. (The controller assumes that as the method
//! progresses, the physicallness of the solution must
//! increase.)
Scalar physicalness_(SolutionFunction &u)
{
	return 1;
}

std::ostringstream endIterMsgStream_;

NewtonMethod *method_;

ConvergenceWriter convergenceWriter_;

Scalar tolerance_;

Scalar maxPhysicalness_;
Scalar curPhysicalness_;
Scalar error_;
Scalar lastError_;
int    probationCount_;

// optimal number of iterations we want to achive
int    targetSteps_;
// maximum number of iterations we do before giving up
int    maxSteps_;
// actual number of steps done so far
int    numSteps_;
};

}


#endif
