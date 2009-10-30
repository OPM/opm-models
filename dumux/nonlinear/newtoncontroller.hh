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
#ifndef DUNE_NEWTON_CONTROLLER_HH
#define DUNE_NEWTON_CONTROLLER_HH

#include "config.h"

#include <dumux/exceptions.hh>

#include <dune/istl/overlappingschwarz.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include "dune/istl/owneroverlapcopy.hh"

#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

#include <dumux/pardiso/pardiso.hh>

#if HAVE_DUNE_PDELAB
#include <dumux/boxmodels/pdelab/preconditionerpdelab.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#endif


namespace Dune
{
/*!
 * \brief Base class for the reference implementation of a newton
 *        controller.
 *
 * If you want to specialize only some methods but are happy with
 * the defaults of the reference controller, derive your
 * controller from this class and simply overload the required
 * methods.
 */
template <class NewtonMethod,
          class Implementation>
class NewtonControllerBase
{
    typedef typename NewtonMethod::Model  Model;
    typedef typename Model::NewtonTraits  NewtonTraits;

public:
    typedef typename NewtonTraits::Scalar  Scalar;
    typedef typename NewtonTraits::Grid    Grid;

    typedef typename NewtonTraits::Function              Function;
    typedef typename NewtonTraits::JacobianAssembler     JacobianAssembler;
#if HAVE_DUNE_PDELAB
    typedef typename JacobianAssembler::GridFunctionSpace GridFunctionSpace;
    typedef typename JacobianAssembler::ConstraintsTrafo ConstraintsTrafo;
#endif
    typedef typename JacobianAssembler::RepresentationType    JacAsmRep;

    NewtonControllerBase(Scalar tolerance, // maximum tolerated deflection between two iterations
                         int targetSteps,
                         int maxSteps)
    {
        assert(maxSteps > targetSteps + 3);
        numSteps_ = 0;
        tolerance_ = tolerance;
        targetSteps_ = targetSteps;
        maxSteps_ = maxSteps;

        curPhysicalness_ = 0;
        maxPhysicalness_ = 0;
    };

    //! Returns true if another iteration should be done.
    bool newtonProceed(Function &u)
    {
        if (numSteps_ < 2)
            return true; // we always do at least two iterations
        else if (numSteps_ > maxSteps_)
            return false; // we have exceeded the allowed number of steps
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

    //! Returns true if the defect of the solution is below the
    //! tolerance
    bool newtonConverged()
    {
        return (defect_ <= tolerance_) && (curPhysicalness_ >= 1.0);
    }

    //! called before the newton method is applied to an equation
    //! system.
    void newtonBegin(NewtonMethod *method, Function &u)
    {
        method_ = method;
        numSteps_ = 0;
        probationCount_ = 0;
        maxPhysicalness_ = 0;
        curPhysicalness_ = 0;
    }

    //! indidicates the beginning of a newton iteration
    void newtonBeginStep()
    {
    }

    //! Returns the number of steps done since newtonBegin() was
    //! called
    int newtonNumSteps()
    { return numSteps_; }

    
    /*!
     * \brief Update the defect of the solution compared to the
     *        previous iteration.
     */
    template <class Function>
    void newtonUpdateRelDefect(const Function &uOld,
                               const Function &deltaU)
    {
        // calculate the relative defect as the maximum relative
        // deflection in any degree of freedom.
        typedef typename Function::BlockType FV;
        defect_ = 0;
        for (int i = 0; i < (*uOld).size(); ++i) {
            for (int j = 0; j < FV::size; ++j) {
                Scalar tmp
                    = 
                    std::abs((*deltaU)[i][j])
                    / std::max(std::abs((*uOld)[i][j]), Scalar(1));
                defect_ = std::max(defect_, tmp);
            }
        };
        model().gridView().comm().max(defect_);

    };

    Scalar relDefect() const
    { return defect_; };

    //! Solve the linear equation system Ax - b = 0 for the
    //! current iteration.
    //! Returns true iff the equation system could be solved.
    template <class Matrix, class Function, class Vector>
    void newtonSolveLinear(Matrix &A,
                           Function &u,
                           Vector &b)
    {
#if HAVE_MPI
        solveParallel_(A, u, b);
#else
        solveSequential_(A, *u, b);
#endif
    }

    //! Indicates that we're done solving one newton step.
    void newtonEndStep(Function &u, Function &uOld)
    {
        ++numSteps_;

/*        Scalar tmp = (*u).two_norm2();
        tmp = sqrt(model().gridView().comm().sum(tmp));
        oneByMagnitude_ = 1.0/std::max(tmp, Scalar(1e-5));
*/

        curPhysicalness_ = asImp_().physicalness_(u);
        if (this->method().verbose())
            std::cout << boost::format("\rNewton iteration %d done: defect=%g, physicalness: %.3f, maxPhysicalness=%.3f\n")
                %numSteps_%defect_%curPhysicalness_%maxPhysicalness_;
    }

    //! Indicates that we're done solving the equation system.
    void newtonEnd()
    {}

    //! Called when the newton method broke down.
    void newtonFail()
    {
        numSteps_ = targetSteps_*2;
    }

    //! Suggest a new time stepsize based on the number of newton
    //! iterations required for the last time step and the old time
    //! step size.
    Scalar suggestTimeStepSize(Scalar oldTimeStep) const
    {
        // be agressive reducing the timestep size but
        // conservative when increasing it. the rationale is
        // that we want to avoid failing in the next newton
        // iteration which would require another linerization
        // of the problem.
        if (numSteps_ > targetSteps_) {
            Scalar percent = ((Scalar) numSteps_ - targetSteps_)/targetSteps_;
            return oldTimeStep/(1 + percent);
        }
        else {
            Scalar percent = ((Scalar) targetSteps_ - numSteps_)/targetSteps_;
            return oldTimeStep*(1 + percent/1.2);
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

    template <class Vector>
    Scalar weightedNorm2(const Vector& vector, const Vector& sol)
    {
    	const unsigned int numEq = vector[0].size;

    	std::vector<Scalar> normsVSquared(numEq, 0);
    	std::vector<Scalar> normsSolSquared(numEq, 0);

    	for (unsigned int block = 0; block < vector.size(); block++)
    		for (unsigned int comp = 0; comp < numEq; comp++)
    		{
    			normsVSquared[comp] += vector[block][comp]*vector[block][comp];
    			normsSolSquared[comp] += sol[block][comp]*sol[block][comp];
    		}

    	Scalar result = 0;
		for (unsigned int comp = 0; comp < numEq; comp++)
			result += normsVSquared[comp]/std::max(normsSolSquared[comp], 1e-10);

		return result;
    }

protected:
    // returns the actual implementation for the cotroller we do
    // it this way in order to allow "poor man's virtual methods",
    // i.e. methods of subclasses which can be called by the base
    // class.
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }


#if HAVE_MPI
    template <class Matrix, class Vector>
    void solveParallel_(Matrix &A,
                        Function &u,
                        Vector &b)
    {
        Vector &x = *u;

        // if the deflection of the newton method is large, we
        // do not need to solve the linear approximation
        // accurately. On the other hand, if this is the first
        // newton step, we don't have a meaningful value for the defect
        // yet, so we use the targeted accurracy for the defect.
        Scalar residTol = tolerance_/1e8;

#if 0//HAVE_PARDISO
    	typedef Dune::SeqPardiso<Matrix,Vector,Vector> SeqPreCond;
    	SeqPreCond seqPreCond(A);
#else
    	typedef Dune::SeqILU0<Matrix,Vector,Vector> SeqPreCond;
    	SeqPreCond seqPreCond(A, 1.0);
#endif

#ifdef HAVE_DUNE_PDELAB
    	typedef Dune::PDELab::ParallelISTLHelper<GridFunctionSpace> ParallelHelper;
    	ParallelHelper parallelHelper(model().jacobianAssembler().gridFunctionSpace());
    	typedef Dune::PDELab::NonoverlappingOperator<GridFunctionSpace,Matrix,Vector,Vector> ParallelOperator;
    	ParallelOperator parallelOperator(model().jacobianAssembler().gridFunctionSpace(), A, parallelHelper);
    	typedef Dune::PDELab::NonoverlappingScalarProduct<GridFunctionSpace,Vector> ParallelScalarProduct;
    	ParallelScalarProduct scalarProduct(model().jacobianAssembler().gridFunctionSpace(), parallelHelper);
    	typedef Dune::PDELab::NonoverlappingWrappedPreconditioner<ConstraintsTrafo, GridFunctionSpace, SeqPreCond> ParPreCond;
    	ParPreCond parPreCond(model().jacobianAssembler().gridFunctionSpace(), seqPreCond,
    								  model().jacobianAssembler().constraintsTrafo(), parallelHelper);
//    	typedef Dune::PDELab::ParallelISTLHelper<GridFunctionSpace> ParallelHelper;
//    	ParallelHelper parallelHelper(model().jacobianAssembler().gridFunctionSpace());
//    	typedef Dune::PDELab::OverlappingOperator<ConstraintsTrafo,Matrix,Vector,Vector> ParallelOperator;
//    	ParallelOperator parallelOperator(model().jacobianAssembler().constraintsTrafo(), A);
//    	typedef Dune::PDELab::OverlappingScalarProduct<GridFunctionSpace,Vector> ParallelScalarProduct;
//    	ParallelScalarProduct scalarProduct(model().jacobianAssembler().gridFunctionSpace(), parallelHelper);
//    	typedef Dune::PDELab::OverlappingWrappedPreconditioner<ConstraintsTrafo, GridFunctionSpace, SeqPreCond> ParPreCond;
//    	ParPreCond parPreCond(model().jacobianAssembler().gridFunctionSpace(), seqPreCond,
//    								  model().jacobianAssembler().constraintsTrafo(), parallelHelper);
#else
        typedef typename Grid::Traits::GlobalIdSet::IdType GlobalId;
        typedef Dune::OwnerOverlapCopyCommunication<GlobalId,int> Communication;
        Dune::IndexInfoFromGrid<GlobalId,int> indexinfo;
        u.fillIndexInfoFromGrid(indexinfo);
        Communication comm(indexinfo, MPIHelper::getCommunicator());
        Dune::OverlappingSchwarzOperator<Matrix,Vector,Vector,Communication> parallelOperator(A, comm);
        Dune::OverlappingSchwarzScalarProduct<Vector,Communication> scalarProduct(comm);
        Dune::BlockPreconditioner<Vector,Vector,Communication> parPreCond(seqPreCond, comm);
#endif

    	Dune::BiCGSTABSolver<Vector>
            solver(parallelOperator,
                   scalarProduct,
                   parPreCond,
                   residTol,
                   1000,
                   0);

        Dune::InverseOperatorResult result;

        solver.apply(x, b, result);

        if (!result.converged)
            DUNE_THROW(Dune::NumericalProblem,
                       "Solving the linear system of equations did not converge.");
    }

#elif !HAVE_MPI
    template <class Matrix, class Vector>
    void solveSequential_(Matrix &A,
                          Vector &x,
                          Vector &b)
    {
        // if the deflection of the newton method is large, we
        // do not need to solve the linear approximation
        // accurately. On the other hand, if this is the first
        // newton step, we don't have a meaningful value for the defect
        // yet, so we use the targeted accurracy for the defect.
        Scalar residTol = tolerance_/1e8;

        typedef Dune::MatrixAdapter<typename JacobianAssembler::RepresentationType,
            typename Function::RepresentationType,
            typename Function::RepresentationType>  MatrixAdapter;
        MatrixAdapter opA(A);

#ifdef HAVE_PARDISO
        SeqPardiso<Matrix,Vector,Vector> pardiso;
        pardiso.factorize(A);
        BiCGSTABSolver<Vector> solver(opA, pardiso, residTol, 100, 2);         // an inverse operator
#else // HAVE_PARDISO
        // initialize the preconditioner
        Dune::SeqILU0<Matrix,Vector,Vector> precond(A, 1.0);
        //                Dune::SeqSSOR<OpAsmRep,FnRep,FnRep> precond(*opAsm, 3, 1.0);
        //                SeqIdentity<OpAsmRep,FnRep,FnRep> precond(*opAsm);
        // invert the linear equation system
        Dune::BiCGSTABSolver<Vector> solver(opA, precond, residTol, 500, 0);
#endif // HAVE_PARDISO

        Dune::InverseOperatorResult result;
        solver.apply(x, b, result);

        if (!result.converged)
            DUNE_THROW(Dune::NumericalProblem,
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
    Scalar physicalness_(Function &u)
    {
        return 1;
    }

    NewtonMethod *method_;

    Scalar tolerance_;

    Scalar maxPhysicalness_;
    Scalar curPhysicalness_;
    Scalar defect_;
    int    probationCount_;

    // optimal number of iterations we want to achive
    int    targetSteps_;
    // maximum number of iterations we do before giving up
    int    maxSteps_;
    // actual number of steps done so far
    int    numSteps_;
};

//! A reference implementation of a newton method controller
//!
//! Basically the only difference to NewtonControllerBase is that
//! this class can be instanciated more easily.
template <class NewtonMethod>
class NewtonController
    : public NewtonControllerBase<NewtonMethod, NewtonController<NewtonMethod> >
{
public:
    typedef NewtonController<NewtonMethod>               ThisType;
    typedef NewtonControllerBase<NewtonMethod, ThisType> ParentType;

    typedef typename ParentType::Scalar            Scalar;
    typedef typename ParentType::Function          Function;
    typedef typename ParentType::JacobianAssembler JacobianAssembler;

    NewtonController(Scalar tolerance = 1e-4, int targetSteps=8, int maxSteps = 12)
        : ParentType(tolerance, targetSteps, maxSteps)
    {};
};
}


#endif
