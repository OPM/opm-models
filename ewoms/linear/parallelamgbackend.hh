// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Ewoms::Linear::ParallelAmgBackend
 */
#ifndef EWOMS_PARALLEL_AMG_BACKEND_HH
#define EWOMS_PARALLEL_AMG_BACKEND_HH

#include <ewoms/linear/paralleliterativebackend.hh>
#include <ewoms/linear/superlubackend.hh>
#include <ewoms/linear/vertexborderlistfromgrid.hh>
#include <ewoms/linear/overlappingbcrsmatrix.hh>
#include <ewoms/linear/overlappingblockvector.hh>
#include <ewoms/linear/overlappingpreconditioner.hh>
#include <ewoms/linear/overlappingscalarproduct.hh>
#include <ewoms/linear/overlappingoperator.hh>
#include <ewoms/linear/solverpreconditioner.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dune/common/parallel/indexset.hh>
#include <dune/common/version.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <iostream>

namespace Ewoms {
namespace Linear {
template <class TypeTag>
class ParallelAmgBackend;
}

namespace Properties {
NEW_TYPE_TAG(ParallelAmgLinearSolver, INHERITS_FROM(ParallelIterativeLinearSolver));

NEW_PROP_TAG(AmgCoarsenTarget);

//! The target number of DOFs per processor for the parallel algebraic
//! multi-grid solver
SET_INT_PROP(ParallelAmgLinearSolver, AmgCoarsenTarget, 5000);

SET_TYPE_PROP(ParallelAmgLinearSolver, LinearSolverBackend,
              Ewoms::Linear::ParallelAmgBackend<TypeTag>);

SET_TYPE_PROP(ParallelAmgLinearSolver, LinearSolverScalar,
              typename GET_PROP_TYPE(TypeTag, Scalar));
} // namespace Properties

namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides a linear solver backend using the parallel
 *        algebraic multi-grid (AMG) linear solver from DUNE-ISTL.
 */
template <class TypeTag>
class ParallelAmgBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverScalar) LinearSolverScalar;
    typedef typename GET_PROP_TYPE(TypeTag, BorderListCreator) BorderListCreator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    typedef Dune::FieldVector<LinearSolverScalar, numEq> VectorBlock;
    typedef Dune::FieldMatrix<LinearSolverScalar, numEq, numEq> MatrixBlock;

    typedef Dune::BCRSMatrix<MatrixBlock> Matrix;
    typedef Dune::BlockVector<VectorBlock> Vector;

    // define the smoother used for the AMG and specify its
    // arguments
    typedef Dune::SeqSOR<Matrix, Vector, Vector> SequentialSmoother;
// typedef Dune::SeqSSOR<Matrix,Vector,Vector> SequentialSmoother;
// typedef Dune::SeqJac<Matrix,Vector,Vector> SequentialSmoother;
// typedef Dune::SeqILU0<Matrix,Vector,Vector> SequentialSmoother;
// typedef Dune::SeqILUn<Matrix,Vector,Vector> SequentialSmoother;

#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<Ewoms::Linear::Index>
    OwnerOverlapCopyCommunication;
    typedef Dune::OverlappingSchwarzOperator<Matrix,
                                             Vector,
                                             Vector,
                                             OwnerOverlapCopyCommunication> FineOperator;
    typedef Dune::OverlappingSchwarzScalarProduct<Vector,
                                                  OwnerOverlapCopyCommunication> FineScalarProduct;
    typedef Dune::BlockPreconditioner<Vector,
                                      Vector,
                                      OwnerOverlapCopyCommunication,
                                      SequentialSmoother> ParallelSmoother;
    typedef Dune::Amg::AMG<FineOperator,
                           Vector,
                           ParallelSmoother,
                           OwnerOverlapCopyCommunication> AMG;
#else
    typedef Dune::MatrixAdapter<Matrix, Vector, Vector> FineOperator;
    typedef Dune::SeqScalarProduct<Vector> FineScalarProduct;
    typedef SequentialSmoother ParallelSmoother;
    typedef Dune::Amg::AMG<FineOperator, Vector, ParallelSmoother> AMG;
#endif

public:
    ParallelAmgBackend(const Simulator& simulator)
        : simulator_(simulator)
    {
        overlappingMatrix_ = nullptr;
        overlappingb_ = nullptr;
        overlappingx_ = nullptr;

        amg_ = nullptr;
    }

    ~ParallelAmgBackend()
    { cleanup_(); }

    static void registerParameters()
    {
        ParallelIterativeSolverBackend<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, int, AmgCoarsenTarget,
                             "The coarsening target for the agglomerations of "
                             "the AMG preconditioner");
    }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     */
    void eraseMatrix()
    { cleanup_(); }

    template <class NativeBCRSMatrix>
    void prepareMatrix(const NativeBCRSMatrix& M)
    {
        if (!overlappingMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        }

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlappingMatrix_->assignAdd(M);
    }

    template <class NativeBCRSMatrix, class NativeVector>
    void prepareRhs(const NativeBCRSMatrix& M, NativeVector& b)
    {
        if (!overlappingMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        }

        overlappingb_->assignAddBorder(b);

        // copy the result back to the non-overlapping vector. This is
        // necessary here as assignAddBorder() might modify the
        // residual vector for the border entities and we need the
        // "globalized" residual in b...
        overlappingb_->assignTo(b);
    }

    /*!
     * \brief Actually solve the linear system of equations.
     *
     * \return true if the residual reduction could be achieved, else false.
     */
    template <class NativeVector>
    bool solve(NativeVector& x)
    {
        int verbosity = 0;
        if (simulator_.gridManager().gridView().comm().rank() == 0)
            verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);

        (*overlappingx_) = 0.0;

        /////////////
        // set-up the AMG preconditioner
        /////////////
        // create the parallel scalar product and the parallel operator
#if HAVE_MPI
        FineOperator fineOperator(*overlappingMatrix_, *istlComm_);
#else
        FineOperator fineOperator(*overlappingMatrix_);
#endif

#if HAVE_MPI
        FineScalarProduct scalarProduct(*istlComm_);
#else
        FineScalarProduct scalarProduct;
#endif

        setupAmg_(fineOperator);

        /////////////
        // set-up the linear solver
        /////////////
        if (verbosity > 1 && simulator_.gridManager().gridView().comm().rank() == 0)
            std::cout << "Creating the solver\n" << std::flush;

        typedef Ewoms::BiCGSTABSolver<Vector> SolverType;
        Scalar linearSolverTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        int maxIterations = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations);
        SolverType solver(fineOperator,
                          scalarProduct,
                          /*preconditioner=*/*amg_,
                          linearSolverTolerance,
                          /*maxSteps=*/maxIterations,
                          verbosity);

        /////
        // create a residual reduction convergence criterion

        // set the weighting of the residuals
        OverlappingVector residWeightVec(*overlappingx_);
        residWeightVec = 0.0;
        const auto& overlap = overlappingMatrix_->overlap();
        for (Index localIdx = 0; localIdx < static_cast<Index>(overlap.numLocal()); ++localIdx) {
            Index nativeIdx = overlap.domesticToNative(localIdx);
            for (unsigned eqIdx = 0; eqIdx < Vector::block_type::dimension; ++eqIdx) {
                residWeightVec[static_cast<unsigned>(localIdx)][eqIdx] =
                    this->simulator_.model().eqWeight(static_cast<unsigned>(nativeIdx), eqIdx);
            }
        }

        Scalar linearSolverAbsTolerance = simulator_.model().newtonMethod().tolerance() / 100.0;
        typedef typename GridView::CollectiveCommunication Comm;
        auto *convCrit =
            new Ewoms::WeightedResidualReductionCriterion<Vector, Comm>(
                simulator_.gridView().comm(),
                residWeightVec,
                /*residualReductionTolerance=*/linearSolverTolerance,
                /*fixPointTolerance=*/0.0,
                /*absoluteResidualTolerance=*/linearSolverAbsTolerance);

        // done creating the convergence criterion
        /////

        // tell the linear solver to use it
        typedef Ewoms::ConvergenceCriterion<Vector> ConvergenceCriterion;
        solver.setConvergenceCriterion(std::shared_ptr<ConvergenceCriterion>(convCrit));

        Dune::InverseOperatorResult result;
        int solverSucceeded = 1;
        try
        {
            solver.apply(*overlappingx_, *overlappingb_, result);
            solverSucceeded = simulator_.gridManager().gridView().comm().min(solverSucceeded);
        }
        catch (const Dune::Exception& )
        {
            solverSucceeded = 0;
            solverSucceeded = simulator_.gridManager().gridView().comm().min(solverSucceeded);
        }

        if (!solverSucceeded)
            return false;

        // copy the result back to the non-overlapping vector
        overlappingx_->assignTo(x);

        return result.converged;
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    template <class NativeBCRSMatrix>
    void prepare_(const NativeBCRSMatrix& M)
    {
        BorderListCreator borderListCreator(simulator_.gridView(),
                                            simulator_.model().dofMapper());

        auto& blackList = borderListCreator.blackList();

        // create the overlapping Jacobian matrix
        unsigned overlapSize = EWOMS_GET_PARAM(TypeTag, unsigned, LinearSolverOverlapSize);
        overlappingMatrix_ = new OverlappingMatrix(M,
                                                   borderListCreator.borderList(),
                                                   blackList,
                                                   overlapSize);

        // create the overlapping vectors for the residual and the
        // solution
        overlappingb_ = new OverlappingVector(overlappingMatrix_->overlap());
        overlappingx_ = new OverlappingVector(*overlappingb_);

        // writeOverlapToVTK_();

#if HAVE_MPI
        // create and initialize DUNE's OwnerOverlapCopyCommunication
        // using the domestic overlap
        istlComm_ = new OwnerOverlapCopyCommunication(MPI_COMM_WORLD);
        setupAmgIndexSet_(overlappingMatrix_->overlap(), istlComm_->indexSet());
        istlComm_->remoteIndices().template rebuild<false>();
#endif
    }

    void cleanup_()
    {
        // create the overlapping Jacobian matrix and vectors
        delete overlappingMatrix_;
        delete overlappingb_;
        delete overlappingx_;
        delete amg_;

        overlappingMatrix_ = nullptr;
        overlappingb_ = nullptr;
        overlappingx_ = nullptr;
        amg_ = nullptr;
    }

#if HAVE_MPI
    template <class ParallelIndexSet>
    void setupAmgIndexSet_(const Overlap& overlap, ParallelIndexSet& istlIndices)
    {
        typedef Dune::OwnerOverlapCopyAttributeSet GridAttributes;
        typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet GridAttributeSet;

        // create DUNE's ParallelIndexSet from a domestic overlap
        istlIndices.beginResize();
        for (Index curIdx = 0; static_cast<size_t>(curIdx) < overlap.numDomestic(); ++curIdx) {
            GridAttributeSet gridFlag =
                overlap.iAmMasterOf(curIdx)
                ? GridAttributes::owner
                : GridAttributes::copy;

            // an index is used by other processes if it is in the
            // domestic or in the foreign overlap.
            bool isShared = overlap.isInOverlap(curIdx);

            assert(curIdx == overlap.globalToDomestic(overlap.domesticToGlobal(curIdx)));
            istlIndices.add(/*globalIdx=*/overlap.domesticToGlobal(curIdx),
                            Dune::ParallelLocalIndex<GridAttributeSet>(static_cast<size_t>(curIdx),
                                                                       gridFlag,
                                                                       isShared));
        }
        istlIndices.endResize();
    }
#endif

    void setupAmg_(FineOperator& fineOperator)
    {
        if (amg_)
            return;

        int verbosity = 0;
        if (simulator_.gridManager().gridView().comm().rank() == 0)
            verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);

        auto rank = simulator_.gridManager().gridView().comm().rank();
        if (verbosity > 1 && rank == 0)
            std::cout << "Setting up the AMG preconditioner\n";

        typedef typename Dune::Amg::SmootherTraits<ParallelSmoother>::Arguments
        SmootherArgs;

        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1.0;

        // specify the coarsen criterion:
        //
        // typedef
        // Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<Matrix,
        //                             Dune::Amg::FirstDiagonal>>
        typedef Dune::Amg::
            CoarsenCriterion<Dune::Amg::SymmetricCriterion<Matrix, Dune::Amg::FrobeniusNorm> >
            CoarsenCriterion;
        int coarsenTarget = EWOMS_GET_PARAM(TypeTag, int, AmgCoarsenTarget);
        CoarsenCriterion coarsenCriterion(/*maxLevel=*/15, coarsenTarget);
        coarsenCriterion.setDefaultValuesAnisotropic(GridView::dimension,
                                                     /*aggregateSizePerDim=*/3);
        if (verbosity > 0)
            coarsenCriterion.setDebugLevel(1);
        else
            coarsenCriterion.setDebugLevel(0); // make the AMG shut up

        // reduce the minium coarsen rate (default is 1.2)
        coarsenCriterion.setMinCoarsenRate(1.05);
        // coarsenCriterion.setAccumulate(Dune::Amg::noAccu);
        coarsenCriterion.setAccumulate(Dune::Amg::atOnceAccu);
        coarsenCriterion.setSkipIsolated(false);

// instantiate the AMG preconditioner
#if HAVE_MPI
        amg_ = new AMG(fineOperator, coarsenCriterion, smootherArgs, *istlComm_);
#else
        amg_ = new AMG(fineOperator, coarsenCriterion, smootherArgs);
#endif
    }

    const Simulator& simulator_;

    AMG *amg_;

#if HAVE_MPI
    OwnerOverlapCopyCommunication *istlComm_;
#endif
    OverlappingMatrix *overlappingMatrix_;
    OverlappingVector *overlappingb_;
    OverlappingVector *overlappingx_;
};

} // namespace Linear
} // namespace Ewoms

#endif
