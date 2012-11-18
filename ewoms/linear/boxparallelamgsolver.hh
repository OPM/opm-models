// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Ewoms::Linear::BoxParallelAmgSolver
 */
#ifndef EWOMS_BOX_PARALLEL_AMG_SOLVER_HH
#define EWOMS_BOX_PARALLEL_AMG_SOLVER_HH

#include <ewoms/linear/seqsolverbackend.hh>
#include <ewoms/linear/vertexborderlistfromgrid.hh>
#include <ewoms/linear/overlappingbcrsmatrix.hh>
#include <ewoms/linear/overlappingblockvector.hh>
#include <ewoms/linear/overlappingpreconditioner.hh>
#include <ewoms/linear/overlappingscalarproduct.hh>
#include <ewoms/linear/overlappingoperator.hh>
#include <ewoms/linear/solverpreconditioner.hh>
#include <ewoms/linear/boxlinearsolver.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parameters.hh>

#include <ewoms/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dune/common/parallel/indexset.hh>
#include <dune/common/mpicollectivecommunication.hh>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(AmgCoarsenTarget);

//! The target number of DOFs per processor for the parallel algebraic
//! multi-grid solver
SET_INT_PROP(BoxLinearSolverTypeTag, AmgCoarsenTarget, 5000);
}

namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides a parallel algebraic multi-grid (AMG) linear solver
 *        from DUNE-ISTL for the vertex-centered finite volume ("box")
 *        method.
 */
template <class TypeTag>
class BoxParallelAmgSolver
{
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;

    enum { dimWorld = GridView::dimensionworld };

    // define the smoother used for the AMG and specify its
    // arguments
    typedef Dune::SeqSOR<Matrix,Vector,Vector> SequentialSmoother;
    //typedef Dune::SeqSSOR<Matrix,Vector,Vector> SequentialSmoother;
    //typedef Dune::SeqJac<Matrix,Vector,Vector> SequentialSmoother;
    //typedef Dune::SeqILU0<Matrix,Vector,Vector> SequentialSmoother;
    //typedef Dune::SeqILUn<Matrix,Vector,Vector> SequentialSmoother;

#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<Ewoms::Linear::Index> OwnerOverlapCopyCommunication;
    typedef Dune::OverlappingSchwarzOperator<Matrix,
                                             Vector,
                                             Vector,
                                             OwnerOverlapCopyCommunication> FineOperator;
    typedef Dune::OverlappingSchwarzScalarProduct<Vector,OwnerOverlapCopyCommunication> FineScalarProduct;
    typedef Dune::BlockPreconditioner<Vector,
                                      Vector,
                                      OwnerOverlapCopyCommunication,
                                      SequentialSmoother> ParallelSmoother;
    typedef Dune::Amg::AMG<FineOperator, Vector, ParallelSmoother, OwnerOverlapCopyCommunication> AMG;
#else
    typedef Dune::MatrixAdapter<Matrix, Vector, Vector> FineOperator;
    typedef Dune::SeqScalarProduct<Vector> FineScalarProduct;
    typedef SequentialSmoother ParallelSmoother;
    typedef Dune::Amg::AMG<FineOperator, Vector, ParallelSmoother> AMG;
#endif

public:
    BoxParallelAmgSolver(const Problem &problem)
        : problem_(problem)
    {
        overlappingMatrix_ = nullptr;
        overlappingb_ = nullptr;
        overlappingx_ = nullptr;

        amg_ = nullptr;
    }

    ~BoxParallelAmgSolver()
    { cleanup_(); }

    static void registerParameters()
    {
        BoxParallelSolver<TypeTag>::registerParameters();

        REGISTER_PARAM(TypeTag, int, AmgCoarsenTarget, "The coarsening target for the agglomerations of the AMG preconditioner");
    }

    /*!
     * \brief Set the structure of the linear system of equations to be solved.
     *
     * This method allocates space an does the necessary
     * communication before actually calling the solve() method.  As
     * long as the structure of the linear system does not change, the
     * solve method can be called arbitrarily often.
     */
    void setStructureMatrix(const Matrix &M)
    {
        cleanup_();
        prepare_();
    }

    /*!
     * \brief Actually solve the linear system of equations.
     *
     * \return true if the residual reduction could be achieved, else false.
     */
    bool solve(const Matrix &M, Vector &x, const Vector &b)
    {
        int verbosity = 0;
        if (problem_.gridView().comm().rank() == 0)
            verbosity = GET_PARAM(TypeTag, int, LinearSolverVerbosity);

        /////////////
        // set-up the overlapping matrix and vector
        /////////////

        if (!overlappingMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        };

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlappingMatrix_->assignAdd(M);
        overlappingMatrix_->resetFront();

        overlappingb_->assignAddBorder(b);
        overlappingb_->resetFront();

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
        if (verbosity > 1 && problem_.gridView().comm().rank() == 0)
            std::cout << "Creating the solver\n";

        typedef Ewoms::BiCGSTABSolver<Vector>  SolverType;
        Scalar linearSolverTolerance = GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        int maxIterations = GET_PARAM(TypeTag, Scalar, LinearSolverMaxIterations);
        SolverType solver(fineOperator,
                          scalarProduct,
                          /*preconditioner=*/*amg_,
                          linearSolverTolerance,
                          /*maxSteps=*/maxIterations,
                          verbosity);

        // use the fixpoint convergence criterion if possible
        Vector weightVec(*overlappingx_);

        // set the default weight of the row to 0 (-> do not consider
        // the row when calculating the error)
        weightVec = 0.0;
        // for rows local to the current peer, ping the model for their
        // relative weight
        const auto &foreignOverlap = overlappingMatrix_->overlap().foreignOverlap();
        for (int localIdx = 0; localIdx < foreignOverlap.numLocal(); ++localIdx) {
            int nativeIdx = foreignOverlap.localToNative(localIdx);
            for (int eqIdx = 0; eqIdx < Vector::block_type::dimension; ++eqIdx) {
                weightVec[localIdx][eqIdx] = this->problem_.model().eqWeight(nativeIdx, eqIdx);
            }
        }

        // create a weighted residual reduction convergence criterion
        Scalar linearSolverAbsTolerance = GET_PARAM(TypeTag, Scalar, LinearSolverAbsTolerance);
        auto *convCrit =
            new Ewoms::WeightedResidReductionCriterion<Vector,
                                                       typename GridView::CollectiveCommunication>
            (problem_.gridView().comm(),
             weightVec,
             linearSolverTolerance,
             linearSolverAbsTolerance);

        typedef Ewoms::ConvergenceCriterion<Vector> ConvergenceCriterion;
        solver.setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(convCrit));

        /////////////
        // run the linear solver
        /////////////
        Ewoms::InverseOperatorResult solverStatistics;
        solver.apply(*overlappingx_,*overlappingb_,solverStatistics);

        // copy the result back to the non-overlapping vector
        overlappingx_->assignTo(x);

        return solverStatistics.converged;
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    void prepare_(const Matrix &M)
    {
        int overlapSize = GET_PARAM(TypeTag, int, LinearSolverOverlapSize);

        int verbosity = 0;
        if (problem_.gridView().comm().rank() == 0) {
            verbosity = GET_PARAM(TypeTag, int, LinearSolverVerbosity);

            if (verbosity > 1)
                std::cout << "Creating algebraic overlap of size " << overlapSize << "\n";
        }

        Linear::VertexBorderListFromGrid<GridView, VertexMapper>
            borderListCreator(problem_.gridView(), problem_.vertexMapper());

        // blacklist all entries that belong to ghost and overlap vertices
        std::set<Ewoms::Linear::Index> blackList;
        auto vIt = problem_.gridView().template begin<dimWorld>();
        const auto &vEndIt = problem_.gridView().template end<dimWorld>();
        for (; vIt != vEndIt; ++vIt) {
            if (vIt->partitionType() != Dune::InteriorEntity &&
                vIt->partitionType() != Dune::BorderEntity)
            {
                // we blacklist everything except interior and border
                // vertices
                blackList.insert(problem_.vertexMapper().map(*vIt));
            }
        }

        // create the overlapping Jacobian matrix
        overlappingMatrix_ = new OverlappingMatrix (M,
                                                borderListCreator.borderList(),
                                                blackList,
                                                overlapSize);

        // create the overlapping vectors for the residual and the
        // solution
        overlappingb_ = new OverlappingVector(overlappingMatrix_->overlap());
        overlappingx_ = new OverlappingVector(*overlappingb_);

#if HAVE_MPI
        // create and initialize DUNE's OwnerOverlapCopyCommunication
        // using the domestic overlap
        istlComm_ = new OwnerOverlapCopyCommunication(MPI_COMM_WORLD);
        setupAmgIndexSet(overlappingMatrix_->overlap(), istlComm_->indexSet());
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
    void setupAmgIndexSet(const Overlap &overlap,
                          ParallelIndexSet &istlIndices)
    {
        typedef Dune::OwnerOverlapCopyAttributeSet GridAttributes;
        typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet GridAttributeSet;

        // create DUNE's ParallelIndexSet from a domestic overlap
        istlIndices.beginResize();
        for (int curIdx = 0; curIdx < overlap.numDomestic(); ++curIdx) {
            GridAttributeSet gridFlag =
                overlap.iAmMasterOf(curIdx)
                ? GridAttributes::owner
                : GridAttributes::copy;

            // an index is used by other processes if it is in the
            // domestic or in the foreign overlap.
            bool isShared = overlap.isInOverlap(curIdx);

            assert(curIdx == int(overlap.globalToDomestic(overlap.domesticToGlobal(curIdx))));
            istlIndices.add(/*globalIdx=*/overlap.domesticToGlobal(curIdx),
                            Dune::ParallelLocalIndex<GridAttributeSet>(curIdx, gridFlag, isShared));

        }
        istlIndices.endResize();
    }
#endif

    void setupAmg_(FineOperator &fineOperator)
    {
        if (amg_)
            return;

        int verbosity = 0;
        if (problem_.gridView().comm().rank() == 0)
            verbosity = GET_PARAM(TypeTag, int, LinearSolverVerbosity);

        int rank = problem_.gridView().comm().rank();
        if (verbosity > 1 && rank == 0)
            std::cout << "Setting up the AMG preconditioner\n";

        typedef typename Dune::Amg::SmootherTraits<ParallelSmoother>::Arguments SmootherArgs;

        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        smootherArgs.relaxationFactor = 1.0;

        // specify the coarsen criterion
        //typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<Matrix,Dune::Amg::FirstDiagonal> >
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<Matrix,Dune::Amg::FrobeniusNorm> >
            CoarsenCriterion;
        int coarsenTarget = GET_PARAM(TypeTag, int, AmgCoarsenTarget);
        CoarsenCriterion coarsenCriterion(/*maxLevel=*/15, coarsenTarget);
        coarsenCriterion.setDefaultValuesAnisotropic(GridView::dimension, /*aggregateSizePerDim=*/3);
        if (verbosity > 0)
            coarsenCriterion.setDebugLevel(1);
        else
            coarsenCriterion.setDebugLevel(0); // make the AMG shut up
        coarsenCriterion.setMinCoarsenRate(1.05); // reduce the minium coarsen rate (default is 1.2)
        //coarsenCriterion.setAccumulate(Dune::Amg::noAccu);
        coarsenCriterion.setAccumulate(Dune::Amg::atOnceAccu);
        coarsenCriterion.setSkipIsolated(false);

        // instantiate the AMG preconditioner
#if HAVE_MPI
        amg_ = new AMG(fineOperator,
                       coarsenCriterion,
                       smootherArgs,
                       *istlComm_);
#else
        amg_ = new AMG(fineOperator,
                       coarsenCriterion,
                       smootherArgs);
#endif
    };

    const Problem &problem_;

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
