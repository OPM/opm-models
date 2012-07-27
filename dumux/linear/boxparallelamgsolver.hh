// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief Provides a parallel algebraic multi-grid (AMG) linear solver
 *        from DUNE-ISTL for the vertex-centered finite volume ("box")
 *        method.
 */
#ifndef DUMUX_BOX_PARALLEL_AMG_SOLVER_HH
#define DUMUX_BOX_PARALLEL_AMG_SOLVER_HH

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/vertexborderlistfromgrid.hh>
#include <dumux/linear/overlappingbcrsmatrix.hh>
#include <dumux/linear/overlappingblockvector.hh>
#include <dumux/linear/overlappingpreconditioner.hh>
#include <dumux/linear/overlappingscalarproduct.hh>
#include <dumux/linear/overlappingoperator.hh>
#include <dumux/linear/solverpreconditioner.hh>
#include <dumux/linear/boxlinearsolver.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

#if HAVE_PARMETIS
// tell the AMG to partition in parallel
#define PARALLEL_PARTITION 1
#endif

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dune/common/parallel/indexset.hh>
#include <dune/common/mpicollectivecommunication.hh>

namespace Dumux {
namespace Properties {
NEW_PROP_TAG(AMGCoarsenTarget);

//! The target number of DOFs per processor for the parallel algebraic
//! multi-grid solver
SET_INT_PROP(BoxLinearSolverTypeTag, AMGCoarsenTarget, 10000);
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

    typedef Dune::OwnerOverlapCopyCommunication<Dumux::Linear::Index> OwnerOverlapCopyCommunication;

public:
    BoxParallelAmgSolver(const Problem &problem)
        : problem_(problem)
    {
        overlappingMatrix_ = 0;
        overlappingb_ = 0;
        overlappingx_ = 0;
    }

    ~BoxParallelAmgSolver()
    { cleanup_(); }

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
        overlappingb_->assignAddBorder(b);

        (*overlappingx_) = 0.0;

        int coarsenTarget = GET_PARAM(TypeTag, int, AMGCoarsenTarget);

        /////////////
        // set-up the AMG preconditioner
        /////////////
        typedef Dune::OverlappingSchwarzOperator<Matrix,
                                                 Vector,
                                                 Vector,
                                                 OwnerOverlapCopyCommunication> FineOperator;
        FineOperator fineOperator(*overlappingMatrix_, *istlComm_);

        // create the parallel scalar product and the parallel operator
        typedef Dune::OverlappingSchwarzScalarProduct<Vector,OwnerOverlapCopyCommunication> FineScalarProduct;
        FineScalarProduct scalarProduct(*istlComm_);

        // define the smoother used for the AMG and specify its
        // arguments
        typedef Dune::SeqSOR<Matrix,Vector,Vector> SequentialSmoother;
        //typedef Dune::SeqSSOR<Matrix,Vector,Vector> SequentialSmoother;
        //typedef Dune::SeqJac<Matrix,Vector,Vector> SequentialSmoother;
        //typedef Dune::SeqILU0<Matrix,Vector,Vector> SequentialSmoother;
        //typedef Dune::SeqILUn<Matrix,Vector,Vector> SequentialSmoother;

        typedef Dune::BlockPreconditioner<Vector,
                                          Vector,
                                          OwnerOverlapCopyCommunication,
                                          SequentialSmoother> ParallelSmoother;
        typedef typename Dune::Amg::SmootherTraits<ParallelSmoother>::Arguments SmootherArgs;

        SmootherArgs smootherArgs;
        smootherArgs.iterations = 1;
        /*smootherArgs.relaxationFactor = 1.0;*/

        // specify the coarsen criterion
        typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<Matrix,Dune::Amg::FirstDiagonal> >
            CoarsenCriterion;
        CoarsenCriterion coarsenCriterion(/*maxLevel=*/15, coarsenTarget);
        coarsenCriterion.setDebugLevel(0); // make the AMG shut up
        coarsenCriterion.setMinCoarsenRate(1.05); // reduce the minium coarsen rate (default is 1.2)
        coarsenCriterion.setDefaultValuesIsotropic(2);
    
        // instantiate the AMG preconditioner
        typedef Dune::Amg::AMG<FineOperator, Vector, ParallelSmoother, OwnerOverlapCopyCommunication> AMG;
        AMG amg(fineOperator, 
                coarsenCriterion, 
                smootherArgs,
                *istlComm_);
        /////////////
        // set-up the linear solver
        /////////////
        int verbosity = 0;
        if (overlappingMatrix_->overlap().myRank() == 0)
            verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);

        typedef Dune::BiCGSTABSolver<Vector>  SolverType;
        Scalar linearSolverTolerance = GET_PARAM_FROM_GROUP(TypeTag, Scalar, LinearSolver, Tolerance);
        SolverType solver(fineOperator, 
                          scalarProduct,
                          /*preconditioner=*/amg,
                          linearSolverTolerance,
                          /*maxSteps=*/300, 
                          verbosity);

#if HAVE_ISTL_CONVERGENCE_CRITERIA
        // use the fixpoint convergence criterion if possible
        Vector weightVec(*overlappingx_);

        // set the default weight of the row to 0 (-> do not consider
        // the row when calculating the error)
        weightVec = 0.0;
        // for rows local to the current peer, ping the model for their
        // relative weight
        for (unsigned i = 0; i < x.size(); ++i) {
            for (int j = 0; j < Vector::block_type::dimension; ++j) {
                weightVec[i][j] = this->problem_.model().primaryVarWeight(i, j);
            }
        }

        auto convCrit = new Dune::FixPointCriterion<Vector>(problem_.gridView().comm());
        convCrit->setWeight(weightVec);
        convCrit->setTolerance(linearSolverTolerance);

        typedef Dune::ConvergenceCriterion<Vector> ConvergenceCriterion;
        solver.setConvergenceCriterion(std::shared_ptr<ConvergenceCriterion>(convCrit));
#endif

        /////////////
        // run the linear solver
        /////////////
        Dune::InverseOperatorResult solverStatistics;
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
        Linear::VertexBorderListFromGrid<GridView, VertexMapper>
            borderListCreator(problem_.gridView(), problem_.vertexMapper());
        
        // blacklist all entries that belong to ghost and overlap vertices
        std::set<Dumux::Linear::Index> blackList;
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
        int overlapSize = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, OverlapSize);
        overlappingMatrix_ = new OverlappingMatrix (M,
                                                borderListCreator.borderList(),
                                                blackList,
                                                overlapSize);

        // create the overlapping vectors for the residual and the
        // solution
        overlappingb_ = new OverlappingVector(overlappingMatrix_->overlap());
        overlappingx_ = new OverlappingVector(*overlappingb_);

        // create and initialize DUNE's OwnerOverlapCopyCommunication
        // using the domestic overlap
        istlComm_ = new OwnerOverlapCopyCommunication(MPI_COMM_WORLD);
        setupAmgIndexSet(overlappingMatrix_->overlap(), istlComm_->indexSet()); 
        istlComm_->remoteIndices().template rebuild<false>();
    }
 
    void cleanup_()
    {
        // create the overlapping Jacobian matrix and vectors
        delete overlappingMatrix_;
        delete overlappingb_;
        delete overlappingx_;

        overlappingMatrix_ = 0;
        overlappingb_ = 0;
        overlappingx_ = 0;
    }
    
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

            istlIndices.add(overlap.domesticToGlobal(curIdx),
                        Dune::ParallelLocalIndex<GridAttributeSet>(curIdx, gridFlag, isShared));

        }
        istlIndices.endResize();
    }

    const Problem &problem_;
    
    OwnerOverlapCopyCommunication *istlComm_;
    OverlappingMatrix *overlappingMatrix_;
    OverlappingVector *overlappingb_;
    OverlappingVector *overlappingx_;
};

} // namespace Linear
} // namespace Dumux

#endif
