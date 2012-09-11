// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
 * \brief Provides the linear solvers for the vertex-centered finite
 *        volume ("box") method.
 */
#ifndef DUMUX_BOX_LINEAR_SOLVER_HH
#define DUMUX_BOX_LINEAR_SOLVER_HH

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/vertexborderlistfromgrid.hh>
#include <dumux/linear/overlappingbcrsmatrix.hh>
#include <dumux/linear/overlappingblockvector.hh>
#include <dumux/linear/overlappingpreconditioner.hh>
#include <dumux/linear/overlappingscalarproduct.hh>
#include <dumux/linear/overlappingoperator.hh>
#include <dumux/linear/solverpreconditioner.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {
namespace Properties {
NEW_TYPE_TAG(BoxLinearSolverTypeTag);

// forward declaration of the required property tags
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(JacobianMatrix);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(GridView);

NEW_PROP_TAG(Overlap);
NEW_PROP_TAG(OverlappingVector);
NEW_PROP_TAG(OverlappingMatrix);
NEW_PROP_TAG(OverlappingScalarProduct);
NEW_PROP_TAG(OverlappingLinearOperator);

//! The type of the linear solver to be used
NEW_PROP_TAG(LinearSolver);
NEW_PROP_TAG(LinearSolverWrapper);
NEW_PROP_TAG(PreconditionerWrapper);

/*!
 * \brief The size of the algebraic overlap of the linear solver.
 *
 * Algebraic overlaps can be thought as being the same as the overlap
 * of a grid, but it is only existant for the linear system of
 * equations.
 */
NEW_PROP_TAG(LinearSolverOverlapSize);

/*!
 * \brief Maximum accepted error of the solution of the linear solver.
 *
 * If fixpoint criterion was patched into DUNE-ISTL, this is the
 * maximum of the weighted difference between two iterations of the
 * linear solver (i.e. the linear solver basically uses the same
 * convergence criterion as the non-linear solver in this case). If
 * the vanilla DUNE-ISTL is used, this is the reduction of the
 * two-norm of the residual of the solution relative to the two-norm
 * of the residual of the initial solution.
 */
NEW_PROP_TAG(LinearSolverTolerance);

/*!
 * \brief Specifies the verbosity of the linear solver
 *
 * By default it is 0, i.e. it doesn't print anything. Setting this
 * property to 1 prints aggregated convergence rates, 2 prints the
 * convergence rate of every iteration of the scheme.
 */
NEW_PROP_TAG(LinearSolverVerbosity);

//! Maximum number of iterations eyecuted by the linear solver
NEW_PROP_TAG(LinearSolverMaxIterations);

//! The order of the sequential preconditioner
NEW_PROP_TAG(PreconditionerOrder);

//! The relaxation factor of the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

//! number of iterations between solver restarts for the GMRES solver
NEW_PROP_TAG(GMResRestart);
}

namespace Linear {
/*!
 * \ingroup Linear
 * \brief The base class of the linear solvers for the vertex-centered
 *        finite volume ("box") method.
 *
 * This solver's intention is to be used in conjunction with the box
 * method, so it assumes that the vertices are the only DOFs.
 */
template <class TypeTag>
class BoxParallelSolver
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

    typedef typename GET_PROP_TYPE(TypeTag, PreconditionerWrapper) PreconditionerWrapper;
    typedef typename PreconditionerWrapper::SequentialPreconditioner SequentialPreconditioner;

    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverWrapper) LinearSolverWrapper;

    typedef Dumux::Linear::OverlappingPreconditioner<SequentialPreconditioner, Overlap> ParallelPreconditioner;
    typedef Dumux::Linear::OverlappingScalarProduct<OverlappingVector, Overlap> ParallelScalarProduct;
    typedef Dumux::Linear::OverlappingOperator<OverlappingMatrix, OverlappingVector, OverlappingVector> ParallelOperator;

    enum { dimWorld = GridView::dimensionworld };

public:
    BoxParallelSolver(const Problem &problem)
        : problem_(problem)
    {
        overlappingMatrix_ = 0;
        overlappingb_ = 0;
        overlappingx_ = 0;
    }

    ~BoxParallelSolver()
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

        // run the solver
        // update sequential preconditioner
        precWrapper_.prepare(*overlappingMatrix_);

        // create the parallel preconditioner
        ParallelPreconditioner parPreCond(precWrapper_.get(), overlappingMatrix_->overlap());

        // create the parallel scalar product and the parallel operator
        ParallelScalarProduct parScalarProduct(overlappingMatrix_->overlap());
        ParallelOperator parOperator(*overlappingMatrix_);

        // run the linear solver and have some fun
        auto &solver = solverWrapper_.get(parOperator, parScalarProduct, parPreCond);

#if HAVE_ISTL_CONVERGENCE_CRITERIA
        // use the fixpoint convergence criterion if possible
        OverlappingVector weightVec(*overlappingx_);

        // set the default weight of the row to 0 (-> do not consider
        // the row when calculating the error)
        weightVec = 0.0;
        // for rows local to the current peer, ping the model for their
        // relative weight
        for (unsigned i = 0; i < x.size(); ++i) {
            for (int j = 0; j < Vector::block_type::dimension; ++j) {
                weightVec[i][j] = this->problem_.model().eqWeight(i, j);
            }
        }

        // create a fixpoint convergence criterion
        Scalar linearSolverTolerance = GET_PARAM_FROM_GROUP(TypeTag, Scalar, LinearSolver, Tolerance);
        auto *convCrit = new Dune::WeightedResidReductionCriterion<OverlappingVector, 
                                                                   typename GridView::CollectiveCommunication>(problem_.gridView().comm(), 
                                                                                                               weightVec,
                                                                                                               linearSolverTolerance);

        // tell the linear solver to use it
        typedef Dune::ConvergenceCriterion<OverlappingVector> ConvergenceCriterion;
        solver.setConvergenceCriterion(std::shared_ptr<ConvergenceCriterion>(convCrit));
#endif

        Dune::InverseOperatorResult result;
        solver.apply(*overlappingx_, *overlappingb_, result);

        // free the unneeded memory of the sequential preconditioner and the linear solver
        solverWrapper_.cleanup();
        precWrapper_.cleanup();

        // copy the result back to the non-overlapping vector
        overlappingx_->assignTo(x);

        // return the result of the solver
        return result.converged;
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
        
        // blacklist all ghost and overlap entries
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

    void writeOverlapToVTK_(int lookedAtRank)
    {
        std::cout << "writing overlap for rank " << lookedAtRank << "\n";
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > VtkField;
        int n = problem_.gridView().size(/*codim=*/dimWorld);
        VtkField isInOverlap(n, 0.0);
        VtkField rankField(n, 0.0);
        auto vIt = problem_.gridView().template begin</*codim=*/dimWorld>();
        const auto &vEndIt = problem_.gridView().template end</*codim=*/dimWorld>();
        const auto &overlap = overlappingMatrix_->overlap();
        for (; vIt != vEndIt; ++vIt) {
            int globalVertexIdx = problem_.vertexMapper().map(*vIt);
            rankField[globalVertexIdx] = problem_.gridView().comm().rank();
            if (overlap.peerHasIndex(lookedAtRank, globalVertexIdx))
                isInOverlap[globalVertexIdx] = 1.0;
        };
        
        typedef Dune::VTKWriter<GridView> VtkWriter;
        VtkWriter writer(problem_.gridView(), Dune::VTK::conforming);
        writer.addVertexData(isInOverlap, "overlap");
        writer.addVertexData(rankField, "rank");
        writer.write("overlap", Dune::VTK::ascii);
    }
    
    const Problem &problem_;

    OverlappingMatrix *overlappingMatrix_;
    OverlappingVector *overlappingb_;
    OverlappingVector *overlappingx_;

    PreconditionerWrapper precWrapper_;
    LinearSolverWrapper solverWrapper_;
};

/*!
 * \brief Macro to create a wrapper around an ISTL solver
 */
#define EWOMS_WRAP_ISTL_SOLVER(SOLVER_NAME, ISTL_SOLVER_TYPE)           \
    template <class TypeTag>                                            \
    class SolverWrapper##SOLVER_NAME                                    \
    {                                                                   \
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;             \
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix; \
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector; \
                                                                        \
    typedef ISTL_SOLVER_TYPE<OverlappingVector> ParallelSolver;         \
                                                                        \
    public:                                                             \
    SolverWrapper##SOLVER_NAME()                                        \
    {}                                                                  \
                                                                        \
    template <class LinearOperator, class ScalarProduct, class Preconditioner> \
    ParallelSolver &get(LinearOperator &parOperator,                    \
                        ScalarProduct &parScalarProduct,                \
                        Preconditioner &parPreCond)                     \
    {                                                                   \
        Scalar tolerance = GET_PARAM_FROM_GROUP(TypeTag, Scalar, LinearSolver, Tolerance); \
        int maxIter = GET_PARAM_FROM_GROUP(TypeTag, Scalar, LinearSolver, MaxIterations); \
                                                                        \
        int verbosity = 0;                                              \
        if (parOperator.overlap().myRank() == 0)                        \
            verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity); \
        solver_ = new ParallelSolver(parOperator,                       \
                                     parScalarProduct,                  \
                                     parPreCond,                        \
                                     tolerance,                         \
                                     maxIter,                           \
                                     verbosity);                        \
                                                                        \
        return *solver_;                                                \
    }                                                                   \
                                                                        \
    void cleanup()                                                      \
    { delete solver_; }                                                 \
                                                                        \
    private:                                                            \
    ParallelSolver *solver_;                                            \
    };

EWOMS_WRAP_ISTL_SOLVER(Loop, Dune::LoopSolver)
//EWOMS_WRAP_ISTL_SOLVER(Gradient, Dune::GradientSolver)
EWOMS_WRAP_ISTL_SOLVER(CG, Dune::CGSolver)
EWOMS_WRAP_ISTL_SOLVER(BiCGStab, Dune::BiCGSTABSolver)
//EWOMS_WRAP_ISTL_SOLVER(MinRes, Dune::MINRESSolver)
//EWOMS_WRAP_ISTL_SOLVER(RestartedGMRes, Dune::RestartedGMResSolver)

#undef EWOMS_WRAP_ISTL_SOLVER
#undef EWOMS_ISTL_SOLVER_TYPDEF

#define EWOMS_WRAP_ISTL_PRECONDITIONER(PREC_NAME, ISTL_PREC_TYPE)       \
    template <class TypeTag>                                            \
    class PreconditionerWrapper##PREC_NAME                              \
    {                                                                   \
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;         \
        typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix; \
        typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector; \
                                                                        \
    public:                                                             \
        typedef ISTL_PREC_TYPE<JacobianMatrix, OverlappingVector, OverlappingVector> SequentialPreconditioner; \
        PreconditionerWrapper##PREC_NAME()                              \
        {}                                                              \
                                                                        \
        void prepare(JacobianMatrix &matrix)                            \
        {                                                               \
            int order = GET_PARAM(TypeTag, int, PreconditionerOrder);   \
            Scalar relaxationFactor = GET_PARAM(TypeTag, Scalar, PreconditionerRelaxation); \
            seqPreCond_ = new SequentialPreconditioner(matrix, order, relaxationFactor); \
        }                                                               \
                                                                        \
        SequentialPreconditioner &get()                                 \
        { return *seqPreCond_; }                                        \
                                                                        \
        void cleanup()                                                  \
        { delete seqPreCond_; }                                         \
                                                                        \
    private:                                                            \
        SequentialPreconditioner *seqPreCond_;                          \
    };

EWOMS_WRAP_ISTL_PRECONDITIONER(Jacobi, Dune::SeqJac)
//EWOMS_WRAP_ISTL_PRECONDITIONER(Richardson, Dune::Richardson)
EWOMS_WRAP_ISTL_PRECONDITIONER(GaussSeidel, Dune::SeqGS)
EWOMS_WRAP_ISTL_PRECONDITIONER(SOR, Dune::SeqSOR)
EWOMS_WRAP_ISTL_PRECONDITIONER(SSOR, Dune::SeqSSOR)
EWOMS_WRAP_ISTL_PRECONDITIONER(ILU, Dune::SeqILUn)
EWOMS_WRAP_ISTL_PRECONDITIONER(Solver, Dumux::Linear::SolverPreconditioner)

#undef EWOMS_WRAP_ISTL_PRECONDITIONER
} // namespace Linear

namespace Properties {
//! make the linear solver shut up by default
SET_INT_PROP(BoxLinearSolverTypeTag, LinearSolverVerbosity, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(BoxLinearSolverTypeTag, PreconditionerRelaxation, 1.0);

//! set the preconditioner order to 0 by default
SET_INT_PROP(BoxLinearSolverTypeTag, PreconditionerOrder, 0);

//! set the GMRes restart parameter to 10 by default
SET_INT_PROP(BoxLinearSolverTypeTag, GMResRestart, 10);

SET_TYPE_PROP(BoxLinearSolverTypeTag, 
              OverlappingMatrix,
              Dumux::Linear::OverlappingBCRSMatrix<typename GET_PROP_TYPE(TypeTag, JacobianMatrix)>);
SET_TYPE_PROP(BoxLinearSolverTypeTag, 
              Overlap,
              typename GET_PROP_TYPE(TypeTag, OverlappingMatrix)::Overlap);
SET_PROP(BoxLinearSolverTypeTag, 
         OverlappingVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> type;
};
SET_PROP(BoxLinearSolverTypeTag, 
         OverlappingScalarProduct)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Dumux::Linear::OverlappingScalarProduct<OverlappingVector, Overlap> type;
};
SET_PROP(BoxLinearSolverTypeTag, 
         OverlappingLinearOperator)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef Dumux::Linear::OverlappingOperator<OverlappingMatrix,
                                               OverlappingVector,
                                               OverlappingVector> type;
};
} // namespace Properties
} // namespace Dumux

#endif
