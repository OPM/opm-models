// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
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
 * \copydoc Dumux::Linear::BoxParallelSolver
 */
#ifndef DUMUX_BOX_LINEAR_SOLVER_HH
#define DUMUX_BOX_LINEAR_SOLVER_HH

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

#include <dumux/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/common/shared_ptr.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

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
 */
NEW_PROP_TAG(LinearSolverTolerance);


/*!
 * \brief Maximum accepted defect of a component for the solution of the linear solver.
 */
NEW_PROP_TAG(LinearSolverAbsTolerance);

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
 *
 * \brief Implements a generic linear solver abstraction for the
 *        vertex-centered finite volume ("box") method.
 *
 * This class' intention is to be used in conjunction with the box
 * method, so it assumes that the vertices are the only degrees of
 * freedom.  It is also capable of parallel executions on arbitrary
 * grids and is generic in the sense that it allows to combine any
 * linear solver implemented by Dune-ISTL with any preconditioner
 * (except the algebraic multigrid preconditioner). To set the linear
 * solver, use
 * \code SET_TYPE_PROP(YourTypeTag, LinearSolverWrapper,Dumux::Linear::DesiredLinearSolver<TypeTag>); \endcode
 *
 * The possible choices for '\c DesiredLinearSolver' are:
 * - \c SolverWrapperLoop: A fixpoint solver (using the Richardson iteration)
 * - \c SolverWrapperGradients: The steepest descent solver
 * - \c SolverWrapperCG: A conjugated gradients solver
 * - \c SolverWrapperBiCGStab: A stabilized bi-conjugated gradients solver
 * - \c SolverWrapperMinRes: A solver based on the  minimized residual algorithm
 * - \c SolverWrapperGMRes: A restarted GMRES solver
 *
 * Chosing the preconditioner works analogous: \code
 * SET_TYPE_PROP(YourTypeTag, PreconditionerWrapper,Dumux::Linear::DesiredPreconditioner<TypeTag>); \endcode
 *
 * Where the choices possible for '\c DesiredPreconditioner' are:
 * - \c PreconditionerWrapperJacobi: A Jacobi preconditioner
 * - \c PreconditionerWrapperGaussSeidel: A Gauss-Seidel preconditioner
 * - \c PreconditionerWrapperSSOR: A symmetric successive overrelaxation (SSOR) preconditioner
 * - \c PreconditionerWrapperSOR: A successive overrelaxation (SOR) preconditioner
 * - \c PreconditionerWrapperILU: An ILU(n) preconditioner
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
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, Scalar, LinearSolverTolerance, "The maximum tolerance of the linear solver");
        REGISTER_PARAM(TypeTag, Scalar, LinearSolverAbsTolerance, "The maximum tolerance of the defect of a component for the linear solver");
        REGISTER_PARAM(TypeTag, int, LinearSolverOverlapSize, "The size of the algebraic overlap for the linear solver");
        REGISTER_PARAM(TypeTag, int, LinearSolverMaxIterations, "The maximum number of iterations of the linear solver");
        REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity, "The verbosity level of the linear solver");
        REGISTER_PARAM(TypeTag, int, PreconditionerOrder, "The order of the preconditioner");
        REGISTER_PARAM(TypeTag, Scalar, PreconditionerRelaxation, "The relaxation factor of the preconditioner");
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
        if (!overlappingMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        };

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlappingMatrix_->assignAdd(M);
        //overlappingMatrix_->resetFront();
        overlappingb_->assignAddBorder(b);
        //overlappingb_->resetFront();

        (*overlappingx_) = 0.0;

        int preconditionerIsReady = 1;
        try {
            // run the solver
            // update sequential preconditioner
            precWrapper_.prepare(*overlappingMatrix_);
        }
        catch (const Dune::Exception &e) {
            std::cout << "Preconditioner threw exception \"" << e.what() << " on rank " << overlappingMatrix_->overlap().myRank() << "\n";
            preconditionerIsReady = 0;
        }

        preconditionerIsReady = problem_.gridView().comm().min(preconditionerIsReady);
        if (!preconditionerIsReady)
            return false;

        // create the parallel preconditioner
        ParallelPreconditioner parPreCond(precWrapper_.get(), overlappingMatrix_->overlap());

        // create the parallel scalar product and the parallel operator
        ParallelScalarProduct parScalarProduct(overlappingMatrix_->overlap());
        ParallelOperator parOperator(*overlappingMatrix_);

        // run the linear solver and have some fun
        auto &solver = solverWrapper_.get(parOperator, parScalarProduct, parPreCond);

        // use the fixpoint convergence criterion if possible
        OverlappingVector weightVec(*overlappingx_);

        // set the default weight of the row to 0 (-> do not consider
        // the row when calculating the error)
        weightVec = 0.0;
        // for rows local to the current peer, ping the model for their
        // relative weight
        const auto &foreignOverlap = overlappingMatrix_->overlap().foreignOverlap();
        for (unsigned localIdx = 0; localIdx < unsigned(foreignOverlap.numLocal()); ++localIdx) {
            int nativeIdx = foreignOverlap.localToNative(localIdx);
            for (int eqIdx = 0; eqIdx < Vector::block_type::dimension; ++eqIdx) {
                weightVec[localIdx][eqIdx] = this->problem_.model().eqWeight(nativeIdx, eqIdx);
            }
        }

        // create a fixpoint convergence criterion
        Scalar linearSolverTolerance = GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        Scalar linearSolverAbsTolerance = GET_PARAM(TypeTag, Scalar, LinearSolverAbsTolerance);
        auto *convCrit = new Dumux::WeightedResidReductionCriterion<OverlappingVector,
                                                                    typename GridView::CollectiveCommunication>
            (problem_.gridView().comm(),
             weightVec,
             linearSolverTolerance,
             linearSolverAbsTolerance);

        // tell the linear solver to use it
        typedef Dumux::ConvergenceCriterion<OverlappingVector> ConvergenceCriterion;
        solver.setConvergenceCriterion(Dune::shared_ptr<ConvergenceCriterion>(convCrit));

        Dumux::InverseOperatorResult result;
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
        int overlapSize = GET_PARAM(TypeTag, int, LinearSolverOverlapSize);
        overlappingMatrix_ = new OverlappingMatrix (M,
                                                    borderListCreator.borderList(),
                                                    blackList,
                                                    overlapSize);

        // create the overlapping vectors for the residual and the
        // solution
        overlappingb_ = new OverlappingVector(overlappingMatrix_->overlap());
        overlappingx_ = new OverlappingVector(*overlappingb_);

        //writeOverlapToVTK_();
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

    void writeOverlapToVTK_()
    {
        for (int lookedAtRank = 0; lookedAtRank < problem_.gridView().comm().size(); ++lookedAtRank) {
            std::cout << "writing overlap for rank " << lookedAtRank << "\n";
            typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > VtkField;
            int n = problem_.gridView().size(/*codim=*/dimWorld);
            VtkField isInOverlap(n);
            VtkField rankField(n);
            isInOverlap = 0.0;
            rankField = 0.0;
            assert(rankField.two_norm() == 0.0);
            assert(isInOverlap.two_norm() == 0.0);
            auto vIt = problem_.gridView().template begin</*codim=*/dimWorld>();
            const auto &vEndIt = problem_.gridView().template end</*codim=*/dimWorld>();
            const auto &overlap = overlappingMatrix_->overlap();
            for (; vIt != vEndIt; ++vIt) {
                int nativeIdx = problem_.vertexMapper().map(*vIt);
                int localIdx = overlap.foreignOverlap().nativeToLocal(nativeIdx);
                if (localIdx < 0)
                    continue;
                rankField[nativeIdx] = problem_.gridView().comm().rank();
                if (overlap.peerHasIndex(lookedAtRank, localIdx))
                    isInOverlap[nativeIdx] = 1.0;
            };

            typedef Dune::VTKWriter<GridView> VtkWriter;
            VtkWriter writer(problem_.gridView(), Dune::VTK::conforming);
            writer.addVertexData(isInOverlap, "overlap");
            writer.addVertexData(rankField, "rank");

            std::ostringstream oss;
            oss << "overlap_rank=" << lookedAtRank;
            writer.write(oss.str().c_str(), Dune::VTK::ascii);
        }
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
        Scalar tolerance = GET_PARAM(TypeTag, Scalar, LinearSolverTolerance); \
        int maxIter = GET_PARAM(TypeTag, Scalar, LinearSolverMaxIterations); \
                                                                        \
        int verbosity = 0;                                              \
        if (parOperator.overlap().myRank() == 0)                        \
            verbosity = GET_PARAM(TypeTag, int, LinearSolverVerbosity); \
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

EWOMS_WRAP_ISTL_SOLVER(Loop, Dumux::LoopSolver)
//EWOMS_WRAP_ISTL_SOLVER(Gradient, Dumux::GradientSolver)
EWOMS_WRAP_ISTL_SOLVER(CG, Dumux::CGSolver)
EWOMS_WRAP_ISTL_SOLVER(BiCGStab, Dumux::BiCGSTABSolver)
//EWOMS_WRAP_ISTL_SOLVER(MinRes, Dumux::MINRESSolver)
//EWOMS_WRAP_ISTL_SOLVER(RestartedGMRes, Dumux::RestartedGMResSolver)

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
