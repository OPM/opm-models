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
 * \copydoc Ewoms::Linear::ParallelIstlSolverBackend
 */
#ifndef EWOMS_PARALLEL_ISTL_BACKEND_HH
#define EWOMS_PARALLEL_ISTL_BACKEND_HH

#include <ewoms/linear/overlappingbcrsmatrix.hh>
#include <ewoms/linear/overlappingblockvector.hh>
#include <ewoms/linear/overlappingpreconditioner.hh>
#include <ewoms/linear/overlappingscalarproduct.hh>
#include <ewoms/linear/overlappingoperator.hh>
#include <ewoms/linear/parallelbicgstabbackend.hh>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/common/shared_ptr.hh>
#include <dune/common/fvector.hh>

#include <sstream>
#include <iostream>

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(ParallelIstlLinearSolver);

// forward declaration of the required property tags
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(JacobianMatrix);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(GridView);

NEW_PROP_TAG(BorderListCreator);
NEW_PROP_TAG(Overlap);
NEW_PROP_TAG(OverlappingVector);
NEW_PROP_TAG(OverlappingMatrix);
NEW_PROP_TAG(OverlappingScalarProduct);
NEW_PROP_TAG(OverlappingLinearOperator);

//! The type of the linear solver to be used
NEW_PROP_TAG(LinearSolverBackend);
NEW_PROP_TAG(LinearSolverWrapper);
NEW_PROP_TAG(PreconditionerWrapper);


//! The floating point type used internally by the linear solver
NEW_PROP_TAG(LinearSolverScalar);

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
} // namespace Properties
} // namespace Ewoms

namespace Ewoms {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides all unmodified linear solvers from dune-istl
 *
 * To set the linear solver, use
 * \code
 * SET_TYPE_PROP(YourTypeTag, LinearSolverWrapper,
 *               Ewoms::Linear::SolverWrapper$SOLVER<TypeTag>);
 * \endcode
 *
 * The possible choices for '\c $SOLVER' are:
 * - \c Richardson: A fixpoint solver using the Richardson iteration
 * - \c SteepestDescent: The steepest descent solver
 * - \c ConjugatedGradients: A conjugated gradients solver
 * - \c BiCGStab: A stabilized bi-conjugated gradients solver
 * - \c MinRes: A solver based on the  minimized residual algorithm
 * - \c RestartedGMRes: A restarted GMRES solver
 *
 * Chosing the preconditioner works in an analogous way:
 * \code
 * SET_TYPE_PROP(YourTypeTag, PreconditionerWrapper,
 *               Ewoms::Linear::PreconditionerWrapper$PRECONDITIONER<TypeTag>);
 * \endcode
 *
 * Where the choices possible for '\c $PRECONDITIONER' are:
 * - \c Jacobi: A Jacobi preconditioner
 * - \c GaussSeidel: A Gauss-Seidel preconditioner
 * - \c SSOR: A symmetric successive overrelaxation (SSOR) preconditioner
 * - \c SOR: A successive overrelaxation (SOR) preconditioner
 * - \c ILUn: An ILU(n) preconditioner
 * - \c ILU0: An ILU(0) preconditioner. The results of this
 *            preconditioner are the same as setting the
 *            PreconditionerOrder property to 0 and using the ILU(n)
 *            preconditioner. The reason for the existence of ILU0 is
 *            that it is computationally cheaper because it does not
 *            need to consider things which are only required for
 *            higher orders
 * - \c Solver: A BiCGSTAB solver wrapped into the preconditioner
 *              interface (may be useful for parallel computations)
 */
template <class TypeTag>
class ParallelIstlSolverBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, BorderListCreator) BorderListCreator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PreconditionerWrapper) PreconditionerWrapper;
    typedef typename PreconditionerWrapper::SequentialPreconditioner SequentialPreconditioner;

    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverWrapper) LinearSolverWrapper;

    typedef Ewoms::Linear::OverlappingPreconditioner<SequentialPreconditioner,
                                                     Overlap> ParallelPreconditioner;
    typedef Ewoms::Linear::OverlappingScalarProduct<OverlappingVector,
                                                    Overlap> ParallelScalarProduct;
    typedef Ewoms::Linear::OverlappingOperator<OverlappingMatrix,
                                               OverlappingVector,
                                               OverlappingVector> ParallelOperator;

    enum { dimWorld = GridView::dimensionworld };

public:
    ParallelIstlSolverBackend(const Simulator& simulator)
        : simulator_(simulator)
        , gridSequenceNumber_( -1 )
    {
        overlappingMatrix_ = nullptr;
        overlappingb_ = nullptr;
        overlappingx_ = nullptr;
    }

    ~ParallelIstlSolverBackend()
    { cleanup_(); }

    /*!
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverTolerance,
                             "The maximum allowed error between of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, LinearSolverOverlapSize,
                             "The size of the algebraic overlap for the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIterations,
                             "The maximum number of iterations of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                             "The verbosity level of the linear solver");

        LinearSolverWrapper::registerParameters();
        PreconditionerWrapper::registerParameters();
    }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     */
    void eraseMatrix()
    { cleanup_(); }

    void prepareMatrix(const Matrix& M)
    {
        // make sure that the overlapping matrix and block vectors
        // have been created
        prepare_(M);

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlappingMatrix_->assignAdd(M);
    }

    void prepareRhs(const Matrix& M, Vector& b)
    {
        // make sure that the overlapping matrix and block vectors
        // have been created
        prepare_(M);

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
    bool solve(Vector& x)
    {
        Scalar oldSingularLimit = Dune::FMatrixPrecision<Scalar>::singular_limit();
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-50);

        (*overlappingx_) = 0.0;

        int preconditionerIsReady = 1;
        try {
            // update sequential preconditioner
            precWrapper_.prepare(*overlappingMatrix_);
        }
        catch (const Dune::Exception& e) {
            std::cout << "Preconditioner threw exception \"" << e.what()
                      << " on rank " << overlappingMatrix_->overlap().myRank()
                      << "\n"  << std::flush;
            preconditionerIsReady = 0;
        }

        // make sure that the preconditioner is also ready on all peer
        // ranks.
        preconditionerIsReady = simulator_.gridView().comm().min(preconditionerIsReady);
        if (!preconditionerIsReady) {
            Dune::FMatrixPrecision<Scalar>::set_singular_limit(oldSingularLimit);
            return false;
        }

        // create the parallel preconditioner
        ParallelPreconditioner parPreCond(precWrapper_.get(),
                                          overlappingMatrix_->overlap());

        // create the parallel scalar product and the parallel operator
        ParallelScalarProduct parScalarProduct(overlappingMatrix_->overlap());
        ParallelOperator parOperator(*overlappingMatrix_);

        // retrieve the linear solver
        auto& solver = solverWrapper_.get(parOperator, parScalarProduct, parPreCond);

        // run the linear solver and have some fun
        Dune::InverseOperatorResult result;
        int solverSucceeded = 1;
        try {
            solver.apply(*overlappingx_, *overlappingb_, result);
            solverSucceeded = simulator_.gridView().comm().min(solverSucceeded);
        }
        catch (const Dune::Exception& ) {
            solverSucceeded = 0;
            solverSucceeded = simulator_.gridView().comm().min(solverSucceeded);
        }

        // free the unneeded memory of the sequential preconditioner
        // and the linear solver
        solverWrapper_.cleanup();
        precWrapper_.cleanup();

        if (!solverSucceeded) {
            Dune::FMatrixPrecision<Scalar>::set_singular_limit(oldSingularLimit);
            return false;
        }

        // copy the result back to the non-overlapping vector
        overlappingx_->assignTo(x);

        // reset the singularity limit to the same value as before the
        // linear solver was invoked.
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(oldSingularLimit);

        // return the result of the solver
        return result.converged;
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void prepare_(const Matrix& M)
    {
        // if grid has changed the sequence number has changed too
        int curSeqNum = simulator_.gridManager().gridSequenceNumber();
        if( gridSequenceNumber_ == curSeqNum && overlappingMatrix_)
            // the grid has not changed since the overlappingMatrix_has been created, so
            // there's noting to do
            return;

        cleanup_();
        gridSequenceNumber_ = curSeqNum;

        BorderListCreator borderListCreator(simulator_.gridView(),
                                            simulator_.model().dofMapper());

        // create the overlapping Jacobian matrix
        unsigned overlapSize = EWOMS_GET_PARAM(TypeTag, unsigned, LinearSolverOverlapSize);
        overlappingMatrix_ = new OverlappingMatrix(M,
                                                   borderListCreator.borderList(),
                                                   borderListCreator.blackList(),
                                                   overlapSize);

        // create the overlapping vectors for the residual and the
        // solution
        overlappingb_ = new OverlappingVector(overlappingMatrix_->overlap());
        overlappingx_ = new OverlappingVector(*overlappingb_);

        // writeOverlapToVTK_();
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
        for (int lookedAtRank = 0;
             lookedAtRank < simulator_.gridView().comm().size(); ++lookedAtRank) {
            std::cout << "writing overlap for rank " << lookedAtRank << "\n"  << std::flush;
            typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > VtkField;
            int n = simulator_.gridView().size(/*codim=*/dimWorld);
            VtkField isInOverlap(n);
            VtkField rankField(n);
            isInOverlap = 0.0;
            rankField = 0.0;
            assert(rankField.two_norm() == 0.0);
            assert(isInOverlap.two_norm() == 0.0);
            auto vIt = simulator_.gridView().template begin</*codim=*/dimWorld>();
            const auto& vEndIt = simulator_.gridView().template end</*codim=*/dimWorld>();
            const auto& overlap = overlappingMatrix_->overlap();
            for (; vIt != vEndIt; ++vIt) {
                int nativeIdx = simulator_.model().vertexMapper().map(*vIt);
                int localIdx = overlap.foreignOverlap().nativeToLocal(nativeIdx);
                if (localIdx < 0)
                    continue;
                rankField[nativeIdx] = simulator_.gridView().comm().rank();
                if (overlap.peerHasIndex(lookedAtRank, localIdx))
                    isInOverlap[nativeIdx] = 1.0;
            }

            typedef Dune::VTKWriter<GridView> VtkWriter;
            VtkWriter writer(simulator_.gridView(), Dune::VTK::conforming);
            writer.addVertexData(isInOverlap, "overlap");
            writer.addVertexData(rankField, "rank");

            std::ostringstream oss;
            oss << "overlap_rank=" << lookedAtRank;
            writer.write(oss.str().c_str(), Dune::VTK::ascii);
        }
    }

    const Simulator& simulator_;
    int gridSequenceNumber_;

    OverlappingMatrix *overlappingMatrix_;
    OverlappingVector *overlappingb_;
    OverlappingVector *overlappingx_;

    PreconditionerWrapper precWrapper_;
    LinearSolverWrapper solverWrapper_;
};

/*!
 * \brief Macro to create a wrapper around an ISTL solver
 */
#define EWOMS_WRAP_ISTL_SOLVER(SOLVER_NAME, ISTL_SOLVER_NAME)                      \
    template <class TypeTag>                                                       \
    class SolverWrapper##SOLVER_NAME                                               \
    {                                                                              \
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;                    \
        typedef typename GET_PROP_TYPE(TypeTag,                                    \
                                       OverlappingMatrix) OverlappingMatrix;       \
        typedef typename GET_PROP_TYPE(TypeTag,                                    \
                                       OverlappingVector) OverlappingVector;       \
                                                                                   \
        typedef ISTL_SOLVER_NAME<OverlappingVector> ParallelSolver;                \
                                                                                   \
    public:                                                                        \
        SolverWrapper##SOLVER_NAME()                                               \
        {}                                                                         \
                                                                                   \
        static void registerParameters()                                           \
        {}                                                                         \
                                                                                   \
        template <class LinearOperator, class ScalarProduct, class Preconditioner> \
        ParallelSolver& get(LinearOperator& parOperator,                           \
                            ScalarProduct& parScalarProduct,                       \
                            Preconditioner& parPreCond)                            \
        {                                                                          \
            Scalar tolerance = EWOMS_GET_PARAM(TypeTag, Scalar,                    \
                                               LinearSolverTolerance);             \
            int maxIter = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations);\
                                                                                   \
            int verbosity = 0;                                                     \
            if (parOperator.overlap().myRank() == 0)                               \
                verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);  \
            solver_ = new ParallelSolver(parOperator, parScalarProduct,            \
                                         parPreCond, tolerance, maxIter,           \
                                         verbosity);                               \
                                                                                   \
            return *solver_;                                                       \
        }                                                                          \
                                                                                   \
        void cleanup()                                                             \
        { delete solver_; }                                                        \
                                                                                   \
    private:                                                                       \
        ParallelSolver *solver_;                                                   \
    };

EWOMS_WRAP_ISTL_SOLVER(Richardson, Dune::LoopSolver)
EWOMS_WRAP_ISTL_SOLVER(SteepestDescent, Dune::GradientSolver)
EWOMS_WRAP_ISTL_SOLVER(ConjugatedGradients, Dune::CGSolver)
EWOMS_WRAP_ISTL_SOLVER(BiCGStab, Dune::BiCGSTABSolver)
EWOMS_WRAP_ISTL_SOLVER(MinRes, Dune::MINRESSolver)

template <class TypeTag>
class SolverWrapperRestartedGMRes
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;

    typedef Dune::RestartedGMResSolver<OverlappingVector> ParallelSolver;

public:
    SolverWrapperRestartedGMRes()
    {}

    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, GMResRestart,
                             "Number of iterations after which the GMRES linear solver is restarted");
    }

    template <class LinearOperator, class ScalarProduct, class Preconditioner>
    ParallelSolver& get(LinearOperator& parOperator,
                        ScalarProduct& parScalarProduct,
                        Preconditioner& parPreCond)
    {
        Scalar tolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        int maxIter = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations);

        int verbosity = 0;
        if (parOperator.overlap().myRank() == 0)
            verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
        int restartAfter = EWOMS_GET_PARAM(TypeTag, int, GMResRestart);
        solver_ = new ParallelSolver(parOperator,
                                     parScalarProduct,
                                     parPreCond,
                                     tolerance,
                                     restartAfter,
                                     maxIter,
                                     verbosity);

        return *solver_;
    }

    void cleanup()
    { delete solver_; }

private:
    ParallelSolver *solver_;
};

#undef EWOMS_WRAP_ISTL_SOLVER
#undef EWOMS_ISTL_SOLVER_TYPDEF
} // namespace Linear
} // namespace Ewoms

namespace Ewoms {
namespace Properties {
//! make the linear solver shut up by default
SET_INT_PROP(ParallelIstlLinearSolver, LinearSolverVerbosity, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(ParallelIstlLinearSolver, PreconditionerRelaxation, 1.0);

//! set the preconditioner order to 0 by default
SET_INT_PROP(ParallelIstlLinearSolver, PreconditionerOrder, 0);

//! set the GMRes restart parameter to 10 by default
SET_INT_PROP(ParallelIstlLinearSolver, GMResRestart, 10);

//! by default use the same kind of floating point values for the linearization and for
//! the linear solve
SET_TYPE_PROP(ParallelIstlLinearSolver,
              LinearSolverScalar,
              typename GET_PROP_TYPE(TypeTag, Scalar));

SET_PROP(ParallelIstlLinearSolver, OverlappingMatrix)
{
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverScalar) LinearSolverScalar;
    typedef Dune::FieldMatrix<LinearSolverScalar, numEq, numEq> MatrixBlock;
    typedef Dune::BCRSMatrix<MatrixBlock> NonOverlappingMatrix;
    typedef Ewoms::Linear::OverlappingBCRSMatrix<NonOverlappingMatrix> type;
};

SET_TYPE_PROP(ParallelIstlLinearSolver,
              Overlap,
              typename GET_PROP_TYPE(TypeTag, OverlappingMatrix)::Overlap);

SET_PROP(ParallelIstlLinearSolver, OverlappingVector)
{
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverScalar) LinearSolverScalar;
    typedef Dune::FieldVector<LinearSolverScalar, numEq> VectorBlock;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Ewoms::Linear::OverlappingBlockVector<VectorBlock, Overlap> type;
};

SET_PROP(ParallelIstlLinearSolver, OverlappingScalarProduct)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Ewoms::Linear::OverlappingScalarProduct<OverlappingVector, Overlap> type;
};

SET_PROP(ParallelIstlLinearSolver, OverlappingLinearOperator)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef Ewoms::Linear::OverlappingOperator<OverlappingMatrix, OverlappingVector,
                                               OverlappingVector> type;
};

SET_TYPE_PROP(ParallelIstlLinearSolver,
              LinearSolverBackend,
              Ewoms::Linear::ParallelIstlSolverBackend<TypeTag>);

SET_TYPE_PROP(ParallelIstlLinearSolver,
              LinearSolverWrapper,
              Ewoms::Linear::SolverWrapperBiCGStab<TypeTag>);

SET_TYPE_PROP(ParallelIstlLinearSolver,
              PreconditionerWrapper,
              Ewoms::Linear::PreconditionerWrapperILU0<TypeTag>);

//! set the default overlap size to 2
SET_INT_PROP(ParallelIstlLinearSolver, LinearSolverOverlapSize, 2);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(ParallelIstlLinearSolver, LinearSolverMaxIterations, 1000);
} // namespace Properties
} // namespace Ewoms

#endif
