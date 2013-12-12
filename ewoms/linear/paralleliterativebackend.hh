// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2008-2013 by Andreas Lauser
  Copyright (C) 2011-2012 by Bernd Flemisch

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
*/
/*!
 * \file
 * \copydoc Ewoms::Linear::ParallelIterativeSolverBackend
 */
#ifndef EWOMS_PARALLEL_ITERATIVE_BACKEND_HH
#define EWOMS_PARALLEL_ITERATIVE_BACKEND_HH

#include <ewoms/linear/overlappingbcrsmatrix.hh>
#include <ewoms/linear/overlappingblockvector.hh>
#include <ewoms/linear/overlappingpreconditioner.hh>
#include <ewoms/linear/overlappingscalarproduct.hh>
#include <ewoms/linear/overlappingoperator.hh>
#include <ewoms/linear/solverpreconditioner.hh>

#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/istl/solvers.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/preconditioners.hh>

#include <dune/common/shared_ptr.hh>
#include <dune/common/fvector.hh>

#include <sstream>
#include <iostream>

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(ParallelIterativeLinearSolver);

// forward declaration of the required property tags
NEW_PROP_TAG(Problem);
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

/*!
 * \brief The size of the algebraic overlap of the linear solver.
 *
 * Algebraic overlaps can be thought as being the same as the overlap
 * of a grid, but it is only existant for the linear system of
 * equations.
 */
NEW_PROP_TAG(LinearSolverOverlapSize);

/*!
 * \brief Maximum accepted error of the solution of the non-linear solver.
 */
NEW_PROP_TAG(NewtonRelativeTolerance);

/*!
 * \brief Maximum accepted error of the solution of the linear solver.
 */
NEW_PROP_TAG(LinearSolverRelativeTolerance);

/*!
 * \brief Maximum accepted defect of a component for the solution of the linear
 * solver.
 */
NEW_PROP_TAG(LinearSolverAbsoluteTolerance);

/*!
 * \brief Maximum difference between two iterations for the linear solver below
 * which the linear solver consideres itself to be converged.
 */
NEW_PROP_TAG(LinearSolverFixPointTolerance);

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

NEW_PROP_TAG(LinearSolverUseTwoNormReductionCriterion);

//! The order of the sequential preconditioner
NEW_PROP_TAG(PreconditionerOrder);

//! The relaxation factor of the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

//! number of iterations between solver restarts for the GMRES solver
NEW_PROP_TAG(GMResRestart);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Implements a generic linear solver abstraction.
 *
 * This class' intention is to be used in conjunction with the
 * vertex-centered finite volume discretization, so it assumes that
 * the vertices are the only degrees of freedom.  It is also capable
 * of parallel executions on arbitrary grids and is generic in the
 * sense that it allows to combine any linear solver implemented by
 * Dune-ISTL with any preconditioner (except the algebraic multigrid
 * preconditioner). To set the linear solver, use
 * \code SET_TAG_PROP(YourTypeTag, LinearSolver, DesiredLinearSolver); \endcode
 *
 * The possible choices for '\c DesiredLinearSolver' are:
 * - \c SolverWrapperRichardson: A fixpoint solver using the Richardson
 *iteration
 * - \c SolverWrapperSteepestDescent: The steepest descent solver
 * - \c SolverWrapperConjugatedGradients: A conjugated gradients solver
 * - \c SolverWrapperBiCGStab: A stabilized bi-conjugated gradients solver
 * - \c SolverWrapperMinRes: A solver based on the  minimized residual algorithm
 * - \c SolverWrapperGMRes: A restarted GMRES solver
 *
 * Chosing the preconditioner works analogous: \code
 * SET_TAG_PROP(YourTypeTag, Preconditioner, DesiredPreconditioner); \endcode
 *
 * Where the choices possible for '\c DesiredPreconditioner' are:
 * - \c PreconditionerJacobi: A Jacobi preconditioner
 * - \c PreconditionerGaussSeidel: A Gauss-Seidel preconditioner
 * - \c PreconditionerSSOR: A symmetric successive overrelaxation (SSOR)
 *preconditioner
 * - \c PreconditionerSOR: A successive overrelaxation (SOR) preconditioner
 * - \c PreconditionerILU: An ILU(n) preconditioner
 */
template <class TypeTag>
class ParallelIterativeSolverBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, BorderListCreator) BorderListCreator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;

    typedef typename GET_PROP_TYPE(TypeTag,
                                   PreconditionerWrapper) PreconditionerWrapper;
    typedef typename PreconditionerWrapper::SequentialPreconditioner
    SequentialPreconditioner;

    typedef typename GET_PROP_TYPE(TypeTag,
                                   LinearSolverWrapper) LinearSolverWrapper;

    typedef Ewoms::Linear::OverlappingPreconditioner<SequentialPreconditioner,
                                                     Overlap> ParallelPreconditioner;
    typedef Ewoms::Linear::OverlappingScalarProduct<OverlappingVector, Overlap>
    ParallelScalarProduct;
    typedef Ewoms::Linear::OverlappingOperator<OverlappingMatrix, OverlappingVector,
                                               OverlappingVector> ParallelOperator;

    enum { dimWorld = GridView::dimensionworld };

public:
    ParallelIterativeSolverBackend(const Problem &problem) : problem_(problem)
    {
        overlappingMatrix_ = 0;
        overlappingb_ = 0;
        overlappingx_ = 0;
    }

    ~ParallelIterativeSolverBackend()
    { cleanup_(); }

    /*!
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverRelativeTolerance, "The maximum allowed weighted difference between two iterations of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverOverlapSize, "The size of the algebraic overlap for the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIterations, "The maximum number of iterations of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity, "The verbosity level of the linear solver");

        LinearSolverWrapper::registerParameters();
        PreconditionerWrapper::registerParameters();
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
        Scalar oldSingularLimit
            = Dune::FMatrixPrecision<Scalar>::singular_limit();
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-50);

        if (!overlappingMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        };

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlappingMatrix_->assignAdd(M);
        // overlappingMatrix_->resetFront();
        overlappingb_->assignAddBorder(b);
        // overlappingb_->resetFront();

        (*overlappingx_) = 0.0;

        int preconditionerIsReady = 1;
        try
        {
            // update sequential preconditioner
            precWrapper_.prepare(*overlappingMatrix_);
        }
        catch (const Dune::Exception &e)
        {
            std::cout << "Preconditioner threw exception \"" << e.what()
                      << " on rank " << overlappingMatrix_->overlap().myRank()
                      << "\n";
            preconditionerIsReady = 0;
        }

        // make sure that the preconditioner is also ready on all peer
        // ranks.
        preconditionerIsReady
            = problem_.gridView().comm().min(preconditionerIsReady);
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

        // run the linear solver and have some fun
        auto &solver = solverWrapper_.get(parOperator, parScalarProduct, parPreCond);
        Dune::InverseOperatorResult result;
        int solverSucceeded = 1;
        try
        {
            solver.apply(*overlappingx_, *overlappingb_, result);
            solverSucceeded = problem_.gridView().comm().min(solverSucceeded);
        }
        catch (const Dune::Exception &e)
        {
            solverSucceeded = 0;
            solverSucceeded = problem_.gridView().comm().min(solverSucceeded);
        }

        // free the unneeded memory of the sequential preconditioner
        // and the linear solver
        solverWrapper_.cleanup();
        precWrapper_.cleanup();

        if (!solverSucceeded)
            return false;

        // copy the result back to the non-overlapping vector
        overlappingx_->assignTo(x);

        // reset the singularity limit to the same value as before the
        // linear solver was invoked.
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(oldSingularLimit);

        // return the result of the solver
        return result.converged;
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void prepare_(const Matrix &M)
    {
        BorderListCreator borderListCreator(problem_.gridView(),
                                            problem_.model().dofMapper());

        // blacklist all ghost and overlap entries
        std::set<Ewoms::Linear::Index> blackList;
        auto vIt = problem_.gridView().template begin<dimWorld>();
        const auto &vEndIt = problem_.gridView().template end<dimWorld>();
        for (; vIt != vEndIt; ++vIt) {
            if (vIt->partitionType() != Dune::InteriorEntity
                && vIt->partitionType() != Dune::BorderEntity) {
                // we blacklist everything except interior and border
                // vertices
                blackList.insert(problem_.vertexMapper().map(*vIt));
            }
        }

        // create the overlapping Jacobian matrix
        int overlapSize = EWOMS_GET_PARAM(TypeTag, int, LinearSolverOverlapSize);
        overlappingMatrix_
            = new OverlappingMatrix(M, borderListCreator.borderList(),
                                    blackList, overlapSize);

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
             lookedAtRank < problem_.gridView().comm().size(); ++lookedAtRank) {
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
            const auto &vEndIt
                = problem_.gridView().template end</*codim=*/dimWorld>();
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
        ParallelSolver &get(LinearOperator &parOperator,                           \
                            ScalarProduct &parScalarProduct,                       \
                            Preconditioner &parPreCond)                            \
        {                                                                          \
            Scalar tolerance = EWOMS_GET_PARAM(TypeTag, Scalar,                    \
                                               LinearSolverRelativeTolerance);     \
            int maxIter                                                            \
                = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations);        \
                                                                                   \
            int verbosity = 0;                                                     \
            if (parOperator.overlap().myRank() == 0)                               \
                verbosity                                                          \
                    = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);        \
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

EWOMS_WRAP_ISTL_SOLVER(Richardson, Ewoms::LoopSolver)
// EWOMS_WRAP_ISTL_SOLVER(SteepestDescent, Ewoms::GradientSolver)
EWOMS_WRAP_ISTL_SOLVER(ConjugatedGradients, Ewoms::CGSolver)
EWOMS_WRAP_ISTL_SOLVER(BiCGStab, Ewoms::BiCGSTABSolver)
// EWOMS_WRAP_ISTL_SOLVER(MinRes, Ewoms::MINRESSolver)
// EWOMS_WRAP_ISTL_SOLVER(RestartedGMRes, Ewoms::RestartedGMResSolver)

#undef EWOMS_WRAP_ISTL_SOLVER
#undef EWOMS_ISTL_SOLVER_TYPDEF

#define EWOMS_WRAP_ISTL_PRECONDITIONER(PREC_NAME, ISTL_PREC_TYPE)               \
    template <class TypeTag>                                                    \
    class PreconditionerWrapper##PREC_NAME                                      \
    {                                                                           \
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;                 \
        typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix; \
        typedef typename GET_PROP_TYPE(TypeTag,                                 \
                                       OverlappingVector) OverlappingVector;    \
                                                                                \
    public:                                                                     \
        typedef ISTL_PREC_TYPE<JacobianMatrix, OverlappingVector,               \
                               OverlappingVector> SequentialPreconditioner;     \
        PreconditionerWrapper##PREC_NAME()                                      \
        {}                                                                      \
                                                                                \
        static void registerParameters()                                        \
        {                                                                       \
            EWOMS_REGISTER_PARAM(TypeTag, int, PreconditionerOrder,             \
                                 "The order of the preconditioner");            \
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, PreconditionerRelaxation,     \
                                 "The relaxation factor of the "                \
                                 "preconditioner");                             \
        }                                                                       \
                                                                                \
        void prepare(JacobianMatrix &matrix)                                    \
        {                                                                       \
            int order = EWOMS_GET_PARAM(TypeTag, int, PreconditionerOrder);     \
            Scalar relaxationFactor                                             \
                = EWOMS_GET_PARAM(TypeTag, Scalar, PreconditionerRelaxation);   \
            seqPreCond_ = new SequentialPreconditioner(matrix, order,           \
                                                       relaxationFactor);       \
        }                                                                       \
                                                                                \
        SequentialPreconditioner &get()                                         \
        { return *seqPreCond_; }                                                \
                                                                                \
        void cleanup()                                                          \
        { delete seqPreCond_; }                                                 \
                                                                                \
    private:                                                                    \
        SequentialPreconditioner *seqPreCond_;                                  \
    };

EWOMS_WRAP_ISTL_PRECONDITIONER(Jacobi, Dune::SeqJac)
// EWOMS_WRAP_ISTL_PRECONDITIONER(Richardson, Dune::Richardson)
EWOMS_WRAP_ISTL_PRECONDITIONER(GaussSeidel, Dune::SeqGS)
EWOMS_WRAP_ISTL_PRECONDITIONER(SOR, Dune::SeqSOR)
EWOMS_WRAP_ISTL_PRECONDITIONER(SSOR, Dune::SeqSSOR)
EWOMS_WRAP_ISTL_PRECONDITIONER(ILU, Dune::SeqILUn)
EWOMS_WRAP_ISTL_PRECONDITIONER(Solver, Ewoms::Linear::SolverPreconditioner)

#undef EWOMS_WRAP_ISTL_PRECONDITIONER
} // namespace Linear
} // namespace Ewoms

namespace Opm {
namespace Properties {
//! make the linear solver shut up by default
SET_INT_PROP(ParallelIterativeLinearSolver, LinearSolverVerbosity, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(ParallelIterativeLinearSolver, PreconditionerRelaxation, 1.0);

//! set the preconditioner order to 0 by default
SET_INT_PROP(ParallelIterativeLinearSolver, PreconditionerOrder, 0);

//! set the GMRes restart parameter to 10 by default
SET_INT_PROP(ParallelIterativeLinearSolver, GMResRestart, 10);

SET_TYPE_PROP(ParallelIterativeLinearSolver, OverlappingMatrix,
              Ewoms::Linear::OverlappingBCRSMatrix<typename GET_PROP_TYPE(
                  TypeTag, JacobianMatrix)>);
SET_TYPE_PROP(ParallelIterativeLinearSolver, Overlap,
              typename GET_PROP_TYPE(TypeTag, OverlappingMatrix)::Overlap);
SET_PROP(ParallelIterativeLinearSolver, OverlappingVector)
{
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Ewoms::Linear::OverlappingBlockVector<typename Vector::block_type,
                                                  Overlap> type;
};
SET_PROP(ParallelIterativeLinearSolver, OverlappingScalarProduct)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef typename GET_PROP_TYPE(TypeTag, Overlap) Overlap;
    typedef Ewoms::Linear::OverlappingScalarProduct<OverlappingVector, Overlap> type;
};
SET_PROP(ParallelIterativeLinearSolver, OverlappingLinearOperator)
{
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingMatrix) OverlappingMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, OverlappingVector) OverlappingVector;
    typedef Ewoms::Linear::OverlappingOperator<OverlappingMatrix, OverlappingVector,
                                               OverlappingVector> type;
};

SET_TYPE_PROP(ParallelIterativeLinearSolver, LinearSolverBackend,
              Ewoms::Linear::ParallelIterativeSolverBackend<TypeTag>);
SET_TYPE_PROP(ParallelIterativeLinearSolver, LinearSolverWrapper,
              Ewoms::Linear::SolverWrapperBiCGStab<TypeTag>);
SET_TYPE_PROP(ParallelIterativeLinearSolver, PreconditionerWrapper,
              Ewoms::Linear::PreconditionerWrapperILU<TypeTag>);

SET_BOOL_PROP(ParallelIterativeLinearSolver,
              LinearSolverUseTwoNormReductionCriterion, false);

//! set the default for the accepted absolute tolerance (by default, we use 0 to
// disable considering the absolute tolerance)
SET_SCALAR_PROP(ParallelIterativeLinearSolver, LinearSolverAbsoluteTolerance,
                0.0);

//! set the default for the accepted fix-point tolerance (by default, we use 0
// to disable considering the fix-point tolerance)
SET_SCALAR_PROP(ParallelIterativeLinearSolver, LinearSolverFixPointTolerance,
                0.0);

//! set the default overlap size to 3
SET_SCALAR_PROP(ParallelIterativeLinearSolver, LinearSolverOverlapSize, 2);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(ParallelIterativeLinearSolver, LinearSolverMaxIterations, 250);
} // namespace Properties
} // namespace Opm

#endif
