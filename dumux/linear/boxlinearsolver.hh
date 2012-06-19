// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \brief Provides the linear solvers for the vertex-centered finite
 *        volume ("box") method.
 */
#ifndef DUMUX_BOX_LINEAR_SOLVER_HH
#define DUMUX_BOX_LINEAR_SOLVER_HH

#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolverproperties.hh>
#include <dumux/linear/vertexborderlistfromgrid.hh>
#include <dumux/linear/overlappingbcrsmatrix.hh>
#include <dumux/linear/overlappingblockvector.hh>
#include <dumux/linear/overlappingpreconditioner.hh>
#include <dumux/linear/overlappingscalarproduct.hh>
#include <dumux/linear/overlappingoperator.hh>
#include <dumux/boxmodels/common/boxproperties.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {
namespace Properties {
// forward declaration of the required property tags
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(JacobianMatrix);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(LinearSolverOverlapSize);
NEW_PROP_TAG(GridView);
}

/*!
 * \ingroup Linear
 * \brief The base class of the linear solvers for the vertex-centered
 *        finite volume ("box") method.
 *
 * This solver's intention is to be used in conjunction with the box
 * method, so it assumes that the vertices are the only DOFs.
 */
template <class TypeTag>
class BoxLinearSolver
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dumux::Linear::OverlappingScalarProduct<OverlappingVector, Overlap> OverlappingScalarProduct;
    typedef Dumux::Linear::OverlappingOperator<OverlappingMatrix, OverlappingVector, OverlappingVector> OverlappingOperator;

    enum { dimWorld = GridView::dimensionworld };

public:
    BoxLinearSolver(const Problem &problem)
        : problem_(problem)
    {
        overlapMatrix_ = 0;
        overlapb_ = 0;
        overlapx_ = 0;
        overlapSize_ = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, OverlapSize);
        assert(overlapSize_ >= 0);
    }

    ~BoxLinearSolver()
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
    template <class PrecBackend, class SolverBackend>
    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        int verbosity = 0;
        if (problem_.gridView().comm().rank() == 0)
            verbosity = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);

        int maxIter = GET_PARAM_FROM_GROUP(TypeTag, Scalar, LinearSolver, MaxIterations);
        Scalar tolerance = GET_PARAM_FROM_GROUP(TypeTag, Scalar, LinearSolver, Tolerance);

        if (!overlapMatrix_) {
            // make sure that the overlapping matrix and block vectors
            // have been created
            prepare_(M);
        };

        // copy the values of the non-overlapping linear system of
        // equations to the overlapping one. On ther border, we add up
        // the values of all processes (using the assignAdd() methods)
        overlapMatrix_->assignAdd(M);
        overlapb_->assignAddBorder(b);

        //overlapMatrix_->resetFront();
        //overlapb_->resetFront();
        (*overlapx_) = 0.0;
        
        // create sequential and overlapping preconditioners
        PrecBackend seqPreCond(*overlapMatrix_);
        typedef typename PrecBackend::Implementation SeqPreconditioner;
        typedef Dumux::OverlappingPreconditioner<SeqPreconditioner, Overlap> OverlappingPreconditioner;
        OverlappingPreconditioner preCond(seqPreCond.imp(), overlapMatrix_->overlap());

        // create the scalar products and linear operators for ISTL
        OverlappingScalarProduct scalarProd(overlapMatrix_->overlap());
        OverlappingOperator opA(*overlapMatrix_);

        // create the actual solver
        SolverBackend solver(opA,
                             scalarProd,
                             preCond,
                             tolerance,
                             maxIter,
                             verbosity);
#if HAVE_ISTL_FIXPOINT_CRITERION
        OverlappingVector weightVec(*overlapx_);

        // set the default weight of the row to 0 (-> do not consider
        // the row when calculating the error)
        weightVec = 0.0;
        // for rows local to the current peer, ping the model for their
        // relative weight
        for (int i = 0; i < x.size(); ++i) {
            for (int j = 0; j < OverlappingVector::block_type::dimension; ++j) {
                weightVec[i][j] = this->problem_.model().primaryVarWeight(i, j);
            }
        }
        solver.imp().convergenceCriterion().setWeight(weightVec);
#endif
        
        // run the solver
        Dune::InverseOperatorResult result;
        solver.imp().apply(*overlapx_, *overlapb_, result);

        // copy the result back to the non-overlapping vector
        overlapx_->assignTo(x);

        // return the result of the solver
        return result.converged;
    }

protected:
    const Problem &problem_;

private:
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
        overlapMatrix_ = new OverlappingMatrix (M,
                                                borderListCreator.borderList(),
                                                blackList,
                                                overlapSize_);

        // create the overlapping vectors for the residual and the
        // solution
        overlapb_ = new OverlappingVector(overlapMatrix_->overlap());
        overlapx_ = new OverlappingVector(*overlapb_);
    }

    void cleanup_()
    {
        // create the overlapping Jacobian matrix and vectors
        delete overlapMatrix_;
        delete overlapb_;
        delete overlapx_;

        overlapMatrix_ = 0;
        overlapb_ = 0;
        overlapx_ = 0;
    }

    int overlapSize_;
    OverlappingMatrix *overlapMatrix_;
    OverlappingVector *overlapb_;
    OverlappingVector *overlapx_;
};

template <class TypeTag, class Imp>
class PrecNoIterBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Imp Implementation;

    template <class Matrix>
    PrecNoIterBackend(Matrix& A)
    // : imp_(A, 1, GET_PARAM(TypeTag, Scalar, PreconditionerRelaxation)) // SSOR
        : imp_(A, GET_PARAM(TypeTag, Scalar, PreconditionerRelaxation)) // ILU
    // : imp_(1.0) // Richardson
    {}

    Imp& imp()
    {
        return imp_;
    }

private:
    Imp imp_;
};

template <class TypeTag, class Imp>
class PrecIterBackend
{
public:
    typedef Imp Implementation;

    template <class Matrix>
    PrecIterBackend(Matrix& A)
        : imp_(A,
               GET_PARAM(TypeTag, int, PreconditionerIterations),
               GET_PARAM(TypeTag, double, PreconditionerRelaxation))
    {}

    Imp& imp()
    {
        return imp_;
    }

private:
    Imp imp_;
};

/*!
 * \brief A standard solver backend.
 */
template <class TypeTag, class Imp>
class StandardSolverBackend
{
public:
    typedef Imp Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    template <class Operator, class ScalarProduct, class Prec>
    StandardSolverBackend(Operator& A, ScalarProduct& sp, Prec& prec,
                          Scalar tolerance, int maxIter, int verbosity)
        : imp_(A, sp, prec, tolerance, maxIter, verbosity)
    {}

    Imp& imp()
    { return imp_; }

private:
    Imp imp_;
};

/*!
 * \brief Backend for an ILU0-preconditioned BiCGSTAB solver.
 */
template <class TypeTag>
class BoxBiCGStabILU0Solver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    //typedef Dune::SeqSSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    //typedef Dune::Richardson<OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef Dune::SeqILU0<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecNoIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
#if HAVE_ISTL_FIXPOINT_CRITERION
    //typedef Dune::ResidReductionCriterion<OverlappingVector> ConvergenceCrit;
    typedef Dune::FixPointCriterion<OverlappingVector> ConvergenceCrit;
    typedef Dune::BiCGSTABSolver<OverlappingVector, ConvergenceCrit> Solver;
#else
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
#endif
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxBiCGStabILU0Solver(const Problem &problem)
        : ParentType(problem)
    { }

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    { return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b); }
};

/*!
 * \brief Backend for an SOR-preconditioned BiCGSTAB solver.
 */
template <class TypeTag>
class BoxBiCGStabSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxBiCGStabSORSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for an SSOR-preconditioned BiCGSTAB solver.
 */
template <class TypeTag>
class BoxBiCGStabSSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxBiCGStabSSORSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for a Jacobi-preconditioned BiCGSTAB solver.
 */
template <class TypeTag>
class BoxBiCGStabJacSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqJac<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxBiCGStabJacSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for a Gauss-Seidel-preconditioned BiCGSTAB solver.
 */
template <class TypeTag>
class BoxBiCGStabGSSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqGS<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::BiCGSTABSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxBiCGStabGSSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for an ILU0-preconditioned CG solver.
 */
template <class TypeTag>
class BoxCGILU0Solver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqILU0<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecNoIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
#if HAVE_ISTL_FIXPOINT_CRITERION
    //typedef Dune::ResidReductionCriterion<OverlappingVector> ConvergenceCrit;
    typedef Dune::FixPointCriterion<OverlappingVector> ConvergenceCrit;
    typedef Dune::CGSolver<OverlappingVector, ConvergenceCrit> Solver;
#else
    typedef Dune::CGSolver<OverlappingVector> Solver;
#endif
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxCGILU0Solver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for an SOR-preconditioned CG solver.
 */
template <class TypeTag>
class BoxCGSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxCGSORSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for an SSOR-preconditioned CG solver.
 */
template <class TypeTag>
class BoxCGSSORSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqSSOR<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxCGSSORSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for a Jacobi-preconditioned CG solver.
 */
template <class TypeTag>
class BoxCGJacSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqJac<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxCGJacSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

/*!
 * \brief Backend for a Gauss-Seidel-preconditioned CG solver.
 */
template <class TypeTag>
class BoxCGGSSolver : public BoxLinearSolver<TypeTag>
{
    typedef BoxLinearSolver<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef Dumux::Linear::OverlappingBCRSMatrix<Matrix> OverlappingMatrix;
    typedef typename OverlappingMatrix::Overlap Overlap;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef Dumux::Linear::OverlappingBlockVector<typename Vector::block_type, Overlap> OverlappingVector;
    typedef Dune::SeqGS<OverlappingMatrix, OverlappingVector, OverlappingVector> SeqPreconditioner;
    typedef PrecIterBackend<TypeTag, SeqPreconditioner> PrecBackend;
    typedef Dune::CGSolver<OverlappingVector> Solver;
    typedef StandardSolverBackend<TypeTag, Solver> SolverBackend;

public:
    BoxCGGSSolver(const Problem &problem)
        : ParentType(problem)
    {}

    bool solve(const Matrix &M,
               Vector &x,
               const Vector &b)
    {
        return ParentType::template solve<PrecBackend, SolverBackend>(M, x, b);
    }
};

} // namespace Dumux

#endif
