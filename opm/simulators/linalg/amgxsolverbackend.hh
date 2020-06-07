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
 * \copydoc Ewoms::Linear::AmgXSolverBackend
 */
#ifndef EWOMS_AMGX_SOLVER_BACKEND_HH
#define EWOMS_AMGX_SOLVER_BACKEND_HH

#include <ewoms/disc/common/fvbaseproperties.hh>

#if USE_AMGX_SOLVERS

#if !HAVE_PETSC || !HAVE_AMGXSOLVER
#error "PETSc and AmgXSolver is needed for the AMGX solver backend"
#endif

#define DISABLE_AMG_DIRECTSOLVER 1
#include <dune/fem/function/petscdiscretefunction.hh>
#include <ewoms/common/genericguard.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/common/fvector.hh>


#include <ewoms/linear/parallelbicgstabbackend.hh>
#include <ewoms/linear/istlsolverwrappers.hh>

#include <AmgXSolver.hpp>

#include <sstream>
#include <memory>
#include <iostream>

namespace Ewoms {
namespace Linear {
template <class TypeTag>
class AmgXSolverBackend;
}} // namespace Linear, Ewoms


BEGIN_PROPERTIES

NEW_TYPE_TAG(AmgXSolverBackend);

SET_TYPE_PROP(AmgXSolverBackend,
              LinearSolverBackend,
              Ewoms::Linear::AmgXSolverBackend<TypeTag>);

NEW_PROP_TAG(LinearSolverTolerance);
NEW_PROP_TAG(LinearSolverAbsTolerance);
NEW_PROP_TAG(LinearSolverMaxIterations);
NEW_PROP_TAG(LinearSolverVerbosity);
NEW_PROP_TAG(LinearSolverMaxError);
NEW_PROP_TAG(LinearSolverOverlapSize);
//! The order of the sequential preconditioner
NEW_PROP_TAG(PreconditionerOrder);

//! The relaxation factor of the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

//! Solver mode for AMGX solver
NEW_PROP_TAG(AmgxSolverMode);

//! Filename for AMGX solver configuration
NEW_PROP_TAG(AmgxSolverConfigFileName);

//! Filename for DUNE-FEM solver parameters
NEW_PROP_TAG(FemSolverParameterFileName);

//! make the linear solver shut up by default
SET_INT_PROP(AmgXSolverBackend, LinearSolverVerbosity, 0);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(AmgXSolverBackend, LinearSolverMaxIterations, 1000);

SET_SCALAR_PROP(AmgXSolverBackend, LinearSolverMaxError, 1e7);

//! set the default overlap size to 2
SET_INT_PROP(AmgXSolverBackend, LinearSolverOverlapSize, 2);

//! set the preconditioner order to 0 by default
SET_INT_PROP(AmgXSolverBackend, PreconditionerOrder, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(AmgXSolverBackend, PreconditionerRelaxation, 1.0);

//! make the linear solver shut up by default
SET_SCALAR_PROP(AmgXSolverBackend, LinearSolverTolerance, 0.01);

//! make the linear solver shut up by default
SET_SCALAR_PROP(AmgXSolverBackend, LinearSolverAbsTolerance, 0.01);

//! set default solver mode for AMGX solver configuratrion
SET_STRING_PROP(AmgXSolverBackend, AmgxSolverMode, "dDDI");

//! set default filename for AMGX solver configuratrion
SET_STRING_PROP(AmgXSolverBackend, AmgxSolverConfigFileName, "amgxconfig.json");

//! set the preconditioner relaxation parameter to 1.0 by default
SET_STRING_PROP(AmgXSolverBackend, FemSolverParameterFileName, "");

END_PROPERTIES

namespace Ewoms {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Provides a linear solver backend that utilizes the AMG-X package.
 */
template <class TypeTag>
class AmgXSolverBackend
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) LinearOperator;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, DiscreteFunctionSpace) DiscreteFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, DiscreteFunction) DiscreteFunction;

    // discrete function to wrap what is used as Vector in eWoms
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace>
        VectorWrapperDiscreteFunction;
    typedef Dune::Fem::PetscDiscreteFunction<DiscreteFunctionSpace>
        PetscDiscreteFunctionType;

    enum { dimWorld = GridView::dimensionworld };

public:
    AmgXSolverBackend(Simulator& simulator)
        : simulator_(simulator)
        , space_(simulator.vanguard().gridPart())
        , amgxSolver_()
        , rhs_(nullptr)
        , iterations_(0)
    {
        std::string paramFileName = EWOMS_GET_PARAM(TypeTag, std::string, FemSolverParameterFileName);
        if (paramFileName != "")
            Dune::Fem::Parameter::append(paramFileName);
    }

    ~AmgXSolverBackend()
    {
        cleanup_();
        if (A_) {
            ::Dune::Petsc::MatDestroy(A_.operator ->());
            A_.reset();
        }
    }

    /*!
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverTolerance,
                             "The maximum allowed error between of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverAbsTolerance,
                             "The maximum allowed error between of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIterations,
                             "The maximum number of iterations of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                             "The verbosity level of the linear solver");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, LinearSolverOverlapSize,
                             "The size of the algebraic overlap for the linear solver");

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverMaxError,
                             "The maximum residual error which the linear solver tolerates"
                             " without giving up");

        EWOMS_REGISTER_PARAM(TypeTag, int, PreconditionerOrder,
                             "The order of the preconditioner");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, PreconditionerRelaxation,
                             "The relaxation factor of the preconditioner");

        EWOMS_REGISTER_PARAM(TypeTag, std::string, AmgxSolverMode,
                             "The name of the solver mode for AMGX");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, AmgxSolverConfigFileName,
                             "The name of the file which contains the AMGX solver configuration");

        EWOMS_REGISTER_PARAM(TypeTag, std::string, FemSolverParameterFileName,
                             "The name of the file which contains the parameters for the DUNE-FEM solvers");

    }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     */
    void eraseMatrix()
    { cleanup_(); }

    void prepare(const LinearOperator& op, Vector& b)
    {
        //Scalar linearSolverTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        //Scalar linearSolverAbsTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverAbsTolerance);

        if (!amgxSolver_) {
            // reset linear solver
            std::string mode = EWOMS_GET_PARAM(TypeTag, std::string, AmgxSolverMode);
            std::string solverconfig = EWOMS_GET_PARAM(TypeTag, std::string, AmgxSolverConfigFileName);

            amgxSolver_.reset(new AmgXSolver());
            amgxSolver_->initialize(PETSC_COMM_WORLD, mode, solverconfig);

            if (!A_) {
                A_.reset(new Mat());
            }

            // convert MATBAIJ to MATAIJ which is needed by AmgXSolver
            ::Dune::Petsc::ErrorCheck(::MatConvert(op.petscMatrix(), MATAIJ, MAT_INITIAL_MATRIX, A_.operator ->()));
            // set up the matrix used by AmgX
            amgxSolver_->setA(*A_);
        }

        // store pointer to right hand side
        rhs_ = &b;
    }

    /*!
     * \brief Actually solve the linear system of equations.
     *
     * \return true if the residual reduction could be achieved, else false.
     */
    bool solve(Vector& x)
    {
        // wrap x into discrete function X (no copy)
        VectorWrapperDiscreteFunction X("FSB::x", space_, x);
        VectorWrapperDiscreteFunction B("FSB::rhs", space_, *rhs_);

        if (!petscRhs_)
            petscRhs_.reset(new PetscDiscreteFunctionType("AMGX::rhs", space_));

        if (!petscX_)
            petscX_.reset(new PetscDiscreteFunctionType("AMGX::X", space_));

        petscRhs_->assign(B);
        petscX_->assign(X);

        // solve with right hand side rhs and store in x
        Vec& p = *(petscX_->petscVec());
        Vec& rhs = *(petscRhs_->petscVec());

        try {
            amgxSolver_->solve(p, rhs);
        }
        catch (...) {
            OPM_THROW(Opm::NumericalIssue, "AmgXSolver: no convergence of linear solver!");
        }

        // copy result to ewoms solution
        X.assign(*petscX_);

        amgxSolver_->getIters(iterations_);

        // return the result of the solver
        return true;
    }

    /*!
     * \brief Return number of iterations used during last solve.
     */
    size_t iterations () const {
        return iterations_;
    }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void cleanup_()
    {
        if (amgxSolver_) {
            amgxSolver_->finalize();
            amgxSolver_.reset();
        }

        rhs_ = nullptr;

        petscRhs_.reset();
        petscX_.reset();
    }

    const Simulator& simulator_;

    DiscreteFunctionSpace space_;
    std::unique_ptr<Mat> A_;

    std::unique_ptr<PetscDiscreteFunctionType> petscRhs_;
    std::unique_ptr<PetscDiscreteFunctionType> petscX_;

    std::unique_ptr<AmgXSolver> amgxSolver_;

    Vector* rhs_;
    int iterations_;
};
}} // namespace Linear, Ewoms

#endif // HAVE_DUNE_FEM

#endif
