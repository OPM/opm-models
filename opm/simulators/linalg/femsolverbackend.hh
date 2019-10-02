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
 * \copydoc Opm::Linear::FemSolverBackend
 */
#ifndef OPM_FEM_SOLVER_BACKEND_HH
#define OPM_FEM_SOLVER_BACKEND_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>

#define DISABLE_AMG_DIRECTSOLVER 1
#include <dune/fem/solver/istlinverseoperators.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#if HAVE_PETSC
#include <dune/fem/solver/petscinverseoperators.hh>
#endif

#if HAVE_VIENNACL
#include <dune/fem/solver/viennacl.hh>
#endif

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>
#include <opm/simulators/linalg/istlsolverwrappers.hh>
#include <opm/simulators/linalg/femsparsematrixadapter.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/common/fvector.hh>

#include <sstream>
#include <memory>
#include <iostream>

namespace Opm {
namespace Linear {
template <class TypeTag>
class FemSolverBackend;
}} // namespace Linear, Opm


BEGIN_PROPERTIES

NEW_TYPE_TAG(FemSolverBackend);

SET_TYPE_PROP(FemSolverBackend,
              LinearSolverBackend,
              Opm::Linear::FemSolverBackend<TypeTag>);

NEW_PROP_TAG(DiscreteFunction);
NEW_PROP_TAG(SparseMatrixAdapter);

//NEW_PROP_TAG(LinearSolverTolerance);
NEW_PROP_TAG(LinearSolverMaxIterations);
NEW_PROP_TAG(LinearSolverVerbosity);
NEW_PROP_TAG(LinearSolverMaxError);
NEW_PROP_TAG(LinearSolverOverlapSize);
//! The order of the sequential preconditioner
NEW_PROP_TAG(PreconditionerOrder);

//! The relaxation factor of the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

//! Filename for DUNE-FEM solver parameters
NEW_PROP_TAG(FemSolverParameterFileName);

//! make the linear solver shut up by default
SET_INT_PROP(FemSolverBackend, LinearSolverVerbosity, 0);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(FemSolverBackend, LinearSolverMaxIterations, 1000);

SET_SCALAR_PROP(FemSolverBackend, LinearSolverMaxError, 1e7);

//! set the default overlap size to 2
SET_INT_PROP(FemSolverBackend, LinearSolverOverlapSize, 2);

//! set the preconditioner order to 0 by default
SET_INT_PROP(FemSolverBackend, PreconditionerOrder, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(FemSolverBackend, PreconditionerRelaxation, 1.0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_STRING_PROP(FemSolverBackend, FemSolverParameterFileName, "");

//! make the linear solver shut up by default
//SET_SCALAR_PROP(FemSolverBackend, LinearSolverTolerance, 0.01);

SET_PROP(FemSolverBackend, DiscreteFunction)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, DiscreteFunctionSpace) DiscreteFunctionSpace;

    //typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
public:
    // discrete function storing solution data
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace> type;
};

//! Set the type of a global jacobian matrix for linear solvers that are based on
//! dune-istl.
SET_PROP(FemSolverBackend, SparseMatrixAdapter)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, DiscreteFunctionSpace) DiscreteFunctionSpace;
    //typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    // discrete function storing solution data
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace> DiscreteFunction;

public:
    typedef Opm::Linear::FemSparseRowMatrixAdapter< DiscreteFunction > type;
};


END_PROPERTIES

namespace Opm {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Uses the dune-fem infrastructure to solve the linear system of equations.
 */
template <class TypeTag>
class FemSolverBackend
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) Vector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, DiscreteFunctionSpace) DiscreteFunctionSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    // discrete function storing solution data
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace, PrimaryVariables> DiscreteFunction;

    // discrete function to wrap what is used as Vector in eWoms
    typedef Dune::Fem::ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpace>
        VectorWrapperDiscreteFunction;

    template <int d, class SparseMatBackend>
    struct SolverSelector;

    template <int d, class DF>
    struct SolverSelector< d, FemSparseRowMatrixAdapter< DF > >
    {
        typedef Dune::Fem::KrylovInverseOperator<VectorWrapperDiscreteFunction>  type;
    };

#if HAVE_PETSC
    template <int d, class DF >
    struct SolverSelector<d, FemPetscMatrixAdapter< DF > >
    {
        typedef Dune::Fem::PetscInverseOperator<VectorWrapperDiscreteFunction>  type;
    };
#endif

    template <int d, class DF >
    struct SolverSelector<d, FemISTLMatrixAdapter< DF > >
    {
        typedef Dune::Fem::ISTLBICGSTABOp<VectorWrapperDiscreteFunction >  type;
    };

    // select solver type depending on linear operator type
    typedef typename SolverSelector<0, SparseMatrixAdapter>::type InverseLinearOperator;

    enum { dimWorld = GridView::dimensionworld };

public:
    FemSolverBackend(Simulator& simulator)
        : simulator_(simulator)
        , invOp_()
        , rhs_(nullptr)
        , space_(simulator.vanguard().gridPart())
    {
        std::string paramFileName;// = EWOMS_GET_PARAM(TypeTag, std::string, FemSolverParameterFileName);
        if (paramFileName != "")
            Dune::Fem::Parameter::append(paramFileName);
        else {
            // default parameters
            Dune::Fem::Parameter::append("fem.solver.errormeasure", "residualreduction");
            Dune::Fem::Parameter::append("fem.solver.verbose", "true");
            Dune::Fem::Parameter::append("fem.verboserank", "0");

            // Krylov solver
            Dune::Fem::Parameter::append("fem.solver.method", "bicgstab");

            // Preconditioner
            Dune::Fem::Parameter::append("fem.solver.preconditioning.method", "ilu");
            Dune::Fem::Parameter::append("fem.solver.preconditioning.level", "0");
            Dune::Fem::Parameter::append("fem.solver.preconditioning.relaxation", "0.9");
        }
    }

    ~FemSolverBackend()
    { cleanup_(); }

    /*!
     * \brief Register all run-time parameters for the linear solver.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, LinearSolverTolerance,
                             "The maximum allowed error between of the linear solver");
        /*
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

        EWOMS_REGISTER_PARAM(TypeTag, std::string, FemSolverParameterFileName,
                             "The name of the file which contains the parameters for the DUNE-FEM solvers");

                             */
    }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     */
    void eraseMatrix()
    { cleanup_(); }

    void prepare(const SparseMatrixAdapter& op, Vector& b)
    {
        Scalar linearSolverTolerance = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        Scalar linearSolverAbsTolerance = this->simulator_.model().newtonMethod().tolerance() / 100000.0;
        //Scalar linearSolverTolerance = 1e-3;//EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
        //Scalar linearSolverAbsTolerance = this->simulator_.model().newtonMethod().tolerance() / 100000.0;

        // reset linear solver
        if (!invOp_)
          invOp_.reset(new InverseLinearOperator(linearSolverTolerance, linearSolverAbsTolerance));

        setMatrix(op);
        setResidual(b);

        // not needed
        asImp_().rescale_();
    }

    /*!
     * \brief Assign values to the internal data structure for the residual vector.
     *
     * This method also cares about synchronizing that vector with the peer processes.
     */
    void setResidual(const Vector& b)
    { rhs_ = &b ; }

    /*!
     * \brief Retrieve the synchronized internal residual vector.
     *
     * This only deals with entries which are local to the current process.
     */
    void getResidual(Vector& b) const
    {
        assert(rhs_);
        b = *rhs_;
    }

    /*!
     * \brief Sets the values of the residual's Jacobian matrix.
     *
     * This method also synchronizes the data structure across the processes which are
     * involved in the simulation run.
     */
    void setMatrix(const SparseMatrixAdapter& op)
    {
        invOp_->bind(op);
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
        assert(rhs_);
        VectorWrapperDiscreteFunction B("FSB::rhs", space_, *rhs_);

        // solve with right hand side rhs and store in x
        (*invOp_)(B, X);

        // return the result of the solver
        return invOp_->iterations() < 0 ? false : true;
    }

    /*!
     * \brief Return number of iterations used during last solve.
     */
    size_t iterations () const
    {
        assert(invOp_);
        return invOp_->iterations();
    }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void rescale_()
    { }

    void cleanup_()
    {
        invOp_.reset();
        rhs_ = nullptr;
    }

    const Simulator& simulator_;
    std::unique_ptr<InverseLinearOperator> invOp_;
    const Vector* rhs_;
    DiscreteFunctionSpace space_;
};
}} // namespace Linear, Opm

#endif
