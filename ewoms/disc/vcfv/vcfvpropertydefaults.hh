// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
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
 * \ingroup VcfvModel
 *
 * \brief Defines defaults for the common properties of the VCVF discretizations.
 */
#ifndef EWOMS_VCFV_PROPERTY_DEFAULTS_HH
#define EWOMS_VCFV_PROPERTY_DEFAULTS_HH

#include "vcfvproperties.hh"
#include "vcfvassembler.hh"
#include "vcfvmodel.hh"
#include "vcfvelementgeometry.hh"
#include "vcfvlocalresidual.hh"
#include "vcfvlocaljacobian.hh"
#include "vcfvelementcontext.hh"
#include "vcfvvolumevariables.hh"
#include "vcfvconstraints.hh"
#include "vcfvconstraintscontext.hh"
#include "vcfvnewtonmethod.hh"

#include <ewoms/linear/vcfvlinearsolver.hh>
//#include <ewoms/linear/vcfvparallelamgsolver.hh>
#include <ewoms/nonlinear/newtonmethod.hh>
#include <ewoms/common/timemanager.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <limits>

namespace Ewoms {

// forward declaration
template<class TypeTag>
class VcfvModel;

namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for the time manager
SET_TYPE_PROP(VcfvModel, TimeManager, Ewoms::TimeManager<TypeTag>);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(VcfvModel,
              GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! Set the default for the ElementGeometry
SET_PROP(VcfvModel, FvElementGeometry)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Ewoms::VcfvElementGeometry<Scalar, GridView> type;
};

//! Mapper for the grid view's vertices.
SET_TYPE_PROP(VcfvModel,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);

//! Mapper for the grid view's elements.
SET_TYPE_PROP(VcfvModel,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(VcfvModel, DofMapper, typename GET_PROP_TYPE(TypeTag, VertexMapper));

//! Set the BaseLocalResidual to VcfvLocalResidual
SET_TYPE_PROP(VcfvModel, BaseLocalResidual, Ewoms::VcfvLocalResidual<TypeTag>);

//! Set the BaseModel to VcfvModel
SET_TYPE_PROP(VcfvModel, BaseModel, Ewoms::VcfvModel<TypeTag>);

//! The local jacobian operator for the VCVF discretization
SET_TYPE_PROP(VcfvModel, LocalJacobian, Ewoms::VcfvLocalJacobian<TypeTag>);

//! The maximum allowed number of timestep divisions for the
//! Newton solver
SET_INT_PROP(VcfvModel, MaxTimeStepDivisions, 10);

/*!
 * \brief A vector of quanties, each for one equation.
 */
SET_TYPE_PROP(VcfvModel,
              EqVector,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief A vector for mass/energy rates.
 *
 * E.g. Neumann fluxes or source terms
 */
SET_TYPE_PROP(VcfvModel,
              RateVector,
              typename GET_PROP_TYPE(TypeTag, EqVector));

/*!
 * \brief Type of object for specifying boundary conditions.
 */
SET_TYPE_PROP(VcfvModel,
              BoundaryRateVector,
              typename GET_PROP_TYPE(TypeTag, RateVector));

/*!
 * \brief The class which represents constraints.
 */
SET_TYPE_PROP(VcfvModel, Constraints, Ewoms::VcfvConstraints<TypeTag>);

/*!
 * \brief The type for storing a residual for an element.
 */
SET_TYPE_PROP(VcfvModel,
              ElementEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief The type for storing a residual for the whole grid.
 */
SET_TYPE_PROP(VcfvModel,
              GlobalEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief An object representing a local set of primary variables.
 */
SET_TYPE_PROP(VcfvModel,
              PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief The type of a solution for the whole grid at a fixed time.
 */
SET_TYPE_PROP(VcfvModel,
              SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

/*!
 * \brief The volume variable class.
 *
 * This should almost certainly be overloaded by the model...
 */
SET_TYPE_PROP(VcfvModel, VolumeVariables, Ewoms::VcfvVolumeVariables<TypeTag>);

/*!
 * \brief An array of secondary variable containers.
 */
SET_TYPE_PROP(VcfvModel, ElementContext, Ewoms::VcfvElementContext<TypeTag>);

/*!
 * \brief Assembler for the global jacobian matrix.
 */
SET_TYPE_PROP(VcfvModel, JacobianAssembler, Ewoms::VcfvAssembler<TypeTag>);

//! use an unlimited time step size by default
#if 0
// requires GCC 4.6 and above to call the constexpr function of
// numeric_limits
SET_SCALAR_PROP(VcfvModel, MaxTimeStepSize, std::numeric_limits<Scalar>::infinity());
#else
SET_SCALAR_PROP(VcfvModel, MaxTimeStepSize, 1e100);
#endif

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(VcfvModel, NumericDifferenceMethod, +1);

//! do not use hints by default
SET_BOOL_PROP(VcfvModel, EnableHints, false);

//! use FE gradients by default
SET_BOOL_PROP(VcfvModel, UseTwoPointGradients, false);

// disable jacobian matrix recycling by default
SET_BOOL_PROP(VcfvModel, EnableJacobianRecycling, false);

// disable partial reassembling by default
SET_BOOL_PROP(VcfvModel, EnablePartialReassemble, false);

// disable constraints by default
SET_BOOL_PROP(VcfvModel, EnableConstraints, false);

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(VcfvModel, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
public:
    typedef typename Dune::BCRSMatrix<MatrixBlock> type;
};

// use a parallel iterative solver by default
SET_TYPE_PROP(VcfvModel, LinearSolver, Ewoms::Linear::VcfvParallelSolver<TypeTag>);

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_SCALAR_PROP(VcfvModel, LinearSolverRelativeTolerance, 1e-6);

// the absolute defect of a component tolerated by the linear solver.
// By default, looking at the absolute defect is "almost" disabled.
SET_SCALAR_PROP(VcfvModel, LinearSolverAbsoluteTolerance, 1e-30);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(VcfvModel, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(VcfvModel, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));

//! set the size of the overlap region of the linear solver
SET_INT_PROP(VcfvModel, LinearSolverOverlapSize, 2);

//! Set the history size of the time discretiuation to 2 (for implicit euler)
SET_INT_PROP(VcfvModel, TimeDiscHistorySize, 2);

/*!
 * \brief Set the algorithm used for the linear solver.
 *
 * Possible choices are:
 * - \c SolverWrapperLoop: A fixpoint solver (using the Richardson iteration)
 * - \c SolverWrapperGradients: The steepest descent solver
 * - \c SolverWrapperCG: A conjugated gradients solver
 * - \c SolverWrapperBiCGStab: A stabilized bi-conjugated gradients solver
 * - \c SolverWrapperMinRes: A solver based on the  minimized residual algorithm
 * - \c SolverWrapperGMRes: A restarted GMRES solver
 */
SET_TYPE_PROP(VcfvModel,
              LinearSolverWrapper,
              Ewoms::Linear::SolverWrapperBiCGStab<TypeTag>);

/*!
 * \brief Set the algorithm used for as the preconditioner for the linear solver.
 *
 * Possible choices are:
 * - \c PreconditionerWrapperJacobi: A simple Jacobi preconditioner
 * - \c PreconditionerWrapperGaussSeidel: A Gauss-Seidel preconditioner
 * - \c PreconditionerWrapperSSOR: A symmetric successive overrelaxation (SSOR) preconditioner
 * - \c PreconditionerWrapperSOR: A successive overrelaxation (SOR) preconditioner
 * - \c PreconditionerWrapperILU: An ILU(n) preconditioner
 */
SET_TYPE_PROP(VcfvModel,
              PreconditionerWrapper,
              Ewoms::Linear::PreconditionerWrapperILU<TypeTag>);

} // namespace Properties
} // namespace Ewoms

#endif
