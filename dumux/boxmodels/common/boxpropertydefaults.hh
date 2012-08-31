// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2011 by Klaus Mosthaf                                     *
 *   Copyright (C) 2012 by Benjamin Faigle                                   *
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
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup BoxModel
 *
 * \brief Default properties for box models
 */
#ifndef DUMUX_BOX_PROPERTY_DEFAULTS_HH
#define DUMUX_BOX_PROPERTY_DEFAULTS_HH

#include "boxproperties.hh"
#include "boxassembler.hh"
#include "boxmodel.hh"
#include "boxfvelementgeometry.hh"
#include "boxlocalresidual.hh"
#include "boxlocaljacobian.hh"
#include "boxelementcontext.hh"
#include "boxvolumevariables.hh"
#include "boxconstraints.hh"
#include "boxconstraintscontext.hh"
#include "boxnewtoncontroller.hh"

#include <dumux/linear/boxlinearsolver.hh>
//#include <dumux/linear/boxparallelamgsolver.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/common/timemanager.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <limits>

namespace Dumux {

// forward declaration
template<class TypeTag>
class BoxModel;

namespace Properties {
//////////////////////////////////////////////////////////////////
// Some defaults for very fundamental properties
//////////////////////////////////////////////////////////////////

//! Set the default type for the time manager
SET_TYPE_PROP(BoxModel, TimeManager, Dumux::TimeManager<TypeTag>);

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(BoxModel,
              GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! Set the default for the FVElementGeometry
SET_PROP(BoxModel, FVElementGeometry)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Dumux::BoxFVElementGeometry<Scalar, GridView> type;
};

//! use the plain newton method for the box scheme by default
SET_TYPE_PROP(BoxModel, NewtonMethod, Dumux::NewtonMethod<TypeTag>);

//! use the plain newton controller for the box scheme by default
SET_TYPE_PROP(BoxModel, NewtonController, Dumux::BoxNewtonController<TypeTag>);

//! Mapper for the grid view's vertices.
SET_TYPE_PROP(BoxModel,
              VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);

//! Mapper for the grid view's elements.
SET_TYPE_PROP(BoxModel,
              ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(BoxModel, DofMapper, typename GET_PROP_TYPE(TypeTag, VertexMapper));

//! Set the BaseLocalResidual to BoxLocalResidual
SET_TYPE_PROP(BoxModel, BaseLocalResidual, Dumux::BoxLocalResidual<TypeTag>);

//! Set the BaseModel to BoxModel
SET_TYPE_PROP(BoxModel, BaseModel, Dumux::BoxModel<TypeTag>);

//! The local jacobian operator for the box scheme
SET_TYPE_PROP(BoxModel, LocalJacobian, Dumux::BoxLocalJacobian<TypeTag>);

/*!
 * \brief A vector of quanties, each for one equation.
 */
SET_TYPE_PROP(BoxModel,
              EqVector,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief A vector for mass/energy rates.
 *
 * E.g. Neumann fluxes or source terms
 */
SET_TYPE_PROP(BoxModel,
              RateVector,
              typename GET_PROP_TYPE(TypeTag, EqVector));

/*!
 * \brief Type of object for specifying boundary conditions.
 */
SET_TYPE_PROP(BoxModel,
              BoundaryRateVector,
              typename GET_PROP_TYPE(TypeTag, RateVector));

/*!
 * \brief The class which represents constraints.
 */
SET_TYPE_PROP(BoxModel, Constraints, Dumux::BoxConstraints<TypeTag>);

/*!
 * \brief The type for storing a residual for an element.
 */
SET_TYPE_PROP(BoxModel,
              ElementEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief The type for storing a residual for the whole grid.
 */
SET_TYPE_PROP(BoxModel,
              GlobalEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief An object representing a local set of primary variables.
 */
SET_TYPE_PROP(BoxModel,
              PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief The type of a solution for the whole grid at a fixed time.
 */
SET_TYPE_PROP(BoxModel,
              SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

/*!
 * \brief The volume variable class.
 *
 * This should almost certainly be overloaded by the model...
 */
SET_TYPE_PROP(BoxModel, VolumeVariables, Dumux::BoxVolumeVariables<TypeTag>);

/*!
 * \brief An array of secondary variable containers.
 */
SET_TYPE_PROP(BoxModel, ElementContext, Dumux::BoxElementContext<TypeTag>);

/*!
 * \brief Assembler for the global jacobian matrix.
 */
SET_TYPE_PROP(BoxModel, JacobianAssembler, Dumux::BoxAssembler<TypeTag>);

//! use an unlimited time step size by default
#if 0
// requires GCC 4.6 and above to call the constexpr function of
// numeric_limits
SET_SCALAR_PROP(BoxModel, MaxTimeStepSize, std::numeric_limits<Scalar>::infinity());
#else
SET_SCALAR_PROP(BoxModel, MaxTimeStepSize, 1e100);
#endif

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(BoxModel, NumericDifferenceMethod, +1);

//! do not use hints by default
SET_BOOL_PROP(BoxModel, EnableHints, false);

//! use FE gradients by default
SET_BOOL_PROP(BoxModel, UseTwoPointGradients, false);

// disable jacobian matrix recycling by default
SET_BOOL_PROP(BoxModel, EnableJacobianRecycling, false);

// disable partial reassembling by default
SET_BOOL_PROP(BoxModel, EnablePartialReassemble, false);

// disable constraints by default
SET_BOOL_PROP(BoxModel, EnableConstraints, false);

//! Set the type of a global jacobian matrix from the solution types
SET_PROP(BoxModel, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
public:
    typedef typename Dune::BCRSMatrix<MatrixBlock> type;
};

// use a parallel iterative solver by default
SET_TYPE_PROP(BoxModel, LinearSolver, Dumux::Linear::BoxParallelSolver<TypeTag>);

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_SCALAR_PROP(BoxModel, LinearSolverTolerance, 1e-6);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(BoxModel, LinearSolverMaxIterations, 250);

//! set number of equations of the mathematical model as default
SET_INT_PROP(BoxModel, LinearSolverBlockSize, GET_PROP_VALUE(TypeTag, NumEq));

//! set the size of the overlap region of the linear solver
SET_INT_PROP(BoxModel, LinearSolverOverlapSize, 2);

//! Set the history size of the time discretiuation to 2 (for implicit euler)
SET_INT_PROP(BoxModel, TimeDiscHistorySize, 2);

/*!
 * \brief Set the algorithm used for the linear solver.
 *
 * Possible choices are:
 * - SolverWrapperLoop: A fixpoint solver
 * - SolverWrapperGradients: Steepest descent
 * - SolverWrapperCG: Conjugated gradients
 * - SolverWrapperBiCGStab: The stabilized bi-conjugated gradients
 * - SolverWrapperMinRes: The minimized residual algorithm
 * - SolverWrapperGMRes: A restarted GMRES solver
 */
SET_TYPE_PROP(BoxModel, 
              LinearSolverWrapper,
              Dumux::Linear::SolverWrapperBiCGStab<TypeTag>);

/*!
 * \brief Set the algorithm used for as the preconditioner for the linear solver.
 *
 * Possible choices are:
 * - PreconditionerWrapperJacobi: A simple Jacobi preconditioner
 * - PreconditionerWrapperGaussSeidel: A Gauss-Seidel preconditioner
 * - PreconditionerWrapperSSOR: A symmetric successive overrelaxation (SSOR) preconditioner
 * - PreconditionerWrapperSOR: A successive overrelaxation (SOR) preconditioner
 * - PreconditionerWrapperILU: An ILU(n) preconditioner
 */
SET_TYPE_PROP(BoxModel, 
              PreconditionerWrapper,
              Dumux::Linear::PreconditionerWrapperILU<TypeTag>);

} // namespace Properties
} // namespace Dumux

#endif
