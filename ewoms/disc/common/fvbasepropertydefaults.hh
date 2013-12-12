// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2009-2013 by Andreas Lauser

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
 * \ingroup Discretization
 *
 * \brief Defines defaults for the common properties of the finite
 *        volume discretizations.
 */
#ifndef EWOMS_FV_BASE_PROPERTY_DEFAULTS_HH
#define EWOMS_FV_BASE_PROPERTY_DEFAULTS_HH

#include <ewoms/disc/common/fvbaseproperties.hh>
#include <ewoms/disc/common/fvbaseassembler.hh>
#include <ewoms/disc/common/fvbaselocaljacobian.hh>
#include <ewoms/disc/common/fvbaselocalresidual.hh>
#include <ewoms/disc/common/fvbaseelementcontext.hh>
#include <ewoms/disc/common/fvbaseboundarycontext.hh>
#include <ewoms/disc/common/fvbaseconstraintscontext.hh>
#include <ewoms/disc/common/fvbaseconstraints.hh>
#include <ewoms/disc/common/fvbasediscretization.hh>
#include <ewoms/disc/common/fvbasegradientcalculator.hh>
#include <ewoms/disc/common/fvbasenewtonmethod.hh>
#include <ewoms/disc/common/fvbasevolumevariables.hh>
#include <ewoms/disc/common/fvbasefluxvariables.hh>

#include <ewoms/linear/nullborderlistcreator.hh>
#include <ewoms/common/timemanager.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <limits>

namespace Opm {
namespace Properties {
//! Set the default type for the time manager
SET_TYPE_PROP(FvBaseDiscretization, TimeManager, Ewoms::TimeManager<TypeTag>);

//! Use the leaf grid view if not defined otherwise
SET_TYPE_PROP(FvBaseDiscretization, GridView,
              typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! Mapper for the grid view's vertices.
SET_TYPE_PROP(FvBaseDiscretization, VertexMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGVertexLayout>);

//! Mapper for the grid view's elements.
SET_TYPE_PROP(FvBaseDiscretization, ElementMapper,
              Dune::MultipleCodimMultipleGeomTypeMapper<typename GET_PROP_TYPE(TypeTag, GridView),
                                                        Dune::MCMGElementLayout>);

//! marks the border indices (required for the algebraic overlap stuff)
SET_PROP(FvBaseDiscretization, BorderListCreator)
{
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Ewoms::Linear::NullBorderListCreator<GridView, DofMapper> type;
};

SET_TYPE_PROP(FvBaseDiscretization, DiscLocalResidual, Ewoms::FvBaseLocalResidual<TypeTag>);

SET_TYPE_PROP(FvBaseDiscretization, DiscVolumeVariables, Ewoms::FvBaseVolumeVariables<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, DiscFluxVariables, Ewoms::FvBaseFluxVariables<TypeTag>);

//! Calculates the gradient of any quantity given the index of a flux approximation point
SET_TYPE_PROP(FvBaseDiscretization, GradientCalculator, Ewoms::FvBaseGradientCalculator<TypeTag>);

SET_TYPE_PROP(FvBaseDiscretization, DiscLocalJacobian, Ewoms::FvBaseLocalJacobian<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, LocalJacobian, typename GET_PROP_TYPE(TypeTag, DiscLocalJacobian));


//! Set the type of a global jacobian matrix from the solution types
SET_PROP(FvBaseDiscretization, JacobianMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef typename Dune::FieldMatrix<Scalar, numEq, numEq> MatrixBlock;
public:
    typedef typename Dune::BCRSMatrix<MatrixBlock> type;
};

//! The maximum allowed number of timestep divisions for the
//! Newton solver
SET_INT_PROP(FvBaseDiscretization, MaxTimeStepDivisions, 10);

/*!
 * \brief A vector of quanties, each for one equation.
 */
SET_TYPE_PROP(FvBaseDiscretization, EqVector,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief A vector for mass/energy rates.
 *
 * E.g. Neumann fluxes or source terms
 */
SET_TYPE_PROP(FvBaseDiscretization, RateVector,
              typename GET_PROP_TYPE(TypeTag, EqVector));

/*!
 * \brief Type of object for specifying boundary conditions.
 */
SET_TYPE_PROP(FvBaseDiscretization, BoundaryRateVector,
              typename GET_PROP_TYPE(TypeTag, RateVector));

/*!
 * \brief The class which represents constraints.
 */
SET_TYPE_PROP(FvBaseDiscretization, Constraints, Ewoms::FvBaseConstraints<TypeTag>);

/*!
 * \brief The type for storing a residual for an element.
 */
SET_TYPE_PROP(FvBaseDiscretization, ElementEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief The type for storing a residual for the whole grid.
 */
SET_TYPE_PROP(FvBaseDiscretization, GlobalEqVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, EqVector)>);

/*!
 * \brief An object representing a local set of primary variables.
 */
SET_TYPE_PROP(FvBaseDiscretization, PrimaryVariables,
              Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                GET_PROP_VALUE(TypeTag, NumEq)>);

/*!
 * \brief The type of a solution for the whole grid at a fixed time.
 */
SET_TYPE_PROP(FvBaseDiscretization, SolutionVector,
              Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

/*!
 * \brief The volume variable class.
 *
 * This should almost certainly be overloaded by the model...
 */
SET_TYPE_PROP(FvBaseDiscretization, VolumeVariables, Ewoms::FvBaseVolumeVariables<TypeTag>);

/*!
 * \brief The element context
 */
SET_TYPE_PROP(FvBaseDiscretization, ElementContext, Ewoms::FvBaseElementContext<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, BoundaryContext, Ewoms::FvBaseBoundaryContext<TypeTag>);
SET_TYPE_PROP(FvBaseDiscretization, ConstraintsContext, Ewoms::FvBaseConstraintsContext<TypeTag>);

/*!
 * \brief Assembler for the global jacobian matrix.
 */
SET_TYPE_PROP(FvBaseDiscretization, JacobianAssembler, Ewoms::FvBaseAssembler<TypeTag>);

//! use an unlimited time step size by default
#if 0
// requires GCC 4.6 and above to call the constexpr function of
// numeric_limits
SET_SCALAR_PROP(FvBaseDiscretization, MaxTimeStepSize, std::numeric_limits<Scalar>::infinity());
#else
SET_SCALAR_PROP(FvBaseDiscretization, MaxTimeStepSize, 1e100);
#endif
//! By default, accept any time step larger than zero
SET_SCALAR_PROP(FvBaseDiscretization, MinTimeStepSize, 0.0);

//! The base epsilon value for finite difference calculations
SET_SCALAR_PROP(FvBaseDiscretization, BaseEpsilon, 0.9123e-10);

//! use forward differences to calculate the jacobian by default
SET_INT_PROP(FvBaseDiscretization, NumericDifferenceMethod, +1);

//! do not use hints by default
SET_BOOL_PROP(FvBaseDiscretization, EnableHints, false);

// disable jacobian matrix recycling by default
SET_BOOL_PROP(FvBaseDiscretization, EnableJacobianRecycling, false);

// disable partial reassembling by default
SET_BOOL_PROP(FvBaseDiscretization, EnablePartialReassemble, false);

// disable constraints by default
SET_BOOL_PROP(FvBaseDiscretization, EnableConstraints, false);

// if the deflection of the newton method is large, we do not
// need to solve the linear approximation accurately. Assuming
// that the initial value for the delta vector u is quite
// close to the final value, a reduction of 6 orders of
// magnitude in the defect should be sufficient...
SET_SCALAR_PROP(FvBaseDiscretization, LinearSolverRelativeTolerance, 1e-6);

// the absolute defect of a component tolerated by the linear solver.
// By default, looking at the absolute defect is "almost" disabled.
SET_SCALAR_PROP(FvBaseDiscretization, LinearSolverAbsoluteTolerance, 1e-30);

//! set the default for the accepted fix-point tolerance (we use 0 to disable considering the fix-point tolerance)
SET_SCALAR_PROP(FvBaseDiscretization, LinearSolverFixPointTolerance, 0.0);

//! Set the history size of the time discretization to 2 (for implicit euler)
SET_INT_PROP(FvBaseDiscretization, TimeDiscHistorySize, 2);

//! Most models don't need the gradients at the center of the SCVs, so
//! we disable them by default.
SET_BOOL_PROP(FvBaseDiscretization, RequireScvCenterGradients, false);

}} // namespace Properties, Opm

#endif
