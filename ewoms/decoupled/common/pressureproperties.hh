// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
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
#ifndef EWOMS_PRESSURE_PROPERTIES_HH
#define EWOMS_PRESSURE_PROPERTIES_HH

#include "decoupledproperties.hh"

//Dune-includes
#include <ewoms/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <ewoms/linear/elementborderlistfromgrid.hh>
#include <ewoms/linear/paralleliterativebackend.hh>
#include <ewoms/linear/superlubackend.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

/*!
 * \ingroup Sequential
 */
/*!
 * \file
 * \brief Base file for properties related to sequential IMPET algorithms
 */
namespace Ewoms
{
namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPET algorithms
 *
 * This class holds properties necessary for the sequential IMPET solution.
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(Pressure, INHERITS_FROM(DecoupledModel));

NEW_PROP_TAG(LinearSolver);
SET_SPLICES(Pressure, LinearSolver);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//Properties for linear solvers
NEW_PROP_TAG(PressureCoefficientMatrix);//!< Type of the coefficient matrix given to the linear solver
NEW_PROP_TAG(PressureRHSVector);//!< Type of the right hand side vector given to the linear solver
NEW_PROP_TAG(PressureSolutionVector);//!Type of solution vector or pressure system
NEW_PROP_TAG( VisitFacesOnlyOnce); //!< Indicates if faces are only regarded from one side
NEW_PROP_TAG(Velocity);
}
}

#include <ewoms/decoupled/common/fv/fvvelocitydefault.hh>

namespace Ewoms
{
namespace Properties
{
//! Faces are only regarded from one side and not from both cells
SET_BOOL_PROP(Pressure, VisitFacesOnlyOnce, false);

//Set defaults
SET_PROP(Pressure, PressureCoefficientMatrix)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;

public:
    typedef Dune::BCRSMatrix<MB> type;
};
SET_PROP(Pressure, PressureRHSVector)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > type;
};

SET_TYPE_PROP(Pressure, PressureSolutionVector, typename GET_PROP(TypeTag, SolutionTypes)::ScalarSolution);

SET_TYPE_PROP(Pressure, JacobianMatrix, typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix));
SET_TYPE_PROP(Pressure, GlobalEqVector, typename GET_PROP_TYPE(TypeTag, PressureSolutionVector));

SET_PROP(Pressure, BorderListCreator)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ElementMapper ElementMapper;

public:
    typedef Linear::ElementBorderListFromGrid<GridView, ElementMapper> type;
};


// use the stabilized BiCG solver preconditioned by the ILU-0 by default
SET_TAG_PROP(Pressure, LinearSolver, ParallelIterativeLinearSolver);

//! For some reason, the semi-implicit models need to use the
//! reduction of the two-norm of the linear residual as their
//! convergence criterion.
SET_BOOL_PROP(Pressure, LinearSolverUseTwoNormReductionCriterion, true);

//! set the default for the reduction of the initial residual
SET_SCALAR_PROP(Pressure, LinearSolverRelativeTolerance, 1e-13);

//! set the default number of maximum iterations for the linear solver
SET_INT_PROP(Pressure, LinearSolverMaxIterations, 500);

SET_TYPE_PROP(Pressure, Velocity, FVVelocityDefault<TypeTag>);

}
}

#endif
