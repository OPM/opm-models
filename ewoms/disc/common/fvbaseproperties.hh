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
 * \brief Declare the properties used by the infrastructure code of
 *        the finite volume discretizations.
 */
#ifndef EWOMS_FV_BASE_PROPERTIES_HH
#define EWOMS_FV_BASE_PROPERTIES_HH

#include "fvbasenewtonmethod.hh"
#include "fvbaseproperties.hh"

#include <ewoms/common/basicproperties.hh>
#include <ewoms/vtk/vtkprimaryvarsmodule.hh>
#include <ewoms/linear/paralleliterativebackend.hh>

namespace Opm {
namespace Properties {
/*!
 * \ingroup Discretization
 */
// \{
//! The type tag for models based on the ECFV-scheme
NEW_TYPE_TAG(FvBaseDiscretization,
             INHERITS_FROM(ImplicitModel,
                           FvBaseNewtonMethod,
                           VtkPrimaryVars));


//! set the splices for the vertex-centered finite volume
//! discretization
NEW_PROP_TAG(LinearSolverSplice);
SET_SPLICES(FvBaseDiscretization, LinearSolverSplice);

//! use a parallel iterative linear solver by default
SET_TAG_PROP(FvBaseDiscretization, LinearSolverSplice, ParallelIterativeLinearSolver);

NEW_PROP_TAG(Grid); //!< The type of the DUNE grid
NEW_PROP_TAG(GridView); //!< The type of the grid view

NEW_PROP_TAG(Stencil); //!< The class describing the stencil of the spatial discretization

NEW_PROP_TAG(Problem); //!< The type of the problem
NEW_PROP_TAG(BaseProblem);   //!< The type of the base class for all problems which use this model
NEW_PROP_TAG(Model); //!< The type of the model
NEW_PROP_TAG(NumEq); //!< Number of equations in the system of PDEs

NEW_PROP_TAG(Discretization); //!< The type of the spatial discretization used by the model
NEW_PROP_TAG(DiscLocalResidual); //!< The discretization specific part of the local residual
NEW_PROP_TAG(LocalResidual); //!< The type of the local residual function
NEW_PROP_TAG(DiscLocalJacobian); //!< The discretization specific part of the local jacobian of the residual
NEW_PROP_TAG(LocalJacobian); //!< The type of the local jacobian operator

NEW_PROP_TAG(BaseJacobianAssembler); //!< Assembles the global jacobian matrix
NEW_PROP_TAG(JacobianMatrix); //!< Type of the global jacobian matrix

NEW_PROP_TAG(EqVector); //!< A vector of holding a quantity for each equation (usually at a given spatial location)
NEW_PROP_TAG(ElementEqVector); //!< A vector of holding a quantity for each equation for each DOF of an element
NEW_PROP_TAG(GlobalEqVector); //!< Vector containing a quantity of for equation for each DOF of the whole grid

NEW_PROP_TAG(RateVector); //!< Vector containing volumetric or areal rates of quantities
NEW_PROP_TAG(BoundaryRateVector); //!< Type of object for specifying boundary conditions
NEW_PROP_TAG(Constraints); //!< The class which represents a constraint degree of freedom

NEW_PROP_TAG(SolutionVector); //!< Vector containing all primary variables of the grid

NEW_PROP_TAG(PrimaryVariables); //!< A vector of primary variables within a sub-control volume
NEW_PROP_TAG(VolumeVariables);  //!< The secondary variables within a sub-control volume
NEW_PROP_TAG(DiscVolumeVariables); //!< The discretization specific part of the volume variables

NEW_PROP_TAG(ElementContext); //!< The secondary variables of all degrees of freedom in an element's stencil
NEW_PROP_TAG(BoundaryContext); //!< The secondary variables of a boundary segment
NEW_PROP_TAG(ConstraintsContext); //!< The secondary variables of a constraint degree of freedom
NEW_PROP_TAG(FluxVariables); //!< Data required to calculate a flux over a face
NEW_PROP_TAG(GradientCalculator); //!< Calculates gradients of arbitrary quantities at flux integration points

NEW_PROP_TAG(DiscBaseVolumeVariables); //!< The part of the volume variables which is specific to the spatial discretization
NEW_PROP_TAG(DiscFluxVariables); //!< The part of the flux variables which is specific to the spatial discretization

NEW_PROP_TAG(DiscVtkOutputModule); //!< The part of the VTK ouput modules which is specific to the spatial discretization

NEW_PROP_TAG(GridCommHandleFactory); //!< The class to create grid communication handles

// high level simulation control
NEW_PROP_TAG(TimeManager);  //!< Manages the simulation time

//! Specify whether the jacobian matrix of the last iteration of a
//! time step should be re-used as the jacobian of the first iteration
//! of the next time step.
NEW_PROP_TAG(EnableJacobianRecycling);

//! Specify whether the jacobian matrix should be only reassembled for
//! elements where at least one degree of freedom is above the
//! specified tolerance
NEW_PROP_TAG(EnablePartialReassemble);

//! Specify whether the some degrees of fredom can be constraint
NEW_PROP_TAG(EnableConstraints);

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
NEW_PROP_TAG(MaxTimeStepSize);

/*!
 * \brief Specify the minimal size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
NEW_PROP_TAG(MinTimeStepSize);

/*!
 * \brief The maximum allowed number of timestep divisions for the
 *        Newton solver.
 */
NEW_PROP_TAG(MaxTimeStepDivisions);

/*!
 * \brief Specify which kind of method should be used to numerically
 * calculate the partial derivatives of the residual.
 *
 * -1 means backward differences, 0 means central differences, 1 means
 * forward differences. By default we use central differences.
 */
NEW_PROP_TAG(NumericDifferenceMethod);

/*!
 * \brief Specify whether to use the already calculated solutions as
 *        starting values of the volume variables.
 *
 * This only makes sense if the calculation of the volume variables is
 * very expensive (e.g. for non-linear fugacity functions where the
 * solver converges faster).
 */
NEW_PROP_TAG(EnableHints);

//! The base epsilon value for finite difference calculations
NEW_PROP_TAG(BaseEpsilon);

// mappers from local to global indices

/*!
 * \brief The mapper to find the global index of a vertex.
 */
NEW_PROP_TAG(VertexMapper);

/*!
 * \brief The mapper to find the global index of an element.
 */
NEW_PROP_TAG(ElementMapper);

/*!
 * \brief The mapper to find the global index of a degree of freedom.
 */
NEW_PROP_TAG(DofMapper);

/*!
 * \brief The class which marks the border indices associated with the
 *        degrees of freedom on a process boundary.
 *
 * This is required for the algebraic overlap stuff.
 */
NEW_PROP_TAG(BorderListCreator);

/*!
 * \brief The history size required by the time discretization
 */
NEW_PROP_TAG(TimeDiscHistorySize);

/*!
 * \brief Specify whether the gradients in the center of the SCVs need
 *        to be updated.
 *
 * Most models don't need this, but the (Navier-)Stokes ones do...
 */
NEW_PROP_TAG(RequireScvCenterGradients);

}} // namespace Properties, Opm

// \}

#endif
