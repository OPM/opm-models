// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2011-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2011-2012 by Klaus Mosthaf                                *
 *   Copyright (C) 2011-2012 by Markus Wolff                                 *
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
 * \ingroup Properties
 * \ingroup Linear
 * \file
 *
 * \brief Defines a type tag and some fundamental properties for
 *        linear solvers
 */
#ifndef DUMUX_LINEAR_SOLVER_PROPERTIES_HH
#define DUMUX_LINEAR_SOLVER_PROPERTIES_HH

#include <dumux/common/basicproperties.hh>

namespace Dumux
{
namespace Properties
{
//! Linear solver type tag for all models.
NEW_TYPE_TAG(LinearSolverTypeTag);

//! The type of the linear solver to be used
NEW_PROP_TAG(LinearSolver);

/*!
 * \brief Specifies the verbosity of the linear solver
 *
 * By default it is 0, i.e. it doesn't print anything. Setting this
 * property to 1 prints aggregated convergence rates, 2 prints the
 * convergence rate of every iteration of the scheme.
 */
NEW_PROP_TAG(LinearSolverVerbosity);

//! target tolerance of the initial residual
NEW_PROP_TAG(LinearSolverTolerance);

//! maximum number of iterations of solver
NEW_PROP_TAG(LinearSolverMaxIterations);

//! relaxation parameter for the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

/*!
 * \brief the order of the preconditioner.
 *
 * for some preconditioners, that means the number of iterations to be
 * performed per iteration of the linear solver, for others it means
 * more sophisticated things (e.g. for the ILU preconditioner).
 */
NEW_PROP_TAG(PreconditionerOrder);

//! restart parameter for GMRes
NEW_PROP_TAG(GMResRestart);

//! Size of the matrix/vector blocks
/*!
 * The number of different types of equations which build the system of equations to solve
 * can differ from the number of equations given by the mathematical/physical model (e.g. IMPES).
 * Thus, the block size does not have to be equal to NumEq.
 * (Especially important for the SuperLU solver!)
 */
NEW_PROP_TAG(LinearSolverBlockSize);

SET_INT_PROP(LinearSolverTypeTag, LinearSolverVerbosity, 0);

//! set the preconditioner relaxation parameter to 1.0 by default
SET_SCALAR_PROP(LinearSolverTypeTag, PreconditionerRelaxation, 1.0);

//! set the preconditioner order to 0 by default
SET_INT_PROP(LinearSolverTypeTag, PreconditionerOrder, 0);

//! set the GMRes restart parameter to 10 by default
SET_INT_PROP(LinearSolverTypeTag, GMResRestart, 10);

} // namespace Properties
} // namespace Dumux

#endif
