// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Klaus Mosthaf                                     *
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \ingroup BoxStokesModel
 *
 * \brief Declares the properties required by the Stokes box model.
 */
#ifndef EWOMS_STOKES_PROPERTIES_HH
#define EWOMS_STOKES_PROPERTIES_HH

#include <ewoms/boxmodels/common/boxproperties.hh>

namespace Ewoms {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the problems using the Stokes equations
NEW_TYPE_TAG(BoxStokes, INHERITS_FROM(BoxModel));

/*!
 * \brief The type tag for the problems using the Navier-Stokes equations.
 *
 * Basically this just takes everything from the Stokes model, but
 * sets the \c EnableNavierTerm property to true by default.
 */
NEW_TYPE_TAG(BoxNavierStokes, INHERITS_FROM(BoxStokes));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(Indices); //!< Enumerations for the Stokes models accessible using a generic name
NEW_PROP_TAG(Fluid);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(HeatConductionLaw);
NEW_PROP_TAG(HeatConductionLawParams);
NEW_PROP_TAG(StokesPhaseIndex); //!< The index of the considered fluid phase
NEW_PROP_TAG(BaseProblem); //!< The type of the base class for problems using the stokes model
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(EnableEnergy); //!< Specify whether the energy equation is enabled or not
NEW_PROP_TAG(EnableNavierTerm); //!< Specify whether the model should feature an inertial term or not (default: false)
}
}

#endif
