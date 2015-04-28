/*
  Copyright (C) 2010-2013 by Andreas Lauser
  Copyright (C) 2011 by Bernd Flemisch

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
 *
 * \brief Defines a type tags and some fundamental properties all models.
 */
#ifndef EWOMS_BASIC_PROPERTIES_HH
#define EWOMS_BASIC_PROPERTIES_HH

#include <dune/common/parametertree.hh>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/io/dgfgridmanager.hh>

namespace Ewoms {
namespace Properties {
///////////////////////////////////
// Type tag definitions:
//
// NumericModel
// |
// +-> ImplicitModel
///////////////////////////////////

//! Type tag for all models.
NEW_TYPE_TAG(NumericModel, INHERITS_FROM(ParameterSystem));

//! Type tag for all fully coupled models.
NEW_TYPE_TAG(ImplicitModel, INHERITS_FROM(NumericModel));

///////////////////////////////////
// Property names which are always available:
//
// Scalar
///////////////////////////////////

//! Property to specify the type of scalar values.
NEW_PROP_TAG(Scalar);

//! Property which provides a Dune::ParameterTree.
NEW_PROP_TAG(ParameterTree);

//! Property which defines the group that is queried for parameters by default
NEW_PROP_TAG(ModelParameterGroup);

//! Property which provides a GridManager (manages grids)
NEW_PROP_TAG(GridManager);

//! Property provides the name of the file from which the grid ought to be
//! loaded from
NEW_PROP_TAG(GridFile);

//! Property which tells the GridManager how often the grid should be refined
//! after creation.
NEW_PROP_TAG(GridGlobalRefinements);

//! Property provides the name of the file from which the additional runtime
//! parameters should to be loaded from
NEW_PROP_TAG(ParameterFile);

/*!
 * \brief Print all properties on startup?
 *
 * 0 means 'no', 1 means 'yes', 2 means 'print only to logfiles'. The
 * default is 2.
 */
NEW_PROP_TAG(PrintProperties);

/*!
 * \brief Print all parameters on startup?
 *
 * 0 means 'no', 1 means 'yes', 2 means 'print only to logfiles'. The
 * default is 2.
 */
NEW_PROP_TAG(PrintParameters);

//! The default value for the simulation's end time
NEW_PROP_TAG(EndTime);

//! The default value for the simulation's initial time step size
NEW_PROP_TAG(InitialTimeStepSize);

//! The default value for the simulation's restart time
NEW_PROP_TAG(RestartTime);

///////////////////////////////////
// Default values for properties:
//
// Scalar -> double
///////////////////////////////////

//! Set the default type of scalar values to double
SET_TYPE_PROP(NumericModel, Scalar, double);

//! Set the ParameterTree property
SET_PROP(NumericModel, ParameterTree)
{
    typedef Dune::ParameterTree type;

    static Dune::ParameterTree &tree()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }
};

//! use the global group as default for the model's parameter group
SET_STRING_PROP(NumericModel, ModelParameterGroup, "");

//! Use the DgfGridManager by default
SET_TYPE_PROP(NumericModel, GridManager, Ewoms::DgfGridManager<TypeTag>);

//! Set a value for the GridFile property
SET_STRING_PROP(NumericModel, GridFile, "");

//! Use the leaf grid view by default.
//!
//! Except for spatial refinement, there is rarly a reason to use
//! anything else...
SET_TYPE_PROP(NumericModel, GridView, typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);

//! Set a value for the ParameterFile property
SET_STRING_PROP(NumericModel, ParameterFile, "");

//! Set the number of refinement levels of the grid to 0. This does not belong
//! here, strictly speaking.
SET_INT_PROP(NumericModel, GridGlobalRefinements, 0);

//! By default, print the properties on startup
SET_INT_PROP(NumericModel, PrintProperties, 2);

//! By default, print the values of the run-time parameters on startup
SET_INT_PROP(NumericModel, PrintParameters, 2);

//! The default value for the simulation's end time
SET_SCALAR_PROP(NumericModel, EndTime, -1e100);

//! The default value for the simulation's initial time step size
SET_SCALAR_PROP(NumericModel, InitialTimeStepSize, -1e100);

//! The default value for the simulation's restart time
SET_SCALAR_PROP(NumericModel, RestartTime, -1e100);

} // namespace Properties
} // namespace Ewoms

#endif
