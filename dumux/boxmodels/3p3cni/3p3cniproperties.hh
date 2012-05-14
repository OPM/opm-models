// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \ingroup BoxProperties
 * \ingroup ThreePThreeCNIModel
 * \file
 *
 * \brief Defines the properties required for the non-isothermal three-phase,
 * three-component BOX model.
 */
#ifndef DUMUX_3P3CNI_PROPERTIES_HH
#define DUMUX_3P3CNI_PROPERTIES_HH

#include <dumux/boxmodels/3p3c/3p3cproperties.hh>
#include <dumux/boxmodels/vtk/boxvtkenergymodule.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the non-isothermal three-phase, three-component problems
NEW_TYPE_TAG(BoxThreePThreeCNI, INHERITS_FROM(BoxThreePThreeC, VtkEnergy));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(Indices); //!< Enumerations used by the model

NEW_PROP_TAG(HeatConductionLaw);   //!< The heat conduction law which ought to be used
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the heat conduction law
}
}

#endif
