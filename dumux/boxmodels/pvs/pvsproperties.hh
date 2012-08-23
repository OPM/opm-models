// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Melanie Darcis                               *
 *   Copyright (C) 2009-2012 by Klaus Mosthaf                                *
 *   Copyright (C) 2009-2012 by Bernd Flemisch                               *
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
 * \ingroup PvsModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the PVS BOX model.
 */
#ifndef DUMUX_PVS_PROPERTIES_HH
#define DUMUX_PVS_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>
#include <dumux/boxmodels/vtk/boxvtkmultiphasemodule.hh>
#include <dumux/boxmodels/vtk/boxvtkcompositionmodule.hh>
#include <dumux/boxmodels/vtk/boxvtkphasepresencemodule.hh>
#include <dumux/boxmodels/vtk/boxvtktemperaturemodule.hh>
#include <dumux/boxmodels/vtk/boxvtkenergymodule.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal single phase problems
NEW_TYPE_TAG(BoxPvs, INHERITS_FROM(BoxModel, VtkPhasePresence, VtkMultiPhase, VtkComposition, VtkTemperature, VtkEnergy));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(PvsVerbosity); //!< The verbosity of the model (0 -> do not print anything, 2 -> spam stdout a lot)
NEW_PROP_TAG(Indices); //!< Enumerations used by the model
NEW_PROP_TAG(FluidSystem); //!< Provides the thermodynamic relations

NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used
NEW_PROP_TAG(MaterialLawParams); //!< The parameters of the material law

NEW_PROP_TAG(HeatConductionLaw);   //!< The heat conduction law which ought to be used
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the heat conduction law

NEW_PROP_TAG(EnableGravity); //!< Specifies whether gravity is considered in the problem
NEW_PROP_TAG(EnableSmoothUpwinding); //!< Specifies whether the smooth upwinding method should be used
NEW_PROP_TAG(EnableEnergy); //!< Specifies whether energy is considered as a conservation quantity or not
}
}

#endif
