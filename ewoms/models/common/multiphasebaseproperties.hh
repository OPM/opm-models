// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2013 by Andreas Lauser

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
 * \ingroup MultiPhaseBaseModel
 *
 * \brief Defines the common properties required by the porous medium
 *        multi-phase models.
 */
#ifndef EWOMS_MULTI_PHASE_BASE_PROPERTIES_HH
#define EWOMS_MULTI_PHASE_BASE_PROPERTIES_HH

#include <ewoms/disc/common/fvbaseproperties.hh>
#include <ewoms/io/vtkmultiphasemodule.hh>
#include <ewoms/io/vtktemperaturemodule.hh>

namespace Ewoms {
namespace Properties {
//! The splice to be used for the spatial discretization
NEW_PROP_TAG(SpatialDiscretizationSplice);
//! Number of fluid phases in the system
NEW_PROP_TAG(NumPhases);
//! Number of chemical species in the system
NEW_PROP_TAG(NumComponents);
//! Enumerations used by the model
NEW_PROP_TAG(Indices);
//! The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLaw);
//! The context material law (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams);
//! The material law for heat conduction
NEW_PROP_TAG(HeatConductionLaw);
//! The parameters of the material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams);
//!The fluid systems including the information about the phases
NEW_PROP_TAG(FluidSystem);
//! Specifies the relation used for velocity
NEW_PROP_TAG(FluxModule);

//! Returns whether gravity is considered in the problem
NEW_PROP_TAG(EnableGravity);
} // namespace Properties
} // namespace Ewoms

#endif
