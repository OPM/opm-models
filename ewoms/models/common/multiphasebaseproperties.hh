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

#include <ewoms/disc/vcfv/vcfvdiscretization.hh>
#include <ewoms/vtk/vtkmultiphasemodule.hh>
#include <ewoms/vtk/vtktemperaturemodule.hh>

namespace Opm {
namespace Properties {
//! The generic type tag for problems using the immiscible multi-phase model
NEW_TYPE_TAG(MultiPhaseBaseModel, INHERITS_FROM(VtkMultiPhase, VtkTemperature));

//! The property for the spatial discretization splice
NEW_PROP_TAG(SpatialDiscretizationSplice);

//! Specify the splices of the MultiPhaseBaseModel type tag
SET_SPLICES(MultiPhaseBaseModel, SpatialDiscretizationSplice);

//! Set the default spatial discretization
//!
//! We use a vertex centered finte volume method by default
SET_TAG_PROP(MultiPhaseBaseModel, SpatialDiscretizationSplice, VcfvDiscretization);

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents);   //!< Number of chemical species in the system
NEW_PROP_TAG(Indices); //!< Enumerations used by the model
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from the spatial parameters)
NEW_PROP_TAG(HeatConductionLaw); //!< The material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the material law for heat conduction
NEW_PROP_TAG(FluidSystem); //!<The fluid systems including the information about the phases
NEW_PROP_TAG(VelocityModule); //! Specifies the relation used for velocity

NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(EnableSmoothUpwinding); //!< Specify whether the "upwinding kink" should be avoided
} // namespace Properties
} // namespace Opm

#endif
