// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2009-2012 by Andreas Lauser

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
 * \ingroup ImmiscibleVcfvModel
 *
 * \brief Defines the properties required for the immiscible
 *        multi-phase VCVF discretization.
 */
#ifndef EWOMS_IMMISCIBLE_PROPERTIES_HH
#define EWOMS_IMMISCIBLE_PROPERTIES_HH

#include <ewoms/disc/vcfv/vcfvproperties.hh>
#include <ewoms/vtk/vcfvvtkmultiphasemodule.hh>
#include <ewoms/vtk/vcfvvtktemperaturemodule.hh>
#include <ewoms/vtk/vcfvvtkenergymodule.hh>

////////////////////////////////
// properties
////////////////////////////////
namespace Opm {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The generic type tag for problems using the immiscible multi-phase model
NEW_TYPE_TAG(VcfvImmiscible,
             INHERITS_FROM(VcfvModel, VtkMultiPhase, VtkTemperature, VtkEnergy));
//! The type tag for single-phase immiscible problems
NEW_TYPE_TAG(VcfvImmiscibleOnePhase, INHERITS_FROM(VcfvImmiscible));
//! The type tag for two-phase immiscible problems
NEW_TYPE_TAG(VcfvImmiscibleTwoPhase, INHERITS_FROM(VcfvImmiscible));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);     //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of chemical species in the system
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the
// problem
NEW_PROP_TAG(Indices);           //!< Enumerations used by the model
NEW_PROP_TAG(MaterialLaw);       //!< The material law which ought to be used
                                 //(extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The context material law (extracted from
// the spatial parameters)
NEW_PROP_TAG(HeatConductionLaw);       //!< The material law for heat conduction
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the material law
// for heat conduction
NEW_PROP_TAG(FluidSystem); //!<The fluid systems including the information about
// the phases
NEW_PROP_TAG(EnableEnergy); //!< Specify whether energy should be considered as
// a conservation quantity or not
NEW_PROP_TAG(EnableSmoothUpwinding); //!< Specifies whether the smooth upwinding
// method should be used

//! Specifies the relation used for velocity
NEW_PROP_TAG(VelocityModule);

// these properties only make sense for the VcfvImmiscibleTwoPhase type tag
NEW_PROP_TAG(WettingPhase);    //!< The wetting phase for two-phase models
NEW_PROP_TAG(NonwettingPhase); //!< The non-wetting phase for two-phase models

// these properties only make sense for the VcfvImmiscibleOnePhase type tag
NEW_PROP_TAG(Fluid); //!< The fluid used by the model

} // namespace Properties
} // namespace Opm

#endif
