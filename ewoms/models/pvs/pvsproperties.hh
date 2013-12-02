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
 * \ingroup PvsModel
 *
 * \brief Declares the properties required for the compositional
 *        multi-phase primary variable switching VCVF discretization.
 */
#ifndef EWOMS_PVS_PROPERTIES_HH
#define EWOMS_PVS_PROPERTIES_HH

#include <ewoms/disc/vcfv/vcfvproperties.hh>
#include <ewoms/vtk/vcfvvtkmultiphasemodule.hh>
#include <ewoms/vtk/vcfvvtkcompositionmodule.hh>
#include <ewoms/vtk/vcfvvtkphasepresencemodule.hh>
#include <ewoms/vtk/vcfvvtktemperaturemodule.hh>
#include <ewoms/vtk/vcfvvtkdiffusionmodule.hh>
#include <ewoms/vtk/vcfvvtkenergymodule.hh>

namespace Opm {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal single phase problems
NEW_TYPE_TAG(VcfvPvs, INHERITS_FROM(VcfvModel, VtkPhasePresence, VtkMultiPhase,
                                    VtkComposition, VtkTemperature, VtkEnergy,
                                    VtkDiffusion));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);     //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents); //!< Number of fluid components in the system
NEW_PROP_TAG(PvsVerbosity);  //!< The verbosity of the model (0 -> do not print
// anything, 2 -> spam stdout a lot)
NEW_PROP_TAG(Indices);     //!< Enumerations used by the model
NEW_PROP_TAG(FluidSystem); //!< Provides the thermodynamic relations

NEW_PROP_TAG(MaterialLaw);       //!< The material law which ought to be used
NEW_PROP_TAG(MaterialLawParams); //!< The parameters of the material law

NEW_PROP_TAG(HeatConductionLaw); //!< The heat conduction law which ought to be
// used
NEW_PROP_TAG(HeatConductionLawParams); //!< The parameters of the heat
// conduction law

NEW_PROP_TAG(VelocityModule); //!< Specifies the relation used for velocity

NEW_PROP_TAG(EnableGravity); //!< Specifies whether gravity is considered in the
// problem
NEW_PROP_TAG(EnableSmoothUpwinding); //!< Specifies whether the smooth upwinding
// method should be used
NEW_PROP_TAG(EnableEnergy); //!< Specifies whether energy is considered as a
// conservation quantity or not
NEW_PROP_TAG(EnableDiffusion); //!< Enable diffusive fluxes?

NEW_PROP_TAG(PvsPressureBaseWeight); //!< The basis value for the weight of the
// pressure primary variable
NEW_PROP_TAG(PvsSaturationsBaseWeight); //!< The basis value for the weight of
// the saturation primary variables
NEW_PROP_TAG(PvsMoleFractionsBaseWeight); //!< The basis value for the weight of
// the mole fraction primary variables
} // namespace Properties
} // namespace Opm

#endif
