// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
 * \ingroup StokesModel
 *
 * \brief Declares the properties required by the Stokes model.
 */
#ifndef EWOMS_STOKES_PROPERTIES_HH
#define EWOMS_STOKES_PROPERTIES_HH

#include <ewoms/disc/vcfv/vcfvdiscretization.hh>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Indices); //!< Enumerations for the Stokes models accessible using
                       //!< a generic name
NEW_PROP_TAG(Fluid);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(HeatConductionLaw);
NEW_PROP_TAG(HeatConductionLawParams);
NEW_PROP_TAG(StokesPhaseIndex); //!< The index of the considered fluid phase
NEW_PROP_TAG(BaseProblem); //!< The type of the base class for problems using
                           //!< the stokes model
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the
                             //!< problem
NEW_PROP_TAG(EnableEnergy); //!< Specify whether the energy equation is enabled
                            //!< or not
NEW_PROP_TAG(EnableNavierTerm); //!< Specify whether the model should feature an
                                //!< inertial term or not (default: false)
} // namespace Properties
} // namespace Opm

#endif
