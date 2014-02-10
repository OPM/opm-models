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
 * \ingroup RichardsModel
 *
 * \brief Contains the property declarations for the Richards model.
 */
#ifndef EWOMS_RICHARDS_PROPERTIES_HH
#define EWOMS_RICHARDS_PROPERTIES_HH

#include <ewoms/models/common/multiphasebaseproperties.hh>

// \{
namespace Opm {
namespace Properties {
//! Index of the fluid which represents the wetting phase
NEW_PROP_TAG(LiquidPhaseIndex);

//! The fluid system used for the problem
NEW_PROP_TAG(FluidSystem);
//! The component which constitutes the wetting phase of the problem
NEW_PROP_TAG(WettingPhase);
//! The component which constitutes the non-wetting phase of the problem
NEW_PROP_TAG(NonwettingPhase);
// \}
}} // namespace Properties, Opm

#endif
