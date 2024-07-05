// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::BlackOilNewtonMethod
 */
#ifndef EWOMS_BLACK_OIL_NEWTON_METHOD_PARAMETERS_HH
#define EWOMS_BLACK_OIL_NEWTON_METHOD_PARAMETERS_HH

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Parameters {

template<class TypeTag, class MyTypeTag>
struct DpMaxRel { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct DsMax { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct PriVarOscilationThreshold { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct ProjectSaturations { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct MaxTemperatureChange { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TemperatureMax { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct TemperatureMin { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct PressureMax { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct PressureMin { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct MaximumWaterSaturation { using type = Properties::UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct WaterOnlyThreshold { using type = Properties::UndefinedProperty; };

} // namespace Opm::Parameters

#endif