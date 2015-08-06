// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2014 by Andreas Lauser

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
 * \brief A general-purpose simulator for ECL decks using the black-oil model.
 */
#include "config.h"

#ifdef HAVE_DUNE_FEM
#warning "HACK: For now, dune-fem cannot be used in conjunction with Dune::CpGrid because CpGrid has some API conformity issues."
#undef HAVE_DUNE_FEM
#endif

#include <ewoms/common/quad.hh>
#include <ewoms/common/start.hh>
#include <ewoms/models/blackoil/blackoilmodel.hh>
#include "eclproblem.hh"

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(EclProblem, INHERITS_FROM(BlackOilModel, EclBaseProblem));
}}

int main(int argc, char **argv)
{
    typedef TTAG(EclProblem) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
