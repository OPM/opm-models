// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Markus Wolff                                 *
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
#ifndef EWOMS_IMPES2PADAPTIVE_PROPERTIES_HH
#define EWOMS_IMPES2PADAPTIVE_PROPERTIES_HH

#include <ewoms/decoupled/common/impetproperties.hh>
#include <ewoms/decoupled/2p/2pproperties.hh>

/*!
 * \ingroup IMPES
 */
/*!
 * \file
 * \brief Properties for adaptive implementations of the sequential IMPES algorithms
 */
namespace Ewoms
{

namespace Properties
{
/*!
 *
 * \brief General properties for adaptive implementations of the sequential IMPES algorithms
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//!  TypeTag for grid-adaptive two-phase IMPES scheme
NEW_TYPE_TAG(IMPESTwoPAdaptive, INHERITS_FROM(IMPET, DecoupledTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
}
}

#include <ewoms/decoupled/common/variableclassadaptive.hh>
#include <ewoms/decoupled/2p/cellData2padaptive.hh>
#include "gridadaptionindicator2p.hh"
#include <ewoms/decoupled/2p/impes/impesproblem2p.hh>
#include <ewoms/decoupled/common/gridadaptinitializationindicator.hh>

namespace Ewoms
{
namespace Properties
{
//! Enable adaptive grid
SET_BOOL_PROP(IMPESTwoPAdaptive, EnableGridAdapt, true);
//! Set variable class for adaptive impet schemes
SET_TYPE_PROP(IMPESTwoPAdaptive, Variables, Ewoms::VariableClassAdaptive<TypeTag>);
//! Set cell data class for adaptive two-phase IMPES schemes
SET_TYPE_PROP(IMPESTwoPAdaptive, CellData, Ewoms::CellData2PAdaptive<TypeTag>);
//! Set the standard indicator class of two-phase models for adaption or coarsening
SET_TYPE_PROP(IMPESTwoPAdaptive, GridAdaptIndicator, GridAdaptionIndicator2P<TypeTag>);
//!Set default class for adaptation initialization indicator
SET_TYPE_PROP(IMPESTwoPAdaptive,  GridAdaptInitializationIndicator, GridAdaptInitializationIndicator<TypeTag>);
}
}

#endif
