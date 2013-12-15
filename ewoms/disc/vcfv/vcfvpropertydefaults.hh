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
 * \ingroup VcfvDiscretization
 *
 * \brief Defines defaults for the common properties of the VCFV discretizations.
 */
#ifndef EWOMS_VCFV_PROPERTY_DEFAULTS_HH
#define EWOMS_VCFV_PROPERTY_DEFAULTS_HH

#include "vcfvproperties.hh"
#include "vcfvdiscretization.hh"
#include "vcfvstencil.hh"
#include "vcfvgradientcalculator.hh"
#include "vcfvgridcommhandlefactory.hh"
#include "vcfvvtkbaseoutputmodule.hh"

#include <ewoms/linear/vertexborderlistfromgrid.hh>
#include <ewoms/disc/common/fvbasepropertydefaults.hh>

namespace Ewoms {
template <class TypeTag>
class VcfvDiscretization;
}

namespace Opm {
namespace Properties {
//! Set the stencil
SET_PROP(VcfvDiscretization, Stencil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::ctype CoordScalar;

public:
    typedef Ewoms::VcfvStencil<CoordScalar, GridView> type;
};

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(VcfvDiscretization, DofMapper, typename GET_PROP_TYPE(TypeTag, VertexMapper));

//! The concrete class which manages the spatial discretization
SET_TYPE_PROP(VcfvDiscretization, Discretization, Ewoms::VcfvDiscretization<TypeTag>);

//! The base class for the VTK output modules (decides whether to write element or vertex based fields)
SET_TYPE_PROP(VcfvDiscretization, DiscVtkBaseOutputModule, Ewoms::VcfvVtkBaseOutputModule<TypeTag>);

//! Calculates the gradient of any quantity given the index of a flux approximation point
SET_TYPE_PROP(VcfvDiscretization, GradientCalculator, Ewoms::VcfvGradientCalculator<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(VcfvDiscretization, GridCommHandleFactory, Ewoms::VcfvGridCommHandleFactory<TypeTag>);

//! Use P1-finite element gradients by default for the vertex centered
//! finite volume scheme.
SET_BOOL_PROP(VcfvDiscretization, UseTwoPointGradients, false);

//! Set the border list creator for vertices
SET_PROP(VcfvDiscretization, BorderListCreator)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Ewoms::Linear::VertexBorderListFromGrid<GridView, VertexMapper> type;
};

}} // namespace Opm, Properties

#endif
