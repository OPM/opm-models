// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \ingroup EcfvDiscretization
 *
 * \brief Defines defaults for the common properties of the ECFV discretizations.
 */
#ifndef EWOMS_ECFV_PROPERTY_DEFAULTS_HH
#define EWOMS_ECFV_PROPERTY_DEFAULTS_HH

#include "ecfvproperties.hh"
#include "ecfvdiscretization.hh"
#include "ecfvstencil.hh"
#include "ecfvgridcommhandlefactory.hh"
#include "ecfvvtkoutputmodule.hh"

#include <ewoms/disc/common/fvbasepropertydefaults.hh>

namespace Ewoms {
template <class TypeTag>
class EcfvDiscretization;
}

namespace Opm {
namespace Properties {
//! Set the stencil
SET_PROP(EcfvDiscretization, Stencil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Ewoms::EcfvStencil<Scalar, GridView> type;
};

//! Mapper for the degrees of freedoms.
SET_TYPE_PROP(EcfvDiscretization, DofMapper, typename GET_PROP_TYPE(TypeTag, ElementMapper));

//! The concrete class which manages the spatial discretization
SET_TYPE_PROP(EcfvDiscretization, Discretization, Ewoms::EcfvDiscretization<TypeTag>);

//! The base class for the VTK output modules (decides whether to write element or vertex based fields)
SET_TYPE_PROP(EcfvDiscretization, DiscVtkOutputModule, Ewoms::EcfvVtkOutputModule<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(EcfvDiscretization, GridCommHandleFactory, Ewoms::EcfvGridCommHandleFactory<TypeTag>);

} // namespace Properties
} // namespace Opm

#endif
