// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2008-2013 by Andreas Lauser

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
 * \copydoc Ewoms::EcfvDiscretization
 */
#ifndef EWOMS_ECFV_DISCRETIZATION_HH
#define EWOMS_ECFV_DISCRETIZATION_HH

#include "ecfvproperties.hh"
#include "ecfvstencil.hh"
#include "ecfvgridcommhandlefactory.hh"
#include "ecfvvtkbaseoutputmodule.hh"

#include <ewoms/linear/elementborderlistfromgrid.hh>
#include <ewoms/disc/common/fvbasediscretization.hh>

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
SET_TYPE_PROP(EcfvDiscretization, DiscVtkBaseOutputModule, Ewoms::EcfvVtkBaseOutputModule<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(EcfvDiscretization, GridCommHandleFactory, Ewoms::EcfvGridCommHandleFactory<TypeTag>);

//! Set the border list creator for vertices
SET_PROP(EcfvDiscretization, BorderListCreator)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef Ewoms::Linear::ElementBorderListFromGrid<GridView, ElementMapper> type;
};

} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief The base class for the element centered finite volume discretization scheme.
 */
template<class TypeTag>
class EcfvDiscretization : public FvBaseDiscretization<TypeTag>
{
    typedef FvBaseDiscretization<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum { dim = GridView::dimension };

    typedef typename GridView::template Codim<dim>::Entity Vertex;

public:
    EcfvDiscretization(Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    std::string discretizationName() const
    { return "ecfv"; }

    /*!
     * \brief Returns the number of global degrees of freedoms (DOFs)
     */
    size_t numDof() const
    { return this->gridView_.size(/*codim=*/0); }

    /*!
     * \brief Mapper to convert the Dune entities of the
     *        discretization's degrees of freedoms are to indices.
     */
    const DofMapper &dofMapper() const
    { return this->problem_.elementMapper(); }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    { res.template serializeEntities</*codim=*/0>(asImp_(), this->gridView_); }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.template deserializeEntities</*codim=*/0>(asImp_(), this->gridView_);
        this->solution_[/*timeIdx=*/1] = this->solution_[/*timeIdx=*/0];
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
} // namespace Ewoms

#endif
