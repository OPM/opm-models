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
 * \copydoc Ewoms::VcfvDiscretization
 */
#ifndef EWOMS_VCFV_DISCRETIZATION_HH
#define EWOMS_VCFV_DISCRETIZATION_HH

#include "vcfvproperties.hh"
#include "vcfvstencil.hh"
#include "vcfvgradientcalculator.hh"
#include "vcfvgridcommhandlefactory.hh"
#include "vcfvbaseoutputmodule.hh"

#include <ewoms/linear/vertexborderlistfromgrid.hh>
#include <ewoms/disc/common/fvbasediscretization.hh>

namespace Ewoms {
template <class TypeTag>
class VcfvDiscretization;
}

namespace Ewoms {
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

//! The base class for the output modules (decides whether to write
//! element or vertex based fields)
SET_TYPE_PROP(VcfvDiscretization, DiscBaseOutputModule,
              Ewoms::VcfvBaseOutputModule<TypeTag>);

//! Calculates the gradient of any quantity given the index of a flux approximation point
SET_TYPE_PROP(VcfvDiscretization, GradientCalculator,
              Ewoms::VcfvGradientCalculator<TypeTag>);

//! The class to create grid communication handles
SET_TYPE_PROP(VcfvDiscretization, GridCommHandleFactory,
              Ewoms::VcfvGridCommHandleFactory<TypeTag>);

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

//! For the vertex centered finite volume method, ghost and overlap elements must _not_
//! be assembled to avoid accounting twice for the fluxes over the process boundary faces
//! of the local process' grid partition
SET_BOOL_PROP(VcfvDiscretization, LinearizeNonLocalElements, false);

}} // namespace Ewoms, Properties

namespace Ewoms {
/*!
 * \ingroup Discretization
 *
 * \brief The base class for the vertex centered finite volume discretization scheme.
 */
template<class TypeTag>
class VcfvDiscretization : public FvBaseDiscretization<TypeTag>
{
    typedef FvBaseDiscretization<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, DofMapper) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    enum { dim = GridView::dimension };


public:
    VcfvDiscretization(Simulator &simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Returns a string of discretization's human-readable name
     */
    static std::string discretizationName()
    { return "vcfv"; }

    /*!
     * \brief Returns the number of global degrees of freedom (DOFs) due to the grid
     */
    size_t numGridDof() const
    { return this->gridView_.size(/*codim=*/dim); }

    /*!
     * \brief Mapper to convert the Dune entities of the
     *        discretization's degrees of freedoms are to indices.
     */
    const DofMapper &dofMapper() const
    { return this->simulator_.problem().vertexMapper(); }

    /*!
     * \brief Serializes the current state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    { res.template serializeEntities</*codim=*/dim>(asImp_(), this->gridView_); }

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
        res.template deserializeEntities</*codim=*/dim>(asImp_(), this->gridView_);
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
