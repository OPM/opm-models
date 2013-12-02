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
 *
 * \copydoc Ewoms::VcfvConstraintsContext
 */
#ifndef EWOMS_VCFV_CONSTRAINTS_CONTEXT_HH
#define EWOMS_VCFV_CONSTRAINTS_CONTEXT_HH

#include "vcfvproperties.hh"
#include "vcfvelementcontext.hh"

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup VcfvModel
 *
 * \brief Represents all quantities which available for calculating constraints
 */
template <class TypeTag>
class VcfvConstraintsContext
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dimWorld = GridView::dimensionworld };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor.
     */
    explicit VcfvConstraintsContext(const ElementContext &elemCtx)
        : elemCtx_(elemCtx)
    {}

    /*!
     * \copydoc Ewoms::VcfvElementContext::problem()
     */
    const Problem &problem() const
    { return elemCtx_.problem(); }

    /*!
     * \copydoc Ewoms::VcfvElementContext::model()
     */
    const Model &model() const
    { return elemCtx_.model(); }

    /*!
     * \copydoc Ewoms::VcfvElementContext::gridView()
     */
    const GridView &gridView() const
    { return elemCtx_.gridView(); }

    /*!
     * \copydoc Ewoms::VcfvElementContext::element()
     */
    const Element &element() const
    { return elemCtx_.element(); }

    /*!
     * \copydoc Ewoms::VcfvElementContext::numScv()
     */
    int numScv() const
    { return elemCtx_.numScv(); }

    /*!
     * \copydoc Ewoms::VcfvElementContext::numScvf()
     */
    int numScvf() const
    { return elemCtx_.numScvf(); }

    /*!
     * \copydoc Ewoms::VcfvElementContext::globalSpaceIndex
     */
    int globalSpaceIndex(int scvIdx, int timeIdx) const
    { return elemCtx_.globalSpaceIndex(scvIdx, timeIdx); }

    /*!
     * \copydoc Ewoms::VcfvElementContext::pos
     */
    GlobalPosition pos(int scvIdx, int timeIdx) const
    { return elemCtx_.pos(scvIdx, timeIdx); }

protected:
    const ElementContext &elemCtx_;
};

} // namespace Ewoms

#endif
