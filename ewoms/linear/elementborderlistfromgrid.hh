// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2013 by Andreas Lauser                               *
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
/*!
 * \file
 * \copydoc Ewoms::Linear::ElementBorderListFromGrid
 */
#ifndef EWOMS_ELEMENT_BORDER_LIST_FROM_GRID_HH
#define EWOMS_ELEMENT_BORDER_LIST_FROM_GRID_HH

#include "overlaptypes.hh"

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include <algorithm>

namespace Ewoms {
namespace Linear {
/*!
 * \brief Uses communication on the grid to find the initial seed list
 *        of indices for methods which use element-based degrees of
 *        freedom.
 */
template <class GridView, class ElementMapper>
class ElementBorderListFromGrid : public Dune::CommDataHandleIF<ElementBorderListFromGrid<GridView, ElementMapper>,
                                                               int >
{
public:
    ElementBorderListFromGrid(const GridView &gridView,
                              const ElementMapper &map)
        : gridView_(gridView), map_(map)
    {
        gridView.communicate(*this,
                             Dune::InteriorBorder_All_Interface,
                             Dune::BackwardCommunication);
    }

    // data handle methods
    bool contains (int dim, int codim) const
    { return codim == 0; }

    bool fixedsize(int dim, int codim) const
    { return true; }

    template<class EntityType>
    size_t size(const EntityType &e) const
    { return 2; }

    template<class MessageBufferImp, class EntityType>
    void gather(MessageBufferImp &buff, const EntityType &e) const
    {
        buff.write(gridView_.comm().rank());
        buff.write(map_.map(e));
    }

    template<class MessageBufferImp, class EntityType>
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        BorderIndex bIdx;

        bIdx.localIdx = map_.map(e);
        { int tmp; buff.read(tmp); bIdx.peerRank = tmp; }
        { int tmp; buff.read(tmp); bIdx.peerIdx = tmp; }
        bIdx.borderDistance = 1;

        borderList_.push_back(bIdx);
    }

    // Access to the border list.
    const BorderList &borderList() const
    { return borderList_; }

private:
    const GridView gridView_;
    const ElementMapper &map_;
    BorderList borderList_;
};

} // namespace Linear
} // namespace Ewoms

#endif
