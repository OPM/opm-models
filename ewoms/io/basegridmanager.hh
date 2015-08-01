// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2014 by Andreas Lauser

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
 * \copydoc Ewoms::BaseGridManager
 */
#ifndef EWOMS_BASE_GRID_MANAGER_HH
#define EWOMS_BASE_GRID_MANAGER_HH

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/common/version.hh>

#include <type_traits>
#include <memory>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridManager);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(GridViewLevel);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(GridGlobalRefinements);
NEW_PROP_TAG(Simulator);
} // namespace Properties

/*!
 * \brief Provides the base class for most (all?) grid managers.
 */
template <class TypeTag>
class BaseGridManager
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, GridManager) Implementation;

public:
    BaseGridManager(Simulator &simulator)
        : simulator_(simulator)
    {}

    /*!
     * \brief Returns a reference to the grid view to be used.
     */
    GridView &gridView()
    { return *gridView_; }

    /*!
     * \brief Returns a reference to the grid view to be used.
     */
    const GridView &gridView() const
    { return *gridView_; }

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    void loadBalance()
    { asImp_().grid().loadBalance(); }

protected:
    // this method should be called after the grid has been allocated
    void finalizeInit_()
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,3)
        gridView_.reset(new GridView(asImp_().grid().leafGridView()));
#else
        gridView_.reset(new GridView(asImp_().grid().leafView()));
#endif
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Simulator &simulator_;
    std::unique_ptr<GridView> gridView_;
};

} // namespace Ewoms

#endif
