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

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/dofmanager.hh>
#endif

#include <type_traits>
#include <memory>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridManager);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(GridPart);
NEW_PROP_TAG(GridViewLevel);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(GridGlobalRefinements);
NEW_PROP_TAG(Simulator);
} // namespace Properties

//! \cond SKIP_THIS

/* Functions to retrieve the GridView regardless of whether the leaf
 * or a level grid view was chosen.  This code is a bit tricky as it
 * uses the SFINAE rule and also has to work around some C++
 * ideosyncracities...
 */
namespace BaseGridManagerHelper {

// leaf grid view retrieval methods
template <class TypeTag, class GridManager>
typename std::enable_if<std::is_same<typename GET_PROP_TYPE(TypeTag, GridView),
                                     typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView>::value,
                        typename GET_PROP_TYPE(TypeTag, GridView)>::type
gimmeGridView_(GridManager &gridManager)
{
    // return the leaf grid view
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    return gridManager.grid().leafGridView();
#else
    return gridManager.grid().leafView();
#endif
}

template <class TypeTag, class GridManager>
const typename std::enable_if<std::is_same<typename GET_PROP_TYPE(TypeTag, GridView),
                                           typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView>::value,
                              typename GET_PROP_TYPE(TypeTag, GridView)>::type
gimmeGridView_(const GridManager &gridManager)
{
    // return the leaf grid view
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    return gridManager.grid().leafGridView();
#else
    return gridManager.grid().leafView();
#endif
}

// level grid view retrieval methods
template <class TypeTag, class GridManager>
typename std::enable_if<std::is_same<typename GET_PROP_TYPE(TypeTag, GridView),
                                           typename GET_PROP_TYPE(TypeTag, Grid)::LevelGridView>::value,
                              typename GET_PROP_TYPE(TypeTag, GridView)>::type
gimmeGridView_(GridManager &gridManager)
{
    // return the chosen level grid view
    int level = GET_PROP_VALUE(TypeTag, GridViewLevel);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    return gridManager.grid().levelGridView(level);
#else
    return gridManager.grid().levelView(level);
#endif
}

template <class TypeTag, class GridManager>
const typename std::enable_if<std::is_same<typename GET_PROP_TYPE(TypeTag, GridView),
                                           typename GET_PROP_TYPE(TypeTag, Grid)::LevelGridView>::value,
                              typename GET_PROP_TYPE(TypeTag, GridView)>::type
gimmeGridView_(const GridManager &gridManager)
{
    // return the chosen level grid view
    int level = GET_PROP_VALUE(TypeTag, GridViewLevel);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    return gridManager.grid().levelGridView(level);
#else
    return gridManager.grid().levelView(level);
#endif
}
} // namespace BaseGridManagerHelper
//! \endcond

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

#if HAVE_DUNE_FEM
    typedef typename GET_PROP_TYPE(TypeTag, GridPart) GridPart;
#endif

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

#if HAVE_DUNE_FEM
    /*!
     * \brief Returns a reference to the grid part to be used.
     */
    const GridPart &gridPart() const
    { return *gridPart_; }

    /*!
     * \brief Returns a reference to the grid part to be used.
     */
    GridPart &gridPart()
    { return *gridPart_; }

    int gridSequenceNumber () const
    {
      typedef Dune::Fem::DofManager< Grid > FemDofManager;
      return FemDofManager::instance( gridPart().grid() ).sequence();
    }
#endif


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
#if HAVE_DUNE_FEM
        gridPart_.reset(new GridPart(asImp_().grid()));
        gridView_.reset(new GridView(gridPart_->gridView()));
#else
        gridView_.reset(new GridView(BaseGridManagerHelper::gimmeGridView_<TypeTag>(asImp_())));
#endif
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Simulator &simulator_;
    std::unique_ptr<GridView> gridView_;
#if HAVE_DUNE_FEM
    std::unique_ptr<GridPart> gridPart_;
#endif
};

} // namespace Ewoms

#endif
