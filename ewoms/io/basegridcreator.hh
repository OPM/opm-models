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
 * \copydoc Ewoms::BaseGridCreator
 */
#ifndef EWOMS_BASE_GRID_CREATOR_HH
#define EWOMS_BASE_GRID_CREATOR_HH

#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <type_traits>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(GridViewLevel);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(GridGlobalRefinements);
}} // namespace Opm, Properties

namespace Ewoms {
template <class TypeTag>
typename std::enable_if<std::is_same<typename GET_PROP_TYPE(TypeTag, GridView),
                                     typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView>::value,
                        typename GET_PROP_TYPE(TypeTag, GridView)>::type
gimmeGridView_()
{
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    // return the leaf grid view
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    return GridCreator::grid().leafGridView();
#else
    return GridCreator::grid().leafView();
#endif
}

template <class TypeTag>
typename std::enable_if<std::is_same<typename GET_PROP_TYPE(TypeTag, GridView),
                                     typename GET_PROP_TYPE(TypeTag, Grid)::LevelGridView>::value,
                        typename GET_PROP_TYPE(TypeTag, GridView)>::type
gimmeGridView_()
{
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    // return the chosen level grid view
    int level = GET_PROP_VALUE(TypeTag, GridViewLevel);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    return GridCreator::grid().levelGridView(level);
#else
    return GridCreator::grid().levelView(level);
#endif
}

/*!
 * \brief Provides a grid creator which reads Dune Grid Format (DGF) files
 */
template <class TypeTag>
class BaseGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) Implementation;

public:
    /*!
     * \brief Returns a reference to the grid view to be used.
     *
     * This code is a bit tricky as it uses the SFINAE rule and also
     * has to work around some C++ ideosyncracities...
     */
    static GridView gridView()
    { return gimmeGridView_<TypeTag>(); }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    { return *Implementation::gridPointer(); }

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    static void loadBalance()
    { Implementation::grid().loadBalance(); }
};

} // namespace Ewoms

#endif
