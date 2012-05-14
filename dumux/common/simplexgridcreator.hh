// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief Provides a grid creator which a regular grid made of
 *        simplices.
 */
#ifndef DUMUX_SIMPLEX_GRID_CREATOR_HH
#define DUMUX_SIMPLEX_GRID_CREATOR_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);

NEW_PROP_TAG(GridSizeX);
NEW_PROP_TAG(GridSizeY);
NEW_PROP_TAG(GridSizeZ);

NEW_PROP_TAG(GridCellsX);
NEW_PROP_TAG(GridCellsY);
NEW_PROP_TAG(GridCellsZ);
}

/*!
 * \brief Provides a grid creator which a regular grid made of
 *        simplices.
 */
template <class TypeTag>
class SimplexGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid)  Grid;
    typedef Dune::shared_ptr<Grid> GridPointer;

    typedef typename Grid::ctype CoordScalar;
    enum { dimWorld = Grid::dimensionworld };
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    
public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::array<unsigned, dimWorld > cellRes;
        GlobalPosition upperRight;
        GlobalPosition lowerLeft;

        lowerLeft[0] = 0.0;
        upperRight[0] = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, SizeX);
        cellRes[0] = GET_PARAM_FROM_GROUP(TypeTag, int, Grid, CellsX);
        if (dimWorld > 1)
        {
            lowerLeft[1] = 0.0;
            upperRight[1] = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, SizeY);
            cellRes[1] = GET_PARAM_FROM_GROUP(TypeTag, int, Grid, CellsY);
        }
        if (dimWorld > 2)
        {
            lowerLeft[2] = 0.0;
            upperRight[2] = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, SizeZ);
            cellRes[2] = GET_PARAM_FROM_GROUP(TypeTag, int, Grid, CellsZ);
        }

        simplexGrid_ = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft, upperRight, cellRes);
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        return *simplexGrid_;
    };

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    {
        simplexGrid_->loadBalance();
    };

private:
    static GridPointer simplexGrid_;
};

template <class TypeTag>
typename SimplexGridCreator<TypeTag>::GridPointer SimplexGridCreator<TypeTag>::simplexGrid_;

}

#endif
