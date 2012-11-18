// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Markus Wolff                                      *
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
 * \copydoc Ewoms::CubeGridCreator
 */
#ifndef EWOMS_CUBE_GRID_CREATOR_HH
#define EWOMS_CUBE_GRID_CREATOR_HH

#include <ewoms/common/basicproperties.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parameters.hh>

#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/fvector.hh>

namespace Ewoms
{

namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);
}

/*!
 * \brief Provides a grid creator which a regular grid made of
 *        quadrilaterals.
 *
 * A quadirlateral is a line segment in 1D, a rectangle in 2D and a
 * cube in 3D.
 */
template <class TypeTag>
class CubeGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Grid)  Grid;
    typedef Dune::shared_ptr<Grid> GridPointer;

    typedef typename Grid::ctype CoordScalar;
    enum { dimWorld = Grid::dimensionworld };
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Register all run-time parameters for the grid creator.
     */
    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, int, GridGlobalRefinements, "The number of global refinements of the grid executed after it was loaded");
        REGISTER_PARAM(TypeTag, Scalar, DomainSizeX, "The size of the domain in x direction");
        REGISTER_PARAM(TypeTag, Scalar, CellsX, "The number of intervalls in x direction");
        if (dimWorld > 1) {
            REGISTER_PARAM(TypeTag, Scalar, DomainSizeY, "The size of the domain in y direction");
            REGISTER_PARAM(TypeTag, Scalar, CellsY, "The number of intervalls in y direction");
        }
        if (dimWorld > 2) {
            REGISTER_PARAM(TypeTag, Scalar, DomainSizeZ, "The size of the domain in z direction");
            REGISTER_PARAM(TypeTag, Scalar, CellsZ, "The number of intervalls in z direction");
        }
    }

    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::array< unsigned int, dimWorld > cellRes;
        GlobalPosition upperRight(0.0);
        GlobalPosition lowerLeft(0.0);

        upperRight[0] = GET_PARAM(TypeTag, Scalar, DomainSizeX);
        cellRes[0] = GET_PARAM(TypeTag, int, CellsX);
        if (dimWorld > 1) {
            upperRight[1] = GET_PARAM(TypeTag, Scalar, DomainSizeY);
            cellRes[1] = GET_PARAM(TypeTag, int, CellsY);
        }
        if (dimWorld > 2)
        {
            upperRight[2] = GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = GET_PARAM(TypeTag, int, CellsZ);
        }
        unsigned numRefinments = GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);

        cubeGrid_ = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, cellRes);
        cubeGrid_->globalRefine(numRefinments);
        initialized_ = true;
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    { return *cubeGrid_; }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    { cubeGrid_->loadBalance(); }
    /*!
     * \brief Destroys the grid
     *
     * This is required to guarantee that the grid is deleted before MPI_Comm_free is called.
     */
    static void deleteGrid()
    { if (initialized_) cubeGrid_.reset(); initialized_ = false; }


protected:
    static GridPointer cubeGrid_;
    static bool initialized_;
};

template <class TypeTag>
typename Ewoms::CubeGridCreator<TypeTag>::GridPointer CubeGridCreator<TypeTag>::cubeGrid_;

template <class TypeTag>
bool CubeGridCreator<TypeTag>::initialized_ = false;

}

#endif
