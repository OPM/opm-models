// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \brief Grid creator for the lens problem used for testing the fully
 *        implicit two-phase model.
 */
#ifndef DUMUX_LENS_GRID_CREATOR_HH
#define DUMUX_LENS_GRID_CREATOR_HH

#include "lensproblem.hh"

namespace Dumux
{

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
// declare the properties required by the for the lens grid creator
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(GridSizeX);
NEW_PROP_TAG(GridSizeY);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
}

/*!
 * \brief Helper class for grid instantiation of the lens problem.
 */
template <class TypeTag, class Grid = typename GET_PROP_TYPE(TypeTag, Grid)>
class LensGridCreator;

#if HAVE_UG
template <class TypeTag>
class LensGridCreator<TypeTag, Dune::UGGrid<2> >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::UGGrid<2> Grid;

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::FieldVector<int, 2> cellRes;
        Dune::FieldVector<Scalar, 2> upperRight;
        Dune::FieldVector<Scalar, 2> lowerLeft;

        lowerLeft[0] = 0.0;
        lowerLeft[1] = 0.0;
        upperRight[0] = GET_PARAM(TypeTag, Scalar, GridSizeX);
        upperRight[1] = GET_PARAM(TypeTag, Scalar, GridSizeY);

        cellRes[0] = GET_PARAM(TypeTag, int, CellsX);
        cellRes[1] = GET_PARAM(TypeTag, int, CellsY);
        
        Dune::GridFactory<Dune::UGGrid<2> > factory;
        for (int i=0; i<=cellRes[0]; i++) {
            for (int j=0; j<=cellRes[1]; j++) {
                Dune::FieldVector<double,2> pos;
                pos[0] = upperRight[0]*double(i)/cellRes[0];
                pos[1] = upperRight[1]*double(j)/cellRes[1];
                factory.insertVertex(pos);
            }
        }

        for (int i=0; i<cellRes[0]; i++) {
            for (int j=0; j<cellRes[1]; j++) {
#if CUBES
                std::vector<unsigned int> v(4);
#else
                std::vector<unsigned int> v(3);
#endif

                int i0 = i*(cellRes[1]+1) + j;
                int i1 = i*(cellRes[1]+1) + j+1;
                int i2 = (i+1)*(cellRes[1]+1) + j;
                int i3 = (i+1)*(cellRes[1]+1) + j+1;

#if CUBES
                v[0] = i0;
                v[1] = i1;
                v[2] = i2;
                v[3] = i3;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
#else
                v[0] = i0;
                v[1] = i1;
                v[2] = i2;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), v);

                v[0] = i1;
                v[1] = i2;
                v[2] = i3;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), v);
#endif
            }
        }

        grid_ = factory.createGrid();
        grid_->loadBalance();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    { return *grid_; };

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    static void loadBalance()
    { grid_->loadBalance(); };
private:
    static Grid *grid_;
};

template <class TypeTag>
Dune::UGGrid<2> *LensGridCreator<TypeTag, Dune::UGGrid<2> >::grid_;
#endif

template <class TypeTag>
class LensGridCreator<TypeTag, Dune::YaspGrid<2> >
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::YaspGrid<2> Grid;

public:
    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::FieldVector<int, 2> cellRes;
        Dune::FieldVector<Scalar, 2> upperRight;
        Dune::FieldVector<Scalar, 2> lowerLeft;

        lowerLeft[0] = 0.0;
        lowerLeft[1] = 0.0;
        upperRight[0] = 6.0;
        upperRight[1] = 4.0;

        cellRes[0] = GET_PARAM(TypeTag, int, CellsX);
        cellRes[1] = GET_PARAM(TypeTag, int, CellsY);

        grid_ = new Dune::YaspGrid<2>(
#ifdef HAVE_MPI
            Dune::MPIHelper::getCommunicator(),
#endif
            upperRight, // upper right
            cellRes, // number of cells
            Dune::FieldVector<bool,2>(false), // periodic
            0); // overlap
    };

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    { return *grid_; };

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    static void loadBalance()
    { grid_->loadBalance(); };

private:
    static Grid *grid_;
};

template <class TypeTag>
Dune::YaspGrid<2> *LensGridCreator<TypeTag, Dune::YaspGrid<2> >::grid_;
}

#endif
