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
 * \brief Grid creator for the 2pni injection problem used for testing
 *        the fully implicit non-isothermal two-phase model.
 */
#ifndef DUMUX_INJECTION_2PNI_GRID_CREATOR_HH
#define DUMUX_INJECTION_2PNI_GRID_CREATOR_HH

#include "injectionproblem2pni.hh"

namespace Dumux
{
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
template <class TypeTag>
class Injection2PNIGridCreator
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
    {
        return *grid_;
    };

private:
    static Grid *grid_;
};

template <class TypeTag>
Dune::YaspGrid<2> *Injection2PNIGridCreator<TypeTag>::grid_;
}

#endif
