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

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif
#include <dune/common/fvector.hh>

#include <vector>

namespace Dumux
{
// some hacky defines for the grid creator
#define LENS_DIM 2
#define LENS_CUBES 1

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
// declare the properties required by the for the lens grid creator
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(GridSizeX);
NEW_PROP_TAG(GridSizeY);
NEW_PROP_TAG(GridSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);
}

/*!
 * \brief Helper class for grid instantiation of the lens problem.
 */
#if HAVE_UG
template <class TypeTag>
class LensGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { dim = LENS_DIM };

public:
    typedef Dune::UGGrid<LENS_DIM> Grid;

    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        grid_ = new Grid;
        
        Dune::FieldVector<int, dim> cellRes;
        Dune::FieldVector<Scalar, dim> upperRight;
        Dune::FieldVector<Scalar, dim> lowerLeft;

        lowerLeft = 0.0;
        upperRight[0] = GET_PARAM(TypeTag, Scalar, GridSizeX);
        upperRight[1] = GET_PARAM(TypeTag, Scalar, GridSizeY);
        
        cellRes[0] = GET_PARAM(TypeTag, int, CellsX);
        cellRes[1] = GET_PARAM(TypeTag, int, CellsY);
        if (dim == 3) {
            upperRight[2] = GET_PARAM(TypeTag, Scalar, GridSizeZ);
            cellRes[2] = GET_PARAM(TypeTag, int, CellsZ);
        }

        Dune::GridFactory<Dune::UGGrid<dim> > factory(grid_);
        
        if (dim == 3) {
            Dune::FieldVector<double,dim> pos;
            for (int k = 0; k <= cellRes[0]; k++) {
                pos[2] = upperRight[2]*double(k)/cellRes[2];
                
                for (int j = 0; j <= cellRes[1]; j++) {
                    pos[1] = upperRight[1]*double(j)/cellRes[1];

                    for (int i = 0; i <= cellRes[0]; i++) {
                        pos[0] = upperRight[0]*double(i)/cellRes[0];
                        factory.insertVertex(pos);
                    }
                }
            }
        }
        else {
            assert(dim == 2);
            Dune::FieldVector<double,dim> pos;
            for (int j = 0; j <= cellRes[1]; j++) {
                pos[1] = upperRight[1]*double(j)/cellRes[1];
                
                for (int i = 0; i <= cellRes[0]; i++) {
                    pos[0] = upperRight[0]*double(i)/cellRes[0];
                    factory.insertVertex(pos);
                }
            }
        }

        for (int i = 0; i < cellRes[0]; ++i) {
            for (int j = 0; j < cellRes[1]; ++j) {
#if LENS_CUBES
                std::vector<unsigned int> v(1 << dim);
#else
                std::vector<unsigned int> v(dim + 1);
#endif
                if (dim == 3) {
                    int m = cellRes[0] + 1;
                    int n = cellRes[1] + 1;
                    for (int k = 0; k < cellRes[2]; ++k) {
                        int i0 = k*m*n + j*m + i;
                        int i1 = k*m*n + j*m + (i+1);
                        int i2 = k*m*n + (j+1)*m + i;
                        int i3 = k*m*n + (j+1)*m + (i+1);
                        int i4 = (k+1)*m*n + j*m + i;
                        int i5 = (k+1)*m*n + j*m + (i+1);
                        int i6 = (k+1)*m*n + (j+1)*m + i;
                        int i7 = (k+1)*m*n + (j+1)*m + (i+1);
                        
#if LENS_CUBES
                        v[0] = i0;
                        v[1] = i1;
                        v[2] = i2;
                        v[3] = i3;
                        v[4] = i4;
                        v[5] = i5;
                        v[6] = i6;
                        v[7] = i7;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,3), v);

#else
                        v[0] = i0;
                        v[1] = i1;
                        v[2] = i2;
                        v[3] = i4;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), v);
                        
                        v[0] = i4;
                        v[1] = i5;
                        v[2] = i6;
                        v[3] = i2;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), v);

                        v[0] = i2;
                        v[1] = i5;
                        v[2] = i4;
                        v[3] = i1;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), v);

                        v[0] = i2;
                        v[1] = i3;
                        v[2] = i7;
                        v[3] = i5;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), v);

                        v[0] = i5;
                        v[1] = i7;
                        v[2] = i6;
                        v[3] = i2;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), v);

                        v[0] = i1;
                        v[1] = i3;
                        v[2] = i5;
                        v[3] = i2;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), v);
#endif
                    }
                }
                else {
                    assert(dim == 2);
                    
                    int m = cellRes[0] + 1;
                    int i0 = j*m + i;
                    int i1 = j*m + (i+1);
                    int i2 = (j+1)*m + i;
                    int i3 = (j+1)*m + (i+1);
#if LENS_CUBES
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
                    v[1] = i3;
                    v[2] = i2;
                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), v);
#endif
                }
            }
        }

        grid_ = factory.createGrid();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    { return *grid_; }

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    static void loadBalance()
    { grid_->loadBalance(); }

private:
    static Grid *grid_;
};

template <class TypeTag>
Dune::UGGrid<LENS_DIM> *LensGridCreator<TypeTag>::grid_;

#else // ! HAVE_UG

template <class TypeTag>
class LensGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { dim = LENS_DIM };

public:
    typedef Dune::YaspGrid<LENS_DIM> Grid;

    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        Dune::FieldVector<int, LENS_DIM> cellRes;
        Dune::FieldVector<Scalar, LENS_DIM> upperRight;
        Dune::FieldVector<Scalar, LENS_DIM> lowerLeft;

        lowerLeft[1] = 0.0;
        upperRight[0] = GET_PARAM(TypeTag, Scalar, GridSizeX);
        upperRight[1] = GET_PARAM(TypeTag, Scalar, GridSizeY);
        
        cellRes[0] = GET_PARAM(TypeTag, int, CellsX);
        cellRes[1] = GET_PARAM(TypeTag, int, CellsY);
        if (dim == 3) {
            upperRight[2] = GET_PARAM(TypeTag, Scalar, GridSizeZ);
            cellRes[2] = GET_PARAM(TypeTag, int, CellsZ);
        }
        
        grid_ = new Dune::YaspGrid<LENS_DIM>(
#ifdef HAVE_MPI
            Dune::MPIHelper::getCommunicator(),
#endif
            upperRight, // upper right
            cellRes, // number of cells
            Dune::FieldVector<bool,LENS_DIM>(false), // periodic
            0); // overlap
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    { return *grid_; }

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    static void loadBalance()
    { grid_->loadBalance(); }

private:
    static Grid *grid_;
};

template <class TypeTag>
Dune::YaspGrid<LENS_DIM> *LensGridCreator<TypeTag>::grid_;

#endif // HAVE_UG

}

#endif
