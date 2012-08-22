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
#include "config.h"

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

#include<dune/common/version.hh>
#include<dune/geometry/quadraturerules.hh>

#include <dumux/boxmodels/common/boxfvelementgeometry.hh>

const int dim = 3;
typedef double Scalar;
typedef Dumux::QuadrialteralQuadratureGeometry<Scalar, dim> QuadratureGeom;
typedef QuadratureGeom::LocalPosition LocalPosition;
typedef QuadratureGeom::GlobalPosition GlobalPosition;

GlobalPosition::field_type f(const GlobalPosition &pos)
{
    GlobalPosition::field_type result = 1;
    for (int i = 0; i < GlobalPosition::dimension; ++i)
        result *= pos[i];
    return result;
}

void testIdenityMapping()
{
    QuadratureGeom foo;
    
    Scalar corners[][3] = {
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 1, 1, 0 },
        { 0, 0, 1 },
        { 1, 0, 1 },
        { 0, 1, 1 },
        { 1, 1, 1 }
    };
    foo.setCorners(corners, 8);

    std::cout << "testing identity mapping...\n";
    int n = 100;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                LocalPosition localPos;
                
                localPos[0] = Scalar(i)/(n - 1);
                localPos[1] = Scalar(j)/(n - 1);
                localPos[2] = Scalar(k)/(n - 1);
                
                GlobalPosition globalPos = foo.global(localPos);
                
                GlobalPosition diff(localPos);
                diff -= globalPos;
                assert(diff.two_norm() < 1e-10);
            }
        }
    }
}

template <class Grid>
void writeSubControlVolumes(const Grid &grid)
{
#if HAVE_ALUGRID
    typedef typename Grid::LeafGridView GridView;

    typedef Dune::ALUCubeGrid<dim, dim> Grid2;
    typedef typename Grid2::LeafGridView GridView2;
    typedef Dune::GridFactory<Grid2> GridFactory2;
    
    // instanciate a FVElementGeometry
    typedef Dumux::BoxFVElementGeometry<Scalar, GridView> FVElementGeometry;
    FVElementGeometry fvElemGeom;

    GridFactory2 gf2;
    const auto &gridView = grid.leafView();
    auto eIt = gridView.template begin<0>();
    const auto &eEndIt = gridView.template end<0>();
    for (; eIt != eEndIt; ++eIt) {
        fvElemGeom.update(gridView, *eIt);
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {      
            const auto &scvLocalGeom = *(fvElemGeom.subContVol[scvIdx].localGeometry);

            for (int i = 0; i < scvLocalGeom.numCorners; ++ i) {
                GlobalPosition pos(eIt->geometry().global(scvLocalGeom.corner(i)));
                gf2.insertVertex(pos);
            }
        }
    }

    int cornerOffset = 0;
    eIt = gridView.template begin<0>();
    for (; eIt != eEndIt; ++eIt) {
        fvElemGeom.update(gridView, *eIt);
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {      
            const auto &scvLocalGeom = *fvElemGeom.subContVol[scvIdx].localGeometry;
            
            std::vector<unsigned int> vertexIndices;
            for (int i = 0; i < scvLocalGeom.numCorners; ++ i) {
                vertexIndices.push_back(cornerOffset);
                ++ cornerOffset;
            }

            gf2.insertElement(Dune::GeometryType(Dune::GeometryType::cube,dim),
                              vertexIndices);
        }
    }

    const auto &grid2 = *gf2.createGrid();
    typedef Dune::VTKWriter<GridView2> VtkWriter;
    VtkWriter writer(grid2.leafView(), Dune::VTK::conforming);
    writer.write("quadrature", Dune::VTK::ascii);
#endif // HAVE_ALUGRID
}

void testTetrahedron()
{
#if HAVE_ALUGRID
    typedef Dune::ALUSimplexGrid<dim, dim> Grid;
    typedef Grid::LeafGridView GridView;
    typedef Dune::GridFactory<Grid> GridFactory;
    GridFactory gf;
    Scalar corners[][3] = {
        { 25, 25, 25 },
        { 50, 25, 25 },
        { 25, 75, 25 },
        { 25, 25, 50 }
    };

    for (unsigned i = 0; i < sizeof(corners)/sizeof(corners[0]); ++i) {
        GlobalPosition pos;
        for (unsigned j = 0; j < dim; ++j)
            pos[j] = corners[i][j];
        gf.insertVertex(pos);
    }
    std::vector<unsigned int> v = { 0, 1, 2, 3 };
    gf.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,dim),
                     v);
    const auto &grid = *gf.createGrid();

    writeSubControlVolumes(grid);
#endif // HAVE_ALUGRID
}

void testQuadrature()
{
    std::cout << "testing SCV quadrature...\n";

    GlobalPosition upperRight(1.0);
    Dune::FieldVector<int, dim> cellRes(10);

    typedef Dune::YaspGrid<dim> Grid;
    typedef Grid::LeafGridView GridView;
    Grid grid(
#ifdef HAVE_MPI
        Dune::MPIHelper::getCommunicator(),
#endif
        upperRight, // upper right
        cellRes, // number of cells
        Dune::FieldVector<bool, dim>(false), // periodic
        0); // overlap

    // compute approximate integral
    auto gridView = grid.leafView();
    auto eIt = gridView.begin<0>();
    const auto eEndIt = gridView.end<0>();
    Scalar result=0;
    // instanciate a FVElementGeometry
    typedef Dumux::BoxFVElementGeometry<Scalar, GridView> FVElementGeometry;
    FVElementGeometry fvElemGeom;
    for (; eIt != eEndIt; ++eIt) {
        const auto &elemGeom = eIt->geometry();

        fvElemGeom.update(gridView, *eIt);

        // loop over all sub-control volumes
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {
            const auto &scvLocalGeom = *fvElemGeom.subContVol[scvIdx].localGeometry;               
            
            Dune::GeometryType geomType = scvLocalGeom.type();
            static const int quadratureOrder = 2;
            const auto &rule = Dune::QuadratureRules<Scalar,dim>::rule(geomType, quadratureOrder);
            
            // integrate f over the sub-control volume
            for (auto it = rule.begin(); it != rule.end(); ++ it)
            {
                auto posScvLocal = it->position();
                auto posElemLocal = scvLocalGeom.global(posScvLocal);
                auto posGlobal = elemGeom.global(posScvLocal);
            
                Scalar fval = f(posGlobal);
                Scalar weight = it->weight();
                Scalar detjac = 
                    scvLocalGeom.integrationElement(posScvLocal) * 
                    elemGeom.integrationElement(posElemLocal);
                
                result += fval * weight * detjac;
            }
        }
    }
    
    std::cout << "result: " << result << " (expected value: " << 1.0/(1 << dim) << ")\n";
}
    
int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    //testIdenityMapping();
    //testTetrahedron();
    testQuadrature();

    return 0;
}
