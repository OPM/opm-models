// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012-2013 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
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
 * \copydoc Ewoms::ArtGridCreator
 */
#ifndef EWOMS_ART_GRID_CREATOR_HH
#define EWOMS_ART_GRID_CREATOR_HH

#include <ewoms/models/discretefracture/fracturemapper.hh>
#include <ewoms/common/parametersystem.hh>
#include <opm/core/utility/PropertySystem.hpp>
#include <opm/material/Valgrind.hpp>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>
#include <memory>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(EnableFractures);
}
}

namespace Ewoms {
/*!
 * \brief Reads in mesh files in the ART format.
 *
 * This file format is used to specify grids with fractures.
 */
template <class TypeTag>
class ArtGridCreator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<double, 2> GlobalPosition;
    typedef Ewoms::FractureMapper<TypeTag> FractureMapper;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef std::shared_ptr<Grid> GridPointer;

public:

    /*!
     * \brief Register all run-time parameters for the grid creator.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, GridFile, "The file name of the DGF file to load");
    }

    /*!
     * \brief Create the Grid
     */
    static void makeGrid()
    {
        enum ParseMode { Vertex, Edge, Element, Finished };
        const std::string artFileName = EWOMS_GET_PARAM(TypeTag, std::string, GridFile);
        std::vector<GlobalPosition> vertexPos;
        std::vector<std::pair<int, int> > edges;
        std::vector<std::pair<int, int> > fractureEdges;
        Dune::GridFactory<Grid> gridFactory;
        std::ifstream inStream(artFileName);
        if (!inStream.is_open()) {
            DUNE_THROW(Dune::IOError,
                       "File '" << artFileName << "' does not exist or is not readable");
        }
        std::string curLine;
        ParseMode curParseMode = Vertex;
        while (inStream) {
            std::getline(inStream, curLine);

            // remove comments
            auto commentPos = curLine.find("%");
            if (commentPos != curLine.npos) {
                curLine = curLine.substr(0, commentPos);
            }

            // remove leading whitespace
            unsigned numLeadingSpaces = 0;
            while (curLine.size() > numLeadingSpaces && std::isspace(curLine[numLeadingSpaces]))
                ++ numLeadingSpaces;
            curLine = curLine.substr(numLeadingSpaces,
                                     curLine.size() - numLeadingSpaces);

            // remove trailing whitespace
            unsigned numTrailingSpaces = 0;
            while (curLine.size() > numTrailingSpaces && std::isspace(curLine[curLine.size() - numTrailingSpaces]))
                ++ numTrailingSpaces;
            curLine = curLine.substr(0, curLine.size() - numTrailingSpaces);

            // a section of the file is finished, go to the next one
            if (curLine == "$") {
                if (curParseMode == Vertex)
                    curParseMode = Edge;
                else if (curParseMode == Edge)
                    curParseMode = Element;
                else if (curParseMode == Element)
                    curParseMode = Finished;
                continue;
            }

            // skip empty lines
            if (curLine.empty())
                continue;

            if (curParseMode == Vertex) {
                GlobalPosition coord;
                std::istringstream iss(curLine);
                // parse only the first two numbers as the vertex
                // coordinate. the last number is the Z coordinate
                // which we ignore (so far)
                iss >> coord[0] >> coord[1];
                gridFactory.insertVertex(coord);
                vertexPos.push_back(coord);
            }
            else if (curParseMode == Edge) {
                // read an edge and update the fracture mapper

                // read the data attached to the edge
                std::istringstream iss(curLine);
                int dataVal;
                std::string tmp;
                iss >> dataVal;
                iss >> tmp;
                assert(tmp == ":");

                // read the vertex indices of an edge
                std::vector<unsigned int> vertIndices;
                while (iss) {
                    unsigned int tmp;
                    iss >> tmp;
                    if (!iss)
                        break;
                    vertIndices.push_back(tmp);
                    assert(tmp < vertexPos.size());
                }
                assert(vertIndices.size() == 2); // an edge always has two indices

                std::pair<int, int> edge(vertIndices[0], vertIndices[1]);
                edges.push_back(edge);

                // add the edge to the fracture mapper if it is a fracture
                if (dataVal < 0) {
                    fractureEdges.push_back(edge);
                }
            }
            else if (curParseMode == Element) {
                // skip the data attached to an element
                std::istringstream iss(curLine);
                int dataVal;
                std::string tmp;
                iss >> dataVal;
                iss >> tmp;
                assert(tmp == ":");

                // read the edge indices of an element
                std::vector<unsigned int> edgeIndices;
                while (iss) {
                    unsigned int tmp;
                    iss >> tmp;
                    if (!iss)
                        break;
                    edgeIndices.push_back(tmp);
                    assert(tmp < edges.size());
                }
                assert(edgeIndices.size() == 3); // so far, we only support triangles

                // extract the vertex indices of the element
                std::vector<unsigned int> vertIndices;
                for (int i = 0; i < 3; ++i) {
                    bool haveFirstVertex = false;
                    for (int j = 0; j < int(vertIndices.size()); ++j) {
                        assert(edgeIndices[i] < edges.size());
                        if (int(vertIndices[j]) == edges[edgeIndices[i]].first) {
                            haveFirstVertex = true;
                            break;
                        }
                    }
                    if (!haveFirstVertex)
                        vertIndices.push_back(edges[edgeIndices[i]].first);

                    bool haveSecondVertex = false;
                    for (unsigned j = 0; j < vertIndices.size(); ++j) {
                        assert(edgeIndices[i] < edges.size());
                        if (int(vertIndices[j]) == edges[edgeIndices[i]].second) {
                            haveSecondVertex = true;
                            break;
                        }
                    }
                    if (!haveSecondVertex)
                        vertIndices.push_back(edges[edgeIndices[i]].second);
                }

                // check whether the element's vertices are given in
                // mathematically positive direction. if not, swap the
                // first two.
                Dune::FieldMatrix<Scalar, 2, 2> mat;
                mat[0] = vertexPos[vertIndices[1]];
                mat[0] -= vertexPos[vertIndices[0]];
                mat[1] = vertexPos[vertIndices[2]];
                mat[1] -= vertexPos[vertIndices[0]];
                assert(mat.determinant() != 0);
                if (mat.determinant() < 0)
                    std::swap(vertIndices[2], vertIndices[1]);

                // insert the element into the dune grid
                gridFactory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2),
                                          vertIndices);
            }
            else if (curParseMode == Finished) {
                assert(curLine.size()==0);
            }
        }

        gridPtr_ = GridPointer(gridFactory.createGrid());

        /////
        // add the fracture edges
        /////

        // first create a map of the dune to ART vertex indices
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView, Dune::MCMGVertexLayout > VertexMapper;
        VertexMapper vertexMapper(gridPtr_->leafView());
        std::vector<int> artToDuneVertexIndex(vertexPos.size());
        const auto gridView = gridPtr_->leafView();
        auto vIt = gridView.template begin</*codim=*/2>();
        const auto &vEndIt = gridView.template end</*codim=*/2>();
        for (; vIt != vEndIt; ++vIt) {
            int duneIdx = vertexMapper.map(*vIt);
            int artIdx = gridFactory.insertionIndex(*vIt);
            artToDuneVertexIndex[artIdx] = duneIdx;
        }

        // add all fracture edges (using DUNE indices)
        auto feIt = fractureEdges.begin();
        const auto &feEndIt = fractureEdges.end();
        for (; feIt != feEndIt; ++ feIt) {
            fractureMapper_.addFractureEdge(artToDuneVertexIndex[feIt->first],
                                            artToDuneVertexIndex[feIt->second]);
        }
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    { return *gridPtr_; };

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     */
    static void loadBalance()
    { gridPtr_->loadBalance(); };

    /*!
     * \brief Destroys the grid
     *
     * This is required to guarantee that the grid is deleted before
     * MPI_Comm_free is called.
     */
    static void deleteGrid()
    { gridPtr_ = GridPointer((Grid*) 0); }

    /*!
     * \brief Returns the fracture mapper
     *
     * The fracture mapper determines the topology of the fractures.
     */
    static FractureMapper &fractureMapper()
    { return fractureMapper_; }

private:
    static GridPointer gridPtr_;
    static FractureMapper fractureMapper_;
};

template <class TypeTag>
typename Ewoms::ArtGridCreator<TypeTag>::GridPointer ArtGridCreator<TypeTag>::gridPtr_;
template <class TypeTag>
typename Ewoms::ArtGridCreator<TypeTag>::FractureMapper ArtGridCreator<TypeTag>::fractureMapper_;
} // end namespace

#endif // EWOMS_ART_GRID_CREATOR_HH
