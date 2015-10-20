// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
 * \copydoc Ewoms::DgfGridManager
 */
#ifndef EWOMS_DGF_GRID_MANAGER_HH
#define EWOMS_DGF_GRID_MANAGER_HH

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <ewoms/models/discretefracture/fracturemapper.hh>

#include <ewoms/io/basegridmanager.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>


#include <type_traits>
#include <string>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(GridManager);
NEW_PROP_TAG(GridGlobalRefinements);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Simulator);
} // namespace Properties

/*!
 * \brief Provides a grid manager which reads Dune Grid Format (DGF) files
 */
template <class TypeTag>
class DgfGridManager : public BaseGridManager<TypeTag>
{
    typedef BaseGridManager<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef Ewoms::FractureMapper<TypeTag> FractureMapper;

    typedef std::unique_ptr< Grid > GridPointer;

public:
    /*!
     * \brief Register all run-time parameters for the grid manager.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, GridFile,
                             "The file name of the DGF file to load");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, GridGlobalRefinements,
                             "The number of global refinements of the grid "
                             "executed after it was loaded");
    }

    /*!
     * \brief Load the grid from the file.
     */
    DgfGridManager(Simulator &simulator)
        : ParentType(simulator)
    {
        const std::string dgfFileName = EWOMS_GET_PARAM(TypeTag, std::string, GridFile);
        unsigned numRefinments = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);

        {
            // create DGF GridPtr from a dgf file
            Dune::GridPtr< Grid > dgfPointer( dgfFileName );

            // this is only implemented for 2d currently
            addFractures( dgfPointer );

            // store pointer to dune grid
            gridPtr_.reset( dgfPointer.release() );
        }

        if (numRefinments > 0)
            gridPtr_->globalRefine(numRefinments);

        this->finalizeInit_();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return *gridPtr_; }

    /*!
     * \brief Returns a reference to the grid.
     */
    const Grid& grid() const
    { return *gridPtr_; }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     *
     * This grid manager plays nice and also distributes the data of
     * the DGF...
     */
    void loadBalance()
    { gridPtr_->loadBalance(); }

    /*!
     * \brief Returns the fracture mapper
     *
     * The fracture mapper determines the topology of the fractures.
     */
    FractureMapper &fractureMapper()
    { return fractureMapper_; }

    /*!
     * \brief Returns the fracture mapper
     *
     * The fracture mapper determines the topology of the fractures.
     */
    const FractureMapper &fractureMapper() const
    { return fractureMapper_; }

protected:
    void addFractures( Dune::GridPtr< Grid >& dgfPointer )
    {
        // check if fractures are available (only 2d currently)
        if( dgfPointer.nofParameters( int(Grid::dimension) ) > 0 )
        {
            typedef typename  Grid::LevelGridView GridView;
            GridView gridView = dgfPointer->levelGridView( 0 );

            // first create a map of the dune to ART vertex indices
            typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
                                                              Dune::MCMGVertexLayout>  ElementMapper;

            const int edgeCodim = Grid::dimension - 1;

            ElementMapper elementMapper( gridView );
            const auto endIt = gridView.template end< 0 > ();
            for( auto eIt = gridView.template begin< 0 >(); eIt != endIt; ++eIt )
            {
                const auto& element = *eIt;
                const Dune::ReferenceElement< Scalar, Grid::dimension > & refElem =
                      Dune::ReferenceElements< Scalar, Grid::dimension >::general( element.type() );

                const int edges = refElem.size( edgeCodim );
                for( int edge = 0; edge < edges; ++edge )
                {
                    const int vertices = refElem.size( edge, edgeCodim, Grid::dimension );
                    std::vector< int > vertexIndices;
                    vertexIndices.reserve( Grid::dimension );
                    for( int vx = 0; vx<vertices; ++vx )
                    {
                        // get local vertex number from edge
                        const int localVx = refElem.subEntity( edge, edgeCodim, vx, Grid::dimension );

#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
                        // get vertex
                        const auto vertex = element.template subEntity< Grid::dimension >( localVx );
#else
                        // get vertex
                        const auto vertexPtr = element.template subEntity< Grid::dimension >( localVx );
                        const auto& vertex = *vertexPtr ;
#endif

                        // if vertex has parameter 1 insert as a fracture vertex
                        if( dgfPointer.parameters( vertex )[ 0 ] > 0 )
                        {
                            vertexIndices.push_back( elementMapper.subIndex( element, localVx, Grid::dimension ) );
                        }
                    }
                    // if 2 vertices have been found with flag 1 insert a fracture edge
                    if( int(vertexIndices.size()) == Grid::dimension )
                    {
                        fractureMapper_.addFractureEdge(vertexIndices[ 0 ], vertexIndices[ 1 ] );
                    }
                }
            }
        }
    }

private:
    GridPointer    gridPtr_;
    FractureMapper fractureMapper_;
};

} // namespace Ewoms

#endif
