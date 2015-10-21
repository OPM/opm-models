// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2014 IRIS AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef EWOMS_POLYHEDRALGRIDCONVERTER_HH
#define EWOMS_POLYHEDRALGRIDCONVERTER_HH

//#if HAVE_DUNE_CORNERPOINT
#include <dune/grid/polyhedralgrid.hh>
// we need dune-cornerpoint for reading the Dune grid.
#include <dune/grid/CpGrid.hpp>
//#else
//#error "This header needs the dune-cornerpoint module"
//#endif

namespace Ewoms
{
  namespace detail{
    // key for generating intersection index
    struct FaceKey : public std::pair< int, int >
    {
      typedef std::pair< int, int > BaseType;
      FaceKey() : BaseType(-1,-1) {}
      FaceKey( const int inside, const int outside )
        : BaseType( inside < outside ? std::make_pair(inside,outside) : std::make_pair(outside,inside) )
      {}
    };
  }


  template <class GridView, class CartesianIndexMapper>
  inline UnstructuredGrid*
  dune2UnstructuredGrid( const GridView& gridView,
                         const CartesianIndexMapper& cartesianIndexMapper,
                         const bool faceTags,
                         const bool onlyInterior = true )
  {
      // if grid is empty return nullptr
      if( gridView.grid().size( 0 ) == 0 )
          return nullptr;

      typedef double ctype;
      typedef typename GridView :: template Codim< 0 > :: template Partition<
          Dune :: All_Partition > :: Iterator Iterator;
      typedef typename GridView :: template Codim< 0 > :: Entity         Element;
      typedef typename Element  :: Geometry                              ElementGeometry;
      typedef typename GridView :: IntersectionIterator           IntersectionIterator;
      typedef typename IntersectionIterator :: Intersection       Intersection;
      typedef typename Intersection :: Geometry                   IntersectionGeometry;
      typedef typename GridView :: IndexSet                       IndexSet;

      typedef typename ElementGeometry :: GlobalCoordinate GlobalCoordinate;

      const IndexSet& indexSet = gridView.indexSet();
      const int dimension = GridView::Grid::dimension;

      const int numCells = indexSet.size( 0 );
      const int numNodes = indexSet.size( dimension );

      const int maxNumVerticesPerFace = 4;
      const int maxNumFacesPerCell    = 6;

      int maxFaceIdx = indexSet.size( 1 );
      const bool validFaceIndexSet = maxFaceIdx > 0;

      typedef detail::FaceKey  FaceKey;

      std::map< FaceKey, int > faceIndexSet;

      if( ! validFaceIndexSet )
      {
          maxFaceIdx = 0;
          const Iterator end = gridView.template end<0, Dune::All_Partition> ();
          for( Iterator it = gridView.template begin<0, Dune::All_Partition> (); it != end; ++it )
          {
              const auto& element = *it;
              const int elIndex = indexSet.index( element );
              const IntersectionIterator endiit = gridView.iend( element );
              for( IntersectionIterator iit = gridView.ibegin( element ); iit != endiit; ++iit)
              {
                  const Intersection& intersection = *iit;
                  int nbIndex = -1;

                  // store face --> cell relation
                  if( intersection.neighbor() )
                  {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
                      const Element& neighbor = intersection.outside();
#else
                      auto  ep = intersection.outside();
                      const Element& neighbor = *ep;
#endif
                      if( ! ( onlyInterior && neighbor.partitionType() != Dune::InteriorEntity ) )
                        nbIndex = indexSet.index( neighbor );
                  }

                  FaceKey faceKey( elIndex, nbIndex );
                  if( faceIndexSet.find( faceKey ) == faceIndexSet.end() )
                      faceIndexSet[ faceKey ] = maxFaceIdx++;
              }
          }
      }
      const int numFaces = maxFaceIdx ;

      // create Unstructured grid struct
      UnstructuredGrid* ug = allocate_grid( dimension, numCells, numFaces,
              numFaces*maxNumVerticesPerFace,
              numCells * maxNumFacesPerCell,
              numNodes );

      std::fill( ug->face_cells, ug->face_cells+(numCells * maxNumFacesPerCell), -1 );

      for( int d=0; d<dimension; ++d )
        ug->cartdims[ d ] = cartesianIndexMapper.cartesianDimensions()[ d ];

      assert( ug->number_of_cells > 0 );
      // allocate data structure for storage of cartesian index
      if( ! ug->global_cell )
          ug->global_cell = (int *) std::malloc( ug->number_of_cells * sizeof(int) );

      int count = 0;
      int cellFace = 0;
      maxFaceIdx = 0;
      const Iterator end = gridView.template end<0, Dune::All_Partition> ();
      for( Iterator it = gridView.template begin<0, Dune::All_Partition> (); it != end; ++it, ++count )
      {
          const Element& element = *it;
          const ElementGeometry geometry = element.geometry();

          // currently only hexahedrons are supported
          // assert( element.type().isHexahedron() );

          const int elIndex = indexSet.index( element );
          assert( indexSet.index( element ) == elIndex );

          const bool isGhost = element.partitionType() != Dune :: InteriorEntity ;

          // make sure that the elements are ordered as before,
          // otherwise the globalCell mapping is invalid
          assert( count == elIndex );

          // store cartesian index
          ug->global_cell[ elIndex ] = cartesianIndexMapper.cartesianIndex( elIndex );
          //std::cout << "global index of cell " << elIndex << " = " <<
          //    ug->global_cell[ elIndex ] << std::endl;

          const GlobalCoordinate center = geometry.center();
          int idx = elIndex * dimension;
          for( int d=0; d<dimension; ++d, ++idx )
              ug->cell_centroids[ idx ] = center[ d ];
          ug->cell_volumes[ elIndex ] = geometry.volume();

          const int vertices = geometry.corners();
          for( int vx=0; vx<vertices; ++vx )
          {
              const GlobalCoordinate vertex = geometry.corner( vx );
              int idx = indexSet.subIndex( element, vx, dimension ) * dimension;
              for( int d=0; d<dimension; ++d, ++idx )
                  ug->node_coordinates[ idx ] = vertex[ d ];
          }

          ug->cell_facepos[ elIndex ] = cellFace;

          Dune::GeometryType geomType = element.type();
          if( geomType.isNone() )
              geomType = Dune::GeometryType( Dune::GeometryType::cube, dimension );

          const Dune::ReferenceElement< ctype, dimension > &refElem
                  = Dune::ReferenceElements< ctype, dimension >::general( geomType );

          int faceCount = 0;
          const IntersectionIterator endiit = gridView.iend( element );
          for( IntersectionIterator iit = gridView.ibegin( element ); iit != endiit; ++iit, ++faceCount )
          {
              const Intersection& intersection = *iit;
              IntersectionGeometry intersectionGeometry = intersection.geometry();
              const double faceVol = intersectionGeometry.volume();

              const int localFace = intersection.indexInInside();
              const int localFaceIdx = isGhost ? 0 : localFace;

              int faceIndex = validFaceIndexSet ? indexSet.subIndex( element, localFace, 1 ) : -1;
              if( ! validFaceIndexSet )
              {
                  int nbIndex = -1;
                  if( intersection.neighbor() )
                  {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
                      const Element& neighbor = intersection.outside();
#else
                      auto  ep = intersection.outside();
                      const Element& neighbor = *ep;
#endif
                      if( ! ( onlyInterior && neighbor.partitionType() != Dune::InteriorEntity ) )
                          nbIndex = indexSet.index( neighbor );
                  }
                  FaceKey faceKey( elIndex, nbIndex );
                  faceIndex = faceIndexSet[ faceKey ];
              }

              maxFaceIdx = std::max( faceIndex, maxFaceIdx );

              ug->face_areas[ faceIndex ] = faceVol;

              // get number of vertices (should be 4)
              const int vxSize = refElem.size( localFace, 1, dimension );
              int faceIdx = faceIndex * maxNumVerticesPerFace ;
              ug->face_nodepos[ faceIndex   ] = faceIdx;
              ug->face_nodepos[ faceIndex+1 ] = faceIdx + maxNumVerticesPerFace;
              for( int vx=0; vx<vxSize; ++vx, ++faceIdx )
              {
                  const int localVx = refElem.subEntity( localFace, 1, vx, dimension );
                  const int vxIndex = indexSet.subIndex( element, localVx, dimension );
                  ug->face_nodes[ faceIdx ] = vxIndex ;
              }

              assert( vxSize    <= maxNumVerticesPerFace );
              assert( localFace <  maxNumFacesPerCell );

              // store cell --> face relation
              ug->cell_faces  [ cellFace + localFaceIdx ] = faceIndex;
              if( faceTags )
              {
                  // fill logical cartesian orientation of the face (here indexInInside)
                  ug->cell_facetag[ cellFace + localFaceIdx ] = localFaceIdx;
              }

              GlobalCoordinate normal = intersection.centerUnitOuterNormal();
              normal *= faceVol;

              // store face --> cell relation
              if( intersection.neighbor() )
              {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,4)
                  const Element& neighbor = intersection.outside();
#else
                  auto  ep = intersection.outside();
                  const Element& neighbor = *ep;
#endif

                  int nbIndex = -1;
                  if( ! ( onlyInterior && neighbor.partitionType() != Dune::InteriorEntity ) )
                      nbIndex = indexSet.index( neighbor );

                  if( nbIndex == -1 || elIndex < nbIndex )
                  {
                      ug->face_cells[ 2*faceIndex     ] = elIndex;
                      ug->face_cells[ 2*faceIndex + 1 ] = nbIndex;
                  }
                  else
                  {
                      ug->face_cells[ 2*faceIndex     ] = nbIndex;
                      ug->face_cells[ 2*faceIndex + 1 ] = elIndex;
                      // flip normal
                      normal *= -1.0;
                  }
              }
              else // domain boundary
              {
                  ug->face_cells[ 2*faceIndex     ] = elIndex;
                  ug->face_cells[ 2*faceIndex + 1 ] = -1; // boundary
              }

              const GlobalCoordinate center = intersectionGeometry.center();
              // store normal
              int idx = faceIndex * dimension;
              for( int d=0; d<dimension; ++d, ++idx )
              {
                  ug->face_normals  [ idx ] = normal[ d ];
                  ug->face_centroids[ idx ] = center[ d ];
              }
          }
          if( faceCount > maxNumFacesPerCell )
              OPM_THROW(std::logic_error,"DuneGrid only supports conforming hexahedral currently");
          cellFace += faceCount;
      }

      // set last entry
      ug->cell_facepos[ numCells ] = cellFace;
      // set number of faces found
      ug->number_of_faces = maxFaceIdx+1;

      // std::cout << cellFace << " " << indexSet.size( 1 ) << " " << maxFaceIdx << std::endl;
      return ug;
  }
} // end namespace Opm

#endif
