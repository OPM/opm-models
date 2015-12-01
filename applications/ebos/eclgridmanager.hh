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
 * \copydoc Ewoms::EclGridManager
 */
#ifndef EWOMS_ECL_GRID_MANAGER_HH
#define EWOMS_ECL_GRID_MANAGER_HH

#include <ewoms/io/basegridmanager.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/common/cartesianindexmapper.hh>
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/polyhedralgrid.hh>

#include <ewoms/io/polyhedralgridconverter.hh>

// set the EBOS_USE_ALUGRID macro. using the preprocessor for this is slightly hacky, but
// the macro is only used by this file...
#if EBOS_USE_ALUGRID
//#define DISABLE_ALUGRID_SFC_ORDERING 1
#if !HAVE_DUNE_ALUGRID || !DUNE_VERSION_NEWER(DUNE_ALUGRID, 2,4)
#warning "ALUGrid was indicated to be used for the ECL black oil simulator, but this "
#warning "requires the presence of dune-alugrid >= 2.4. Falling back to Dune::CpGrid"
#undef EBOS_USE_ALUGRID
#define EBOS_USE_ALUGRID 0
#endif
#else
#define EBOS_USE_ALUGRID 0
#endif

#if EBOS_USE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/common/fromtogridfactory.hh>
#endif

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseMode.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/checkDeck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>

namespace Ewoms {
template <class TypeTag>
class EclProblem;

template <class TypeTag>
class EclGridManager;

namespace Properties {
NEW_TYPE_TAG(EclGridManager);

// declare the properties required by the for the ecl grid manager
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(EclDeckFileName);

SET_STRING_PROP(EclGridManager, EclDeckFileName, "data/ecl.DATA");

// set the Grid and GridManager properties
#if EBOS_USE_ALUGRID
SET_TYPE_PROP(EclGridManager, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>);
#else
SET_TYPE_PROP(EclGridManager, Grid, Dune::CpGrid);
#endif

SET_TYPE_PROP(EclGridManager, GridManager, Ewoms::EclGridManager<TypeTag>);
} // namespace Properties

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 */
template <class TypeTag>
class EclGridManager : public BaseGridManager<TypeTag>
{
    typedef BaseGridManager<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

protected:
    typedef std::unique_ptr<Grid> GridPointer;

#if EBOS_USE_ALUGRID
    typedef Dune::PolyhedralGrid< Grid::dimension, Grid::dimensionworld > EquilGrid;
    struct EquilGridDeleter
    {
        void operator()(EquilGrid* polyhedralGrid)
        {
            if( polyhedralGrid )
            {
                const UnstructuredGrid& ug = static_cast< const UnstructuredGrid& > (*polyhedralGrid);

                UnstructuredGrid* ugPtr = (UnstructuredGrid *) &ug;
                delete polyhedralGrid;
                destroy_grid( ugPtr );
            }
        }
    };

    typedef std::unique_ptr< EquilGrid, EquilGridDeleter > EquilGridPointer;
#else
    typedef Dune::CpGrid EquilGrid;
    typedef std::unique_ptr< EquilGrid >  EquilGridPointer;
#endif

    typedef Dune :: CartesianIndexMapper< Grid > CartesianIndexMapper ;
    typedef std::unique_ptr< CartesianIndexMapper > CartesianIndexMapperPointer;


    static const int dimension = Grid :: dimension;
public:
    /*!
     * \brief Register all run-time parameters for the grid manager.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, EclDeckFileName,
                             "The name of the file which contains the ECL deck to be simulated");
    }

    /*!
     * \brief Create the grid for problem data files which use the ECL file format.
     *
     * This is the file format used by the commercial ECLiPSE simulator. Usually it uses
     * a cornerpoint description of the grid.
     */
    EclGridManager(Simulator &simulator)
        : ParentType(simulator)
    {
        std::string fileName = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);

        // compute the base name of the input file name
        const char directorySeparator = '/';
        long int i;
        for (i = fileName.size(); i >= 0; -- i)
            if (fileName[i] == directorySeparator)
                break;
        std::string baseName = fileName.substr(i + 1, fileName.size());

        // remove the extension from the input file
        for (i = baseName.size(); i >= 0; -- i)
            if (baseName[i] == '.')
                break;
        std::string rawCaseName;
        if (i < 0)
            rawCaseName = baseName;
        else
            rawCaseName = baseName.substr(0, i);

        // transform the result to ALL_UPPERCASE
        caseName_ = "";
        for (size_t i = 0; i < rawCaseName.size(); ++i)
            caseName_ += std::toupper(rawCaseName[i]);

        Opm::ParserPtr parser(new Opm::Parser());
        typedef std::pair<std::string, Opm::InputError::Action> ParseModePair;
        typedef std::vector<ParseModePair> ParseModePairs;
        ParseModePairs tmp;
        tmp.push_back(ParseModePair(Opm::ParseMode::PARSE_RANDOM_SLASH , Opm::InputError::IGNORE));
        Opm::ParseMode parseMode(tmp);
        std::cout << "Reading the deck file ('" << fileName << "')" << std::endl;
        deck_ = parser->parseFile(fileName , parseMode);
        eclState_.reset(new Opm::EclipseState(deck_, parseMode));

        Dune::CpGrid* cpgrid = new Dune::CpGrid();
        std::vector<double> porv = eclState_->getDoubleGridProperty("PORV")->getData();
        cpgrid->processEclipseFormat(eclState_->getEclipseGrid(),
                                      /*isPeriodic=*/false,
                                      /*flipNormals=*/false,
                                      /*clipZ=*/false,
                                      porv);

#if EBOS_USE_ALUGRID
        // copy cartesian ids
        std::vector<int> cartesianCellId( cpgrid->globalCell() );
        std::array<int,dimension> cartesianDimension;

        for (unsigned i = 0; i < dimension; ++i)
            cartesianDimension[i] = cpgrid->logicalCartesianSize()[i];

        Dune::FromToGridFactory< Grid > factory;
        grid_ = GridPointer(factory.convert(*cpgrid, cartesianCellId));
        // store cpgrid for equil initialization
        cartesianIndexMapper_.reset( new CartesianIndexMapper( *grid_, cartesianDimension, cartesianCellId ) );
#else
        grid_ = GridPointer(cpgrid);
        cartesianIndexMapper_.reset( new CartesianIndexMapper( *grid_ ) );
#endif

        this->finalizeInit_();
    }

    void loadBalance()
    {
#if EBOS_USE_ALUGRID
        releaseEquilGrid();

        auto gridView = grid().leafGridView();
        auto dataHandle = cartesianIndexMapper_->dataHandle( gridView );
        grid().loadBalance( *dataHandle );

        // communicate non-interior cells values
        grid().communicate( *dataHandle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication );

#else
        // distribute the grid and switch to the distributed view
        grid().loadBalance();
        grid().switchToDistributedView();
#endif
    }

    /*!
     * \brief Return a reference to the grid.
     */
    Grid& grid()
    { return *grid_; }

    /*!
     * \brief Return a reference to the grid.
     */
    const Grid& grid() const
    { return *grid_; }

    /*!
     * \brief Return a pointer to the parsed ECL deck
     */
    Opm::DeckConstPtr deck() const
    { return deck_; }

    /*!
     * \brief Return a pointer to the internalized ECL deck
     */
    Opm::EclipseStateConstPtr eclState() const
    { return eclState_; }

    /*!
     * \brief Return a pointer to the internalized schedule of the ECL deck
     */
    Opm::ScheduleConstPtr schedule() const
    { return eclState_->getSchedule(); }

    /*!
     * \brief Return a pointer to the EclipseGrid object
     *
     * The EclipseGrid class is provided by the opm-parser module and is used to
     * internalize the cornerpoint grid representation and, amongst others, can be used
     * to write EGRID files (which tends to be difficult with a plain Dune::CpGrid)
     */
    Opm::EclipseGridConstPtr eclGrid() const
    { return eclState_->getEclipseGrid(); }

    /*!
     * \brief Returns the name of the case.
     *
     * i.e., the all-uppercase version of the file name from which the
     * deck is loaded with the ".DATA" suffix removed.
     */
    const std::string& caseName() const
    { return caseName_; }

    /*!
     * \brief Returns the number of Cartesian cells in each direction
     */
    const std::array<int, dimension>& cartesianDimensions() const
    {
        return cartesianIndexMapper_->cartesianDimensions();
    }

    /*!
     * \brief Returns the overall number of Cartesian cells
     */
    int cartesianSize() const
    {
        return cartesianIndexMapper_->cartesianSize();
    }

    /*!
     * \brief Returns the Cartesian cell id for identifaction with Ecl data
     */
    unsigned cartesianIndex(unsigned compressedCellIdx) const
    {
        return cartesianIndexMapper_->cartesianIndex( compressedCellIdx );
    }

    /** \brief return index of the cells in the logical Cartesian grid */
    unsigned cartesianIndex( const std::array<int,dimension>& coords ) const
    {
#if EBOS_USE_ALUGRID
        return cartesianIndexMapper_->cartesianIndex( coords );
#else
        unsigned cartIndex = coords[ 0 ];
        int factor = cartesianDimensions()[ 0 ];
        for( unsigned i=1; i<dimension; ++i )
        {
            cartIndex += coords[ i ] * factor;
            factor *= cartesianDimensions()[ i ];
        }
        return cartIndex;
#endif
    }


    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void cartesianCoordinate(unsigned cellIdx, std::array<int,3>& ijk) const
    {
        return cartesianIndexMapper_->cartesianCoordinate( cellIdx, ijk );
    }

    const EquilGrid& equilGrid() const
    {
#if EBOS_USE_ALUGRID
        if( ! equilgrid_ )
            createEquilGrid();

        assert( equilgrid_ );
        return *equilgrid_;
#else
        return grid();
#endif
    }

    void releaseEquilGrid()
    {
        equilgrid_.reset();
    }

    void createEquilGrid() const
    {
#if EBOS_USE_ALUGRID
        UnstructuredGrid* ug = dune2UnstructuredGrid( grid_->leafGridView(),
                                                      *cartesianIndexMapper_,
                                                      true, false );
        if( ug )
        {
            equilgrid_ = EquilGridPointer( new EquilGrid( *ug ) );
        }
#endif
    }
private:
    std::string caseName_;
    GridPointer   grid_;
    mutable EquilGridPointer equilgrid_;
    CartesianIndexMapperPointer  cartesianIndexMapper_;
    Opm::DeckConstPtr deck_;
    Opm::EclipseStateConstPtr eclState_;

};

} // namespace Ewoms

#endif
