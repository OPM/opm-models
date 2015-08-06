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

#include <dune/grid/CpGrid.hpp>

// set the EBOS_USE_ALUGRID macro. using the preprocessor for this is slightly hacky, but
// the macro is only used by this file...
#if EBOS_USE_ALUGRID
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
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    typedef std::unique_ptr<Grid> GridPointer;

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
        : ParentType(simulator),
          cartesianCellId_(),
          cartesianSize_()
    {
        std::string fileName = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);

        // compute the base name of the input file name
        const char directorySeparator = '/';
        int i;
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

        std::cout << "Reading the deck file ('" << fileName << "')" << std::endl;
        Opm::ParseMode parseMode;
        deck_ = parser->parseFile(fileName, parseMode);
        eclState_.reset(new Opm::EclipseState(deck_, parseMode));

#if EBOS_USE_ALUGRID
        std::unique_ptr< Dune::CpGrid > cpgrid(new Dune::CpGrid());
#else
        Grid* cpgrid = new Grid();
#endif
        cpgrid->processEclipseFormat(eclState_->getEclipseGrid(),
                                     /*isPeriodic=*/false,
                                     /*flipNormals=*/false,
                                     /*clipZ=*/false);

        for (int i = 0; i < dimension; ++i)
            cartesianSize_[i] = cpgrid->logicalCartesianSize()[i];

#if EBOS_USE_ALUGRID
        Dune::FromToGridFactory< Grid > factory;
        std::vector< int > ordering;
        grid_ = GridPointer(factory.convert(*cpgrid, ordering));
        if (ordering.empty())
            // copy cartesian cell index from cp grid
            cartesianCellId_ = cpgrid->globalCell();
        else {
            const int size = ordering.size();
            cartesianCellId_.reserve(size);
            const std::vector<int>& globalCell = cpgrid->globalCell();
            for (int i = 0; i < size; ++i)
                cartesianCellId_.push_back(globalCell[ordering[i]]);
        }
#else
        grid_ = GridPointer(cpgrid);
#endif

        this->finalizeInit_();
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
     * \brief Returns the logical Cartesian size
     */
    const std::array<int, dimension>& logicalCartesianSize() const
    { return cartesianSize_; }

    /*!
     * \brief Returns the logical Cartesian size
     */
    int numLogicalCartesianCells() const
    {
#if EBOS_USE_ALUGRID
        int n = cartesianCellId_.size();
#else
        int n = grid_->globalCell().size();
#endif
        assert(cartesianSize_[0]*cartesianSize_[1]*cartesianSize_[2] == n);
        return n;
    }

    /*!
     * \brief Returns the Cartesian cell id for identifaction with Ecl data
     */
    int cartesianCellId(int compressedCellIdx) const
    {
#if EBOS_USE_ALUGRID
        return cartesianCellId_[compressedCellIdx];
#else
        return grid_->globalCell()[compressedCellIdx];
#endif
    }

    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void getIJK(int cellIdx, std::array<int,3>& ijk) const
    {
        assert(cellIdx < int(numLogicalCartesianCells()));
        int cartesianCellIdx = cartesianCellId(cellIdx);

        ijk[0] = cartesianCellIdx % cartesianSize_[0];
        cartesianCellIdx /= cartesianSize_[0];

        ijk[1] = cartesianCellIdx % cartesianSize_[1];

        ijk[2] = cartesianCellIdx / cartesianSize_[1];

#if !defined NDEBUG && !EBOS_USE_ALUGRID
        // make sure ijk computation is the same as in CpGrid
        std::array<int,3> checkIjk;
        grid_->getIJK(cellIdx, checkIjk);
        for (int i=0; i < 3; ++i)
            assert(checkIjk[i] == ijk[i]);
#endif
    }

private:
    std::string caseName_;
    GridPointer grid_;
    Opm::DeckConstPtr deck_;
    Opm::EclipseStateConstPtr eclState_;

    std::vector<int> cartesianCellId_;
    std::array<int,dimension> cartesianSize_;
};

} // namespace Ewoms

#endif
