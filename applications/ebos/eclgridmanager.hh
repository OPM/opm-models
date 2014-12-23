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

#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/io/basegridmanager.hh>
#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/CpGrid.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Log/Logger.hpp>
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
} // namespace Ewoms

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(EclGridManager);

// declare the properties required by the for the ecl grid manager
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(EclDeckFileName);

SET_STRING_PROP(EclGridManager, EclDeckFileName, "data/ecl.DATA");

// set the Grid and GridManager properties
SET_TYPE_PROP(EclGridManager, Grid, Dune::CpGrid);
SET_TYPE_PROP(EclGridManager, GridManager, Ewoms::EclGridManager<TypeTag>);
}} // namespace Opm, Properties

namespace Ewoms {
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
        Opm::LoggerPtr opmLog(new Opm::Logger());

        try {
            deck_ = parser->parseFile(fileName, opmLog);
            Opm::checkDeck(deck_, opmLog);
            eclState_.reset(new Opm::EclipseState(deck_, opmLog));
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "Non-recoverable error encountered while parsing the deck file:"
                      << e.what() << "\n";

            if (opmLog->size() > 0) {
                std::cerr << "Issues found while parsing the deck file:\n";
                opmLog->printAll(std::cerr);
            }

            throw e;
        }

        if (opmLog->size() > 0) {
            std::cout << "Issues found while parsing the deck file:\n";
            opmLog->printAll(std::cout);
        }

        grid_ = GridPointer(new Grid());
        grid_->processEclipseFormat(eclState_->getEclipseGrid(),
                                    /*isPeriodic=*/false,
                                    /*flipNormals=*/false,
                                    /*clipZ=*/false);

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

private:
    std::string caseName_;
    GridPointer grid_;
    Opm::DeckConstPtr deck_;
    Opm::EclipseStateConstPtr eclState_;
};

} // namespace Ewoms

#endif
