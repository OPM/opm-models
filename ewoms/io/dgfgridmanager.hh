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

#include <ewoms/parallel/mpihelper.hh>

#include <dune/grid/io/file/dgfparser.hh>

#include <ewoms/io/basegridmanager.hh>
#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <type_traits>
#include <string>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(GridGlobalRefinements);
}} // namespace Opm, Properties

namespace Ewoms {
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

    typedef Dune::GridPtr<Grid> GridPointer;

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

        gridPtr_ = GridPointer(dgfFileName.c_str(), Dune::MPIHelper::getCommunicator());

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
     * \brief Returns a reference to the DGF grid pointer.
     */
    GridPointer& dgfGridPointer()
    { return gridPtr_; }

    /*!
     * \brief Returns a constant reference to the DGF grid pointer.
     */
    const GridPointer& dgfGridPointer() const
    { return gridPtr_; }

    /*!
     * \brief Distributes the grid on all processes of a parallel
     *        computation.
     *
     * This grid manager plays nice and also distributes the data of
     * the DGF...
     */
    void loadBalance()
    { gridPtr_.loadBalance(); }

private:
    GridPointer gridPtr_;
};

} // namespace Ewoms

#endif
