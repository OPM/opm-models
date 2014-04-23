/*
  Copyright (C) 2014 by Andreas Lauser

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
 *
 * \copydoc Ewoms::BlackOilProblem
 */
#ifndef EWOMS_BLACKOIL_PROBLEM_HH
#define EWOMS_BLACKOIL_PROBLEM_HH

#include "blackoilproperties.hh"

#include <ewoms/models/common/multiphasebaseproblem.hh>
#include <ewoms/io/eclipsewriter.hh>

namespace Ewoms {

/*!
 * \ingroup Discretization
 *
 * \brief Base class for all problems which use the black-oil model.
 */
template<class TypeTag>
class BlackOilProblem : public MultiPhaseBaseProblem<TypeTag>
{
private:
    typedef MultiPhaseBaseProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     *
     * \param simulator The manager object of the simulation
     */
    BlackOilProblem(Simulator &simulator)
        : ParentType(simulator)
        , eclipseWriter_(simulator)
    {}

    /*!
     * \brief Registers all available parameters for the problem and
     *        the model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEclipseOutput,
                             "Write binary output which is compatible with the commercial "
                             "Eclipse simulator");
    }

    /*!
     * \brief Write the relevant secondary variables of the current
     *        solution into the output files.
     *
     * \param verbose If true, then a message will be printed to stdout if a file is written
     */
    void writeOutput(bool verbose = true)
    {
        if (!asImp_().shouldWriteOutput())
            return;

        // calculate the time _after_ the time was updated
        Scalar t = this->simulator().time() + this->simulator().timeStepSize();

        // prepare the Eclipse and the VTK writers
        if (enableEclipseOutput_())
            eclipseWriter_.beginWrite(t);

        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        ParentType::writeOutput(verbose);

        if (enableEclipseOutput_()) {
            this->model().appendOutputFields(eclipseWriter_);
            eclipseWriter_.endWrite();
        }
    }

private:
    static bool enableEclipseOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclipseOutput); }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    EclipseWriter<TypeTag> eclipseWriter_;
};

} // namespace Ewoms

#endif
