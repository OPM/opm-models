// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 2012 by Philipp Nuske                                     *
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
 * \brief Source file of a test for the decoupled one-phase model.
 */
#include "config.h"

#include "test_1pproblem.hh"
#include "benchmarkresult.hh"
#include <dumux/common/start.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe List of Mandatory arguments for this program is:\n"
                                        "\t-refine                        The refinement level of the grid. [-] \n"
                                        "\t-tEnd                          The end of the simulation. [s] \n"
                                       "\t-Grid.NumRefine        The refinement level of the grid. [-] \n";
                    errorMessageOut += "\n\nAdditionaly the following arguments can be specified:\n"
                                       "\t-Problem.Delta         Anisotropy of permeability tensor. Value out"
                                       "\t                       of (0, 1], with 1 being isotrop. Default: 1e-3.\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    typedef TTAG(TestProblemOneP) ProblemTypeTag;

    return Dumux::start<ProblemTypeTag>(argc, argv, usage);
}
