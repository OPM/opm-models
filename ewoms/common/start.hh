/*
  Copyright (C) 2010-2013 by Andreas Lauser
  Copyright (C) 2011 by Philipp Nuske
  Copyright (C) 2012 by Bernd Flemisch

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
 * \brief Provides convenience routines to bring up the simulation at runtime.
 */
#ifndef EWOMS_START_HH
#define EWOMS_START_HH

#include <opm/core/utility/PropertySystem.hpp>
#include "parametersystem.hh"
#include <opm/material/Valgrind.hpp>

#include <ewoms/version.hh>
#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/simulator.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/common/parametertreeparser.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <locale>

#include <stdio.h>
#include <unistd.h>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Opm {
// forward declaration of property tags
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(EndTime);
NEW_PROP_TAG(RestartTime);
NEW_PROP_TAG(InitialTimeStepSize);
NEW_PROP_TAG(PrintProperties);
NEW_PROP_TAG(PrintParameters);
NEW_PROP_TAG(ParameterFile);
} // namespace Properties
} // namespace Opm
//! \cond SKIP_THIS

namespace Ewoms {
/*!
 * \brief Register all runtime parameters, parse the command line
 *        arguments and the parameter file.
 *
 * \param argc The number of command line arguments
 * \param argv Array with the command line argument strings
 */
template <class TypeTag>
int setupParameters_(int argc, char **argv)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP(TypeTag, ParameterMetaData) ParameterMetaData;

    // first, get the MPI rank of the current process
    int myRank = 0;
#if HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

    ////////////////////////////////////////////////////////////
    // Register all parameters
    ////////////////////////////////////////////////////////////
    EWOMS_REGISTER_PARAM(TypeTag, std::string, ParameterFile,
                         "An .ini file which contains a set of run-time "
                         "parameters");
    EWOMS_REGISTER_PARAM(TypeTag, int, PrintProperties,
                         "Print the values of the compile time properties at "
                         "the start of the simulation");
    EWOMS_REGISTER_PARAM(TypeTag, int, PrintParameters,
                         "Print the values of the run-time parameters at the "
                         "start of the simulation");
    EWOMS_REGISTER_PARAM(TypeTag, Scalar, EndTime,
                         "The time at which the simulation is finished [s]");
    EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitialTimeStepSize,
                         "The size of the initial time step [s]");
    EWOMS_REGISTER_PARAM(TypeTag, Scalar, RestartTime,
                         "The time time at which a simulation restart should "
                         "be attempted [s]");

    GridCreator::registerParameters();
    Simulator::registerParameters();

    ////////////////////////////////////////////////////////////
    // set the parameter values
    ////////////////////////////////////////////////////////////

    // check whether the user wanted to see the help message
    for (int i = 1; i < argc; ++i) {
        if (std::string("--help") == argv[i] || std::string("-h") == argv[i]) {
            if (myRank == 0)
                Parameters::printUsage<TypeTag>(argv[0], "", /*handleHelp=*/true, std::cout);
            return /*status=*/2;
        }
    }

    // fill the parameter tree with the options from the command line
    std::string s = Parameters::parseCommandLineOptions<TypeTag>(argc, argv);
    if (!s.empty()) {
        return /*status=*/1;
    }

    std::string paramFileName
        = EWOMS_GET_PARAM_(TypeTag, std::string, ParameterFile);
    if (paramFileName != "") {
        ////////////////////////////////////////////////////////////
        // add the parameters specified using an .ini file
        ////////////////////////////////////////////////////////////

        // check whether the parameter file is readable.
        std::ifstream tmp;
        tmp.open(paramFileName.c_str());
        if (!tmp.is_open()) {
            std::ostringstream oss;
            if (myRank == 0) {
                oss << "Parameter file \"" << paramFileName
                    << "\" does not exist or is not readable.";
                Parameters::printUsage<TypeTag>(argv[0], oss.str());
            }
            return /*status=*/1;
        }

        // read the parameter file.
        Dune::ParameterTreeParser::readINITree(paramFileName,
                                               ParameterMetaData::tree(),
                                               /*overwrite=*/false);
    }

    EWOMS_END_PARAM_REGISTRATION(TypeTag);

    return /*status=*/0;
}

//! \endcond

/*!
 * \ingroup Startup
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc The number of command line arguments
 * \param argv The array of the command line arguments
 */
template <class TypeTag>
int start(int argc, char **argv)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    // initialize MPI, finalize is done automatically on exit
    const Ewoms::MpiHelper mpiHelper(argc, argv);

    int myRank = mpiHelper.rank();

    try
    {
        int paramStatus = setupParameters_<TypeTag>(argc, argv);
        if (paramStatus == 1)
            return 1;
        if (paramStatus == 2)
            return 0;

        // read the initial time step and the end time
        double endTime;
        double initialTimeStepSize;

        endTime = EWOMS_GET_PARAM(TypeTag, Scalar, EndTime);
        if (endTime < -1e50) {
            if (myRank == 0)
                Parameters::printUsage<TypeTag>(argv[0],
                                                "Mandatory parameter '--end-time' not specified!");
            return 1;
        }

        initialTimeStepSize = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
        if (initialTimeStepSize < -1e50) {
            if (myRank == 0)
                Parameters::printUsage<TypeTag>(argv[0],
                                                "Mandatory parameter '--initial-time-step-size' "
                                                "not specified!");
            return 1;
        }

        // deal with the restart stuff
        bool doRestart = false;
        Scalar startTime = EWOMS_GET_PARAM(TypeTag, Scalar, RestartTime);
        if (startTime > -1e100)
            doRestart = true;
        else
            startTime = 0.0;

        if (myRank == 0)
            std::cout << "eWoms " << EWOMS_VERSION
#ifdef EWOMS_CODENAME
                      << " (\"" << EWOMS_CODENAME << "\")"
#endif
                      << " will now start the trip. "
                      << "Please sit back, relax and enjoy the ride.\n"
                      << std::flush;

        // print the parameters if requested
        int printParams = EWOMS_GET_PARAM(TypeTag, int, PrintParameters);
        if (myRank == 0) {
            std::string endParametersSeparator("# [end of parameters]\n");
            if (printParams) {
                bool printSeparator = false;
                if (printParams == 1 || !isatty(fileno(stdout))) {
                    Ewoms::Parameters::printValues<TypeTag>();
                    printSeparator = true;
                }
                else
                    // always print the list of specified but unused parameters
                    printSeparator =
                        printSeparator ||
                        Ewoms::Parameters::printUnused<TypeTag>();
                if (printSeparator)
                    std::cout << endParametersSeparator;
            }
            else
                // always print the list of specified but unused parameters
                if (Ewoms::Parameters::printUnused<TypeTag>())
                    std::cout << endParametersSeparator;
        }

        // print the properties if requested
        int printProps = EWOMS_GET_PARAM(TypeTag, int, PrintProperties);
        if (printProps && myRank == 0) {
            if (printProps == 1 || !isatty(fileno(stdout)))
                Opm::Properties::printValues<TypeTag>();
        }

        // instantiate and run the concrete problem. make sure to
        // deallocate the problem and before the time manager and the
        // grid
        Simulator simulator;
        simulator.init(startTime, initialTimeStepSize, endTime, doRestart);
        simulator.run();

        if (myRank == 0) {
            std::cout << "eWoms reached the destination. If it took a wrong "
                      << "corner, change your booking and try again.\n"
                      << std::flush;
        }
        return 0;
    }
    catch (std::exception &e)
    {
        if (myRank == 0)
            std::cout << e.what() << ". Abort!\n" << std::flush;
        return 1;
    }
    catch (Dune::Exception &e)
    {
        if (myRank == 0)
            std::cout << "Dune reported an error: " << e.what() << std::endl  << std::flush;
        return 2;
    }
    catch (...)
    {
        if (myRank == 0)
            std::cout << "Unknown exception thrown!\n" << std::flush;
        return 3;
    }
}

} // namespace Ewoms

#endif
