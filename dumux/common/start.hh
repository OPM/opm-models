// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \brief Provides a few default main functions for convenience.
 */
#ifndef DUMUX_START_HH
#define DUMUX_START_HH

#include <iostream>

#include "propertysystem.hh"
#include "parameters.hh"
#include "valgrind.hh"

#include <dune/common/mpihelper.hh>
#include <dune/common/version.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/dgfparser.hh>

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridCreator);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(TimeManager);
}

//! \cond 0
/*!
 * \ingroup Start
 *
 * \brief Prints the default usage message for the startWithParameters() function
 */
void printDefaultUsage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::cout << errorMsg << "\n"
                  << "\n";
    }
    std::cout
        << 
        "Usage: " << progName << " [options]\n"
        "Mandatory options include:\n"
        "\t--t-end=ENDTIME                  The time of the end of the simlation [s]\n"
        "\t--dt-initial=STEPSIZE            The initial time step size [s]\n"
        "\n"
        "Alternative syntax:\n"
        "\t-tEnd ENDTIME                    The time of the end of the simlation [s]\n"
        "\t-dtInitial STEPSIZE              The initial time step size [s]\n"
        "\n"
        "If --parameter-file is specified parameters can also be defined there. In this case,\n"
        "camel case is used for the parameters (e.g.: --t-end becomes tEnd). Parameters\n"
        "specified on the command line have priority over those in the parameter file.\n"
        "Important optional options include:\n"
        "\t--help,-h                        Print this usage message and exit\n"
        "\t--print-parameters[=true|false]  Print the run-time modifiable parameters _after_ \n"
        "\t                                 the simulation [default: true]\n"
        "\t--print-properties[=true|false]  Print the compile-time parameters _before_ \n"
        "\t                                 the simulation [default: false]\n"
        "\t--parameter-file=FILENAME        File with parameter definitions\n"
        "\t--restart=RESTARTTIME            Restart simulation from a restart file\n"
        "\n"
        "If the --parameter-file option is not specified, the program tries to load the input\n"
        "parameter file '"<<progName <<".input'\n"
        "\n";
}

/*!
 * \ingroup Start
 * \brief Read the command line arguments and write them into the parameter tree.
 *        Do some syntax checks.
 *
 * \param   argc      The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
 * \param   argv      The 'argv' argument of the main function: array of pointers to the argument strings
 * \param   paramTree The parameterTree. It can be filled from an input file or the command line.
 * \return            Empty string if everything worked out. Otherwise the thing that could not be read.
 */
std::string readOptions_(int argc, char **argv, Dune::ParameterTree &paramTree)
{
    // All command line options need to start with '-'
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            std::ostringstream oss;
            oss << "Command line argument " << i << " (='" << argv[i] << "') is invalid.";
            return oss.str();
        }

        std::string paramName, paramValue;

        // read a --my-opt=abc option. This gets transformed
        // into the parameter "myOpt" with the value being
        // "abc"
        if (argv[i][1] == '-') {
            std::string s(argv[i] + 2);
            // There is nothing after the '-'
            if (s.size() == 0 || s[0] == '=')
            {
                std::ostringstream oss;
                oss << "Parameter name of argument " << i << " ('" << argv[i] << "')"
                    << " is empty.";
                return oss.str();
            }

            // parse argument
            int j = 0;
            while (true) {
                if (j >= s.size()) {
                    // encountered the end of the string, i.e. we
                    // have a parameter where the argument is empty
                    paramName = s;
                    paramValue = "";
                    break;
                }
                else if (s[j] == '=') {
                    // we encountered a '=' character. everything
                    // before is the name of the parameter,
                    // everything after is the value.
                    paramName = s.substr(0, j);
                    paramValue = s.substr(j+1);
                    break;
                }
                else if (s[j] == '-') {
                    // remove all "-" characters and capitalize the
                    // character after them
                    s.erase(j, 1);
                    if (s.size() == j)
                    {
                        std::ostringstream oss;
                        oss << "Parameter name of argument " << i << " ('" << argv[i] << "')"
                            << " is invalid (ends with a '-' character).";
                        return oss.str();
                    }
                    else if (s[j] == '-')
                    {
                        std::ostringstream oss;
                        oss << "Malformed parameter name name in argument " << i << " ('" << argv[i] << "'): "
                            << "'--' in parameter name.";
                        return oss.str();
                    }
                    s[j] = std::toupper(s[j]);
                }
                else if (s[j] == '.') {
                    if (j + 1 >= s.size() || 
                        !std::isalpha(s[j+1]))
                    {
                        std::ostringstream oss;
                        oss << "Parameter name of argument " << i << " ('" << argv[i] << "')"
                            << " is invalid (character after '.' is not a letter).";
                        return oss.str();
                    }
                        
                    s[j + 1] = std::toupper(s[j + 1]);
                }
                
                ++j;
            }
        }
        else {
            // read a -myOpt abc option
            paramName = argv[i] + 1;

            if (argc == i + 1 || argv[i+1][0] == '-') {
                std::ostringstream oss;
                oss << "No argument given for parameter '" << argv[i] << "'!";
                return oss.str();
            }

            paramValue = argv[i+1];
            ++i; // In the case of '-myOpt abc' each pair counts as two arguments
        }

        // capitalize first letter of parameter name
        paramName[0] = std::toupper(paramName[0]);

        // Put the key=value pair into the parameter tree
        paramTree[paramName] = paramValue;
    }
    return "";
}

/*!
 * \ingroup Start
 *
 * \brief Provides a main function which reads in parameters from the
 *        command line and a parameter file.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param   argc    The 'argc' argument of the main function: count of arguments (1 if there are no arguments)
 * \param   argv    The 'argv' argument of the main function: array of pointers to the argument strings
 * \param   usage   Callback function for printing the usage message
 */
template <class TypeTag>
int start(int argc,
          char **argv,
          void (*usage)(const char *, const std::string &) = 0)
{
    // initialize MPI, finalize is done automatically on exit
    const Dune::MPIHelper &mpiHelper = Dune::MPIHelper::instance(argc, argv);
    
    // make sure that we always have a meaningful usage function
    if (!usage)
        usage = printDefaultUsage;
    
    int myRank = mpiHelper.rank();

    try {
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator; 
        typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////

        // check whether the user wanted to see the help message
        for (int i = 1; i < argc; ++i) {
            if (std::string("--help") == argv[i] || std::string("-h") == argv[i])
            {
                if (myRank == 0)
                    usage(argv[0], "");
                return 0;
            }
        }

        // fill the parameter tree with the options from the command line
        typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
        std::string s = readOptions_(argc, argv, ParameterTree::tree());
        if (!s.empty()) {
            if (myRank == 0)
                usage(argv[0], s);
            return 1;
        }

        if (!ParameterTree::tree().hasKey("ParameterFile")) {
            std::string paramFileName;
            paramFileName += argv[0];
            paramFileName += ".input";
            
            std::ifstream ifs(paramFileName.c_str());
            if (ifs) {
                ifs.close();
                if (myRank == 0)
                    std::cout << "Using fallback parameter file '" << paramFileName << "'\n";
                ParameterTree::tree()["ParameterFile"] = paramFileName;
            }
        }

        if (ParameterTree::tree().hasKey("ParameterFile")) {
            std::string inputFileName =
                GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile);

            // check whether the parameter file is readable.
            std::ifstream tmp;
            tmp.open(inputFileName.c_str());
            if (!tmp.is_open()){
                std::ostringstream oss;
                if (myRank == 0) {
                    oss << "Parameter file \"" << inputFileName << "\" is does not exist or is not readable.";
                    usage(argv[0], oss.str());
                }
                return 1;
            }

            // read the parameter file.
            Dune::ParameterTreeParser::readINITree(inputFileName,
                                                   ParameterTree::tree(),
                                                   /*overwrite=*/false);
        }

        bool printProps = true;
        if (ParameterTree::tree().hasKey("PrintProperties"))
            printProps = GET_RUNTIME_PARAM(TypeTag, bool, PrintProperties);

        if (printProps && myRank == 0) {
            Dumux::Properties::print<TypeTag>();
        }

        // deal with the restart stuff
        bool restart = false;
        Scalar restartTime = 0;
        if (ParameterTree::tree().hasKey("RestartTime")) {
            restart = true;
            restartTime = GET_RUNTIME_PARAM(TypeTag, Scalar, RestartTime);
        }

        // read the PrintParams parameter
        bool printParams = true;
        if (ParameterTree::tree().hasKey("PrintParameters"))
            printParams = GET_RUNTIME_PARAM(TypeTag, bool, PrintParameters);

        // try to create a grid (from the given grid file)
        try { GridCreator::makeGrid(); }
        catch (...) { if (myRank == 0) usage(argv[1], "Creation of the grid failed!"); throw; }
        GridCreator::loadBalance();

        // read the initial time step and the end time
        double tEnd;
        double dt;

        try { tEnd = GET_RUNTIME_PARAM(TypeTag, Scalar, TEnd); }
        catch (...) { if (myRank == 0) usage(argv[1], "Mandatory parameter '--t-end' not specified!"); throw; }

        try { dt = GET_RUNTIME_PARAM(TypeTag, Scalar, DtInitial); }
        catch (...) { if (myRank == 0) usage(argv[1], "Mandatory parameter '--dt-initial' not specified!"); throw; }

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager);
        timeManager.init(problem, restartTime, dt, tEnd, restart);
        timeManager.run();

        if (printParams && myRank == 0) {
            Dumux::Parameters::print<TypeTag>();
        }
        return 0;
    }
    catch (Dumux::ParameterException &e) {
        if (myRank == 0)
            std::cerr << e << ". Abort!\n";
        return 1;
    }
    catch (Dune::Exception &e) {
        if (myRank == 0)
            std::cerr << "Dune reported error: " << e << std::endl;
        return 2;
    }
    catch (...) {
        if (myRank == 0)
            std::cerr << "Unknown exception thrown!\n";
        return 3;
    }
}

} // namespace Dumux

#endif
