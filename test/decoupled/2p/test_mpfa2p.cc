// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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

#define STRUCTUREDGRID 1

#include "config.h"

#if HAVE_ALUGRID

#include "test_mpfa2pproblem.hh"
#include <dumux/common/start.hh>

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    Dune::ParameterTree paramTree;
    std::string s(Dumux::readOptions_(argc, argv, paramTree));
    if (s.empty())
    {
        if (paramTree.hasKey("ModelType"))
        {
            std::string modelType(paramTree.get<std::string>("ModelType"));

            if (modelType == "FV")
            {
                typedef TTAG(FVTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::tree();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
                std::cout<<"######################################################\n";
                std::cout<<"Used standard finite volume model\n";
                return startReturn;
            }
            else if (modelType == "FVAdaptive")
            {
                typedef TTAG(FVAdaptiveTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::tree();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
                std::cout<<"######################################################\n";
                std::cout<<"Used adaptive finite volume model\n";
                return startReturn;
            }
            else if (modelType == "MPFAO")
            {
                typedef TTAG(MPFAOTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::tree();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
                std::cout<<"######################################################\n";
                std::cout<<"Used finite volume MPFA o-method model\n";
                return startReturn;
            }
            else if (modelType == "MPFAL")
            {
                typedef TTAG(MPFALTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::tree();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
                std::cout<<"######################################################\n";
                std::cout<<"Used finite volume MPFA l-method model\n";
                return startReturn;
            }
            else if (modelType == "MPFALAdaptive")
            {
                typedef TTAG(MPFALAdaptiveTwoPTestProblem) ProblemTypeTag;
                typedef GET_PROP(ProblemTypeTag, ParameterTree) ParamTree;
                Dune::ParameterTree &rt = ParamTree::tree();
                rt["ModelType"]=modelType;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
                std::cout<<"######################################################\n";
                std::cout<<"Used adaptive finite volume MPFA l-method model\n";
                return startReturn;
            }
            else
            {
                typedef TTAG(MPFAOTwoPTestProblem) ProblemTypeTag;
                int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
                std::cout<<"######################################################\n";
                std::cout<<"Unknwon model type "<<modelType<<" specified\n";
                std::cout<<"Default to finite volume MPFA o-method model\n";
                return startReturn;
            }
        }
        else
        {
            typedef TTAG(MPFAOTwoPTestProblem) ProblemTypeTag;
            int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
            std::cout<<"######################################################\n";
            std::cout<<"No model type specified\n";
            std::cout<<"Default to finite volume MPFA o-method model\n";
            return startReturn;
        }
    }
    else
    {
        typedef TTAG(MPFAOTwoPTestProblem) ProblemTypeTag;
        int startReturn =  Dumux::start<ProblemTypeTag>(argc, argv);
        std::cout<<"######################################################\n";
        std::cout<<s<<" is not a valid model type specification!\n";
        std::cout<<"Default to finite volume MPFA o-method model\n";
        return startReturn;
    }
}

#else

#warning You need to have ALUGrid installed to run this test

#include <iostream>

int main()
{
    std::cerr << "You need to have ALUGrid installed to run this test\n";
    return 1;
}

#endif
