// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
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
#include <ewoms/common/start.hh>

int main(int argc, char** argv)
{
/*
    typedef TTAG(Pressure) TypeTag;
    std::cout << "Children: "
              << Dune::className<
                  TypeTag::ChildrenTuple
                     >() << "\n";
    std::cout << "EffTypeTag: "
              << Dune::className<
                  Ewoms::Properties::GetProperty<TypeTag, Ewoms::Properties::PTag::LinearSolverBackend>::GetEffectiveTypeTag_<TypeTag>::type
                     >() << "\n";
    std::cout << "Splices: "
              << Dune::className<
                  Ewoms::Properties::PTag::Splices<TypeTag>::tuple
                     >() << "\n";
*/

    typedef TTAG(TestProblemOneP) ProblemTypeTag;

    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
