// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Benjamin Faigle                                   *
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
 * \brief Source file of a test for the adaptive 2p2c model
 */
#include "config.h"

#include "test_adaptive2p2cproblem.hh"
#include <ewoms/common/start.hh>

// The main function using the standard start procedure
int main(int argc, char** argv)
{
#if HAVE_ALUGRID
    typedef TTAG(Adaptive2p2c) TypeTag;
    return Ewoms::start<TypeTag>(argc, argv);
#else
#warning ALUGrid needed for this test.
    std::cout << "ALUGrid needed for this test. Aborting." << std::endl;
#endif
}
