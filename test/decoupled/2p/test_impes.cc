// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Philipp Nuske                                     *
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Benjamin Faigle                                   *
 *   Copyright (C) 2010 by Bernd Flemisch                                    *
 *   Copyright (C) 2010-2012 by Markus Wolff                                 *
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
 * \brief Source file of a test for the sequential 2p model
 */
#include "config.h"

#include "test_impesproblem.hh"
#include <dumux/common/start.hh>

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    typedef TTAG(IMPESTestProblem) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv);
}
