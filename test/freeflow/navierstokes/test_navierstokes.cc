/*****************************************************************************
 *   Copyright (C) 2012 by Christoph Grueninger                              *
 *   Copyright (C) 2009-2012 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \brief Test for the isothermal Navier-Stokes box model; this test case is
 *        known as lid-driven cavity-flow in literature.
 */
#include "config.h"
#include "navierstokestestproblem.hh"
#include <dumux/common/start.hh>

int main(int argc, char** argv)
{
    typedef TTAG(NavierStokesTestProblem) ProblemTypeTag;
    return Dumux::start<ProblemTypeTag>(argc, argv);
}
