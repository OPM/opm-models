/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Test for the two-phase flow with discrete fracture-matrix (discrete
 *        fracture model) and box model scheme.
 */
#include "config.h"

#include "2pdfmtestproblem.hh"
#include <dumux/common/start.hh>

int main(int argc, char** argv)
{
    typedef TTAG(TwoPDFMTestProblem) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv, usage);
}
