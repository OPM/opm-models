// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012-2013 by Markus Wolff                                 *
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
 * \copydoc Ewoms::DecoupledOnePCommonIndices
 */
#ifndef EWOMS_DECOUPLED_1P_INDICES_HH
#define EWOMS_DECOUPLED_1P_INDICES_HH

namespace Ewoms
{
/*!
 * \ingroup OnePhase
 */
// \{

/*!
 * \brief The common indices for the 1-p models.
 */
struct DecoupledOnePCommonIndices
{
    // Formulations
    static const int pressureEqIdx = 0;//!< Index of the pressure equation
};

// \}
} // namespace Ewoms


#endif
