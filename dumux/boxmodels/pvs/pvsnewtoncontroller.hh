// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Melanie Darcis                                    *
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
 * \copydoc Dumux::PvsNewtonController
 */
#ifndef DUMUX_PVS_NEWTON_CONTROLLER_HH
#define DUMUX_PVS_NEWTON_CONTROLLER_HH

#include "pvsproperties.hh"

#include <dumux/boxmodels/common/boxnewtoncontroller.hh>

namespace Dumux {

/*!
 * \ingroup Newton
 * \ingroup PvsModel
 *
 * \brief A controller for the newton solver which is specific to the
 *        compositional multi-phase PVS box model.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class PvsNewtonController : public BoxNewtonController<TypeTag>
{
    typedef BoxNewtonController<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

public:
    PvsNewtonController(Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \copydoc NewtonController::newtonEndStep
     */
    void newtonEndStep(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ParentType::newtonEndStep(uCurrentIter, uLastIter);
        this->problem().model().switchPrimaryVars_();
    }

    /*!
     * \copydoc NewtonController::newtonConverged
     */
    bool newtonConverged()
    {
        if (this->problem().model().switched())
            return false;

        return ParentType::newtonConverged();
    }
};
}

#endif
