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
 * \copydoc Ewoms::PvsNewtonMethod
 */
#ifndef EWOMS_PVS_NEWTON_METHOD_HH
#define EWOMS_PVS_NEWTON_METHOD_HH

#include "pvsproperties.hh"

#include <ewoms/disc/vcfv/vcfvnewtonmethod.hh>

namespace Ewoms {

/*!
 * \ingroup Newton
 * \ingroup PvsModel
 *
 * \brief A newton solver which is specific to the compositional
 *        multi-phase PVS VCVF discretization.
 */
template <class TypeTag>
class PvsNewtonMethod : public VcfvNewtonMethod<TypeTag>
{
    typedef VcfvNewtonMethod<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

public:
    PvsNewtonMethod(Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \copydoc NewtonMethod::converged
     */
    bool converged() const
    {
        if (this->problem().model().switched())
            return false;

        return ParentType::converged();
    }

protected:
    friend class NewtonMethod<TypeTag>;
    friend class VcfvNewtonMethod<TypeTag>;

    /*!
     * \copydoc NewtonMethod::endIteration_
     */
    void endIteration_(SolutionVector &uCurrentIter,
                       const SolutionVector &uLastIter)
    {
        ParentType::endIteration_(uCurrentIter, uLastIter);
        this->problem().model().switchPrimaryVars_();
    }
};
}

#endif
