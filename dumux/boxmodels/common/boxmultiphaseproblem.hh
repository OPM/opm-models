// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 * \ingroup BoxModel
 *
 * \brief The base class for the problems of box models which deal
 *        with a multi-phase flow through a porous medium.
 */
#ifndef DUMUX_BOX_MULTI_PHASE_PROBLEM_HH
#define DUMUX_BOX_MULTI_PHASE_PROBLEM_HH

#include "boxporousproblem.hh"

#include <dumux/common/math.hh>

namespace Dumux {
namespace Properties {
NEW_PROP_TAG(MaterialLawParams);
}

/*!
 * \ingroup BoxModel
 *
 * \brief The base class for the problems of box models which deal
 *        with a multi-phase flow through a porous medium.
 */
template<class TypeTag>
class BoxMultiPhaseProblem : public BoxPorousProblem<TypeTag>
{
    typedef Dumux::BoxPorousProblem<TypeTag> ParentType;
    
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    
public:
    BoxMultiPhaseProblem(TimeManager &timeManager, const GridView &gv)
        : ParentType(timeManager, gv)
    { }

    template <class Context>
    const MaterialLawParams& 
    materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Problem::materialLawParams()");
    }
};

} // namespace Dumux

#endif
