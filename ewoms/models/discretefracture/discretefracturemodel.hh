// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::DiscreteFractureModel
 */
#ifndef EWOMS_DISCRETE_FRACTURE_MODEL_HH
#define EWOMS_DISCRETE_FRACTURE_MODEL_HH

#include "discretefractureproperties.hh"

#include <ewoms/models/immiscible/immisciblemodel.hh>

#include <string>

namespace Ewoms {
/*!
 * \ingroup DiscreteFractureVcfvModel
 * \brief A fully-implicit multi-phase flow model which assumes
 *        immiscibility of the phases and is able to include fractures
 *        in the domain.
 *
 * This model implements multi-phase flow of \f$M > 0\f$ immiscible
 * fluids \f$\alpha\f$. It also can consider edges of the
 * computational grid as fractures i.e. as a porous medium with
 * different higher permeability than the rest of the domain.
 *
 * \todo So far, the discrete fracture model only works for 2D grids
 *       and without energy. Also only the Darcy velocity model is
 *       supported for the fractures.
 *
 * \sa ImmiscibleModel
 */
template <class TypeTag>
class DiscreteFractureModel : public ImmiscibleModel<TypeTag>
{
    typedef ImmiscibleModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

public:
    DiscreteFractureModel(Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible VCVF
     * discretization.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VcfvVtkDiscreteFractureModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc VcfvModel::name
     */
    const char *name() const
    { return "discretefracture"; }

protected:
    friend class VcfvModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        this->vtkOutputModules_.push_back(
            new Ewoms::VcfvVtkDiscreteFractureModule<TypeTag>(this->problem_));
    }
};
} // namespace Ewoms

#include "discretefracturepropertydefaults.hh"

#endif
