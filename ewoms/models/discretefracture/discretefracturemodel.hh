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
#include "discretefractureprimaryvariables.hh"
#include "discretefracturevolumevariables.hh"
#include "discretefracturefluxvariables.hh"
#include "discretefracturelocalresidual.hh"

#include <ewoms/models/immiscible/immisciblemodel.hh>
#include <ewoms/vtk/vtkdiscretefracturemodule.hh>

#include <string>

namespace Ewoms {
template <class TypeTag>
class DiscreteFractureModel;
}

namespace Opm {
namespace Properties {
//! The generic type tag for problems using the immiscible multi-phase model
NEW_TYPE_TAG(DiscreteFractureModel, INHERITS_FROM(ImmiscibleTwoPhaseModel, VtkDiscreteFracture));

//! The class for the model
SET_TYPE_PROP(DiscreteFractureModel, Model, Ewoms::DiscreteFractureModel<TypeTag>);

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(DiscreteFractureModel, LocalResidual, Ewoms::DiscreteFractureLocalResidual<TypeTag>);

// The type of the base base class for actual problems.
// TODO!?
//SET_TYPE_PROP(DiscreteFractureModel BaseProblem, DiscreteFractureBaseProblem<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(DiscreteFractureModel, PrimaryVariables,
              Ewoms::DiscreteFracturePrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(DiscreteFractureModel, VolumeVariables,
              Ewoms::DiscreteFractureVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(DiscreteFractureModel, FluxVariables, Ewoms::DiscreteFractureFluxVariables<TypeTag>);

//! For the discrete fracture model, we need to use two-point flux
//! appoximation or it will converge very poorly
SET_BOOL_PROP(DiscreteFractureModel, UseTwoPointGradients, true);

// The volume variable cache cannot be used by the discrete fracture
// model, because the volume variables of a control sub-control volume
// are not identical to the volume variables of the other volume
// variables of the same of the same control volume. This is because
// the fracture properties (volume, permeability, etc) are specific
// for each sub control volume...
SET_BOOL_PROP(DiscreteFractureModel, EnableVolumeVariablesCache, false);
}} // namespace Properties, Opm

namespace Ewoms {
/*!
 * \ingroup DiscreteFractureModel
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
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkDiscreteFractureModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "discretefracture"; }

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        this->vtkOutputModules_.push_back(
            new Ewoms::VtkDiscreteFractureModule<TypeTag>(this->problem_));
    }
};
} // namespace Ewoms

#endif
