// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
 *
 * \brief Base class for all models which use the one-phase,
 *        box model.
 */
#ifndef DUMUX_1P_MODEL_HH
#define DUMUX_1P_MODEL_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include "1plocalresidual.hh"
#include "1pproblem.hh"

namespace Dumux
{
/*!
 * \ingroup OnePBoxModel
 * \brief A single-phase, isothermal flow model using the box scheme.
 *
 * Single-phase, isothermal flow model, which solves the mass
 * continuity equation
 * \f[
 \phi \frac{\partial \varrho}{\partial t} + \text{div} (- \varrho \frac{\textbf K}{\mu} ( \textbf{grad}\, p -\varrho {\textbf g})) = q,
 * \f]
 * discretized using a vertex-centered finite volume (box) scheme as
 * spatial and the implicit Euler method as time discretization.  The
 * model supports compressible as well as incompressible fluids.
 */
template<class TypeTag >
class OnePBoxModel : public BoxModel<TypeTag>
{
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    enum { dim = GridView::dimension };

public:
    /*!
     * \brief Returns a string with the model's human-readable name
     */
    std::string name() const
    { return "1p"; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    std::string primaryVarName(int pvIdx) const
    { 
        assert(pvIdx == 0);
        std::ostringstream oss;
        oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    std::string eqName(int eqIdx) const
    { 
        assert(eqIdx == 0);
        std::ostringstream oss;
        oss << "conti_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        return oss.str();
    }

protected:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules available on all model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
    };
};
}

#include "1ppropertydefaults.hh"

#endif
