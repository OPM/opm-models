// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Base class for all models which use the single-phase,
 *        two-component box model
 *        Adaption of the BOX scheme to the one-phase two-component flow model.
 */

#ifndef DUMUX_ONEP_TWOC_MODEL_HH
#define DUMUX_ONEP_TWOC_MODEL_HH

#include "1p2cproperties.hh"
#include "1p2cproblem.hh"
#include "1p2clocalresidual.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModels
 * \defgroup OnePTwoCBoxModel One-phase Two-component box model
 */

/*!
 * \ingroup OnePTwoCBoxModel
 * \brief Adaption of the BOX scheme to the one-phase two-component flow model.
 *
 * This model implements an one-phase flow of an incompressible fluid, that consists of two components,
 * using a standard Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v_{D} = - \frac{\textbf K}{\mu}
 \left(\text{grad} p - \varrho {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the Property system.
 * By inserting this into the continuity equation, one gets
 \f[
 \Phi \frac{\partial \varrho}{\partial t} - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\text{grad}\, p - \varrho {\textbf g} \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components is described by the following equation:
 \f[
 \Phi \frac{ \partial \varrho x}{\partial t} - \text{div} \left( \varrho \frac{{\textbf K} x}{\mu} \left( \text{grad}\, p -
 \varrho {\textbf g} \right) + \varrho \tau \Phi D \text{grad} x \right) = q.
 \f]
 *
 * All equations are discretized using a fully-coupled vertex-centered
 * finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * The primary variables are the pressure \f$p\f$ and the mole fraction of dissolved component \f$x\f$.
 */

template<class TypeTag >
class OnePTwoCBoxModel : public BoxModel<TypeTag>
{
    typedef BoxModel<TypeTag> ParentType;

protected:
    friend class BoxModel<TypeTag>;
    
    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules meaninful for the model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
    };
};
}

#include "1p2cpropertydefaults.hh"

#endif
