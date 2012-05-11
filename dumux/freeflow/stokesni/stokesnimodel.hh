// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \brief Adaption of the BOX scheme to the non-isothermal
 *        compositional stokes model (with two components).
 */
#ifndef DUMUX_STOKES_NI_MODEL_HH
#define DUMUX_STOKES_NI_MODEL_HH

#include "stokesnilocalresidual.hh"
#include "stokesniproperties.hh"

#include <dumux/freeflow/stokes/stokesmodel.hh>

#include <dune/common/fvector.hh>

#include <sstream>

namespace Dumux {
/*!
 * \ingroup BoxStokesNIModel
 * \brief Adaption of the BOX scheme to the non-isothermal compositional Stokes model.
 *
 * This model implements a non-isothermal two-component Stokes flow of a fluid
 * solving a momentum balance, a mass balance, a conservation equation for one component,
 * and one balance equation for the energy.
 *
 * Momentum Balance:
 * \f[
\frac{\partial \left(\varrho_g {\boldsymbol{v}}_g\right)}{\partial t}
+ \boldsymbol{\nabla} \boldsymbol{\cdot} \left(p_g {\bf {I}}
- \mu_g \left(\boldsymbol{\nabla} \boldsymbol{v}_g
+ \boldsymbol{\nabla} \boldsymbol{v}_g^T\right)\right)
- \varrho_g {\bf g} = 0,
 * \f]
 *
 * Mass balance equation:
 * \f[
\frac{\partial \varrho_g}{\partial t} + \boldsymbol{\nabla}\boldsymbol{\cdot}\left(\varrho_g {\boldsymbol{v}}_g\right) - q_g = 0
 * \f]
 *
 * Component mass balance equation:
 * \f[
 \frac{\partial \left(\varrho_g X_g^\kappa\right)}{\partial t}
 + \boldsymbol{\nabla} \boldsymbol{\cdot} \left( \varrho_g {\boldsymbol{v}}_g X_g^\kappa
 - D^\kappa_g \varrho_g \boldsymbol{\nabla} X_g^\kappa \right)
 - q_g^\kappa = 0
 * \f]
 *
 * Energy balance equation:
 * \f[
\frac{\partial (\varrho_g  u_g)}{\partial t}
+ \boldsymbol{\nabla} \boldsymbol{\cdot} \varrho_g h_g {\boldsymbol{v}}_g
- \lambda_g \boldsymbol{\nabla} T - q_T = 0
 * \f]
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 */
template<class TypeTag>
class StokesNIModel : public StokesModel<TypeTag>
{
    typedef StokesModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    
public:
    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::ostringstream oss;
        if (pvIdx == Indices::temperatureIdx) {
            oss << "temperature";
            return oss.str();
        }

        return ParentType::primaryVarName(pvIdx);
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    std::string eqName(int eqIdx) const
    {
        std::ostringstream oss;
        if (eqIdx == Indices::energyEqIdx) {
            oss << "energy";
            return oss.str();
        }

        return ParentType::eqName(eqIdx);
    }
};

}

#include "stokesnipropertydefaults.hh"

#endif
