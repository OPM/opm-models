// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Katherina Baber, Klaus Mosthaf                    *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
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
 * \brief Contains the quantities which are constant within a
 *        finite volume in the compositional Stokes model.
 */
#ifndef DUMUX_STOKES2C_VOLUME_VARIABLES_HH
#define DUMUX_STOKES2C_VOLUME_VARIABLES_HH

#include <dumux/freeflow/stokes/stokesvolumevariables.hh>
#include "stokes2cproperties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxVolumeVariables
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-component Stokes box model.
 */
template <class TypeTag>
class Stokes2cVolumeVariables : public StokesVolumeVariables<TypeTag>
{
    typedef StokesVolumeVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cIndices) Indices;

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { compIdx = GET_PROP_VALUE(TypeTag, StokesComponentIndex) };
    enum { transportIdx = Indices::transportIdx };

    static_assert(FluidSystem::numComponents == 2,
                  "Only fluid systems with two components are supported by the stokes2c model");
    static_assert(0 <= phaseIdx && phaseIdx < FluidSystem::numPhases,
                  "Invalid phase index");
    static_assert(0 <= compIdx && compIdx < FluidSystem::numComponents,
                  "Invalid component index");

public:
    /*!
     * \copydoc BoxVolumeVariables::update()
     */
    void update(const ElementContext &elemCtx, int scvIdx, int timeIdx)
    {
        const auto &priVars = elemCtx.primaryVars(scvIdx, timeIdx);

        Scalar X1 = priVars[transportIdx];
        Scalar X2 = 1.0 - X1;

        // calculate average molar mass of the gas phase
        Scalar M1 = FluidSystem::molarMass(compIdx);
        Scalar M2 = FluidSystem::molarMass(1 - compIdx);
        Scalar avgMolarMass = M1*M2/(M2 + X2*(M1 - M2));

        // convert mass to mole fractions and set the fluid state
        this->fluidState_.setMoleFraction(phaseIdx, compIdx, X1*avgMolarMass/M1);
        this->fluidState_.setMoleFraction(phaseIdx, 1 - compIdx, X2*avgMolarMass/M2);

        ParentType::update(elemCtx, scvIdx, timeIdx);

        // Second instance of a parameter cache.
        // This should be avoided if possible.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(this->fluidState_);

        diffCoeff_ = FluidSystem::binaryDiffusionCoefficient(this->fluidState_,
                                                             paramCache,
                                                             phaseIdx,
                                                             compIdx,
                                                             1 - compIdx);

        Valgrind::CheckDefined(diffCoeff_);
    };

    /*!
     * \brief Returns the binary (mass) diffusion coefficient
     */
    Scalar diffusionCoeff() const
    { return diffCoeff_; }

protected:
    Scalar diffCoeff_; //!< Binary diffusion coefficient
};

} // end namespace

#endif
