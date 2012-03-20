// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \brief   This file contains the data which is required to calculate
 *          all fluxes of components over a face of a finite volume for
 *          the two-phase, two-component model.
 */
/*!
 * \ingroup ThreePThreeCModel
 */
#ifndef DUMUX_3P3C_FLUX_VARIABLES_HH
#define DUMUX_3P3C_FLUX_VARIABLES_HH

#include "3p3cproperties.hh"

#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>
#include <dumux/common/math.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class ThreePThreeCFluxVariables : public BoxMultiPhaseFluxVariables<TypeTag>
{
    typedef BoxMultiPhaseFluxVariables<TypeTag> MultiPhaseFluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        MultiPhaseFluxVariables::update(elemCtx, scvfIdx, timeIdx);

        calculateGradients_(elemCtx, scvfIdx, timeIdx);
        calculateDiffCoeffPM_(elemCtx, scvfIdx, timeIdx);
    };

    /*!
     * \brief The effective diffusion coefficient in the porous
     *        medium.
     */
    Scalar porousDiffusionCoefficient(int phaseIdx, int compIdx) const
    { return porousDiffCoeff_[phaseIdx][compIdx]; };

    /*!
     * \brief Returns the gradient of the mole fraction of a component in a phase.
     */
    const Vector &moleFracGrad(int phaseIdx, int compIdx) const
    { return moleFracGrad_[phaseIdx][compIdx]; };

    /*!
     * \brief Returns the molar density of a phase at the integration point
     */
    Scalar molarDensity(int phaseIdx) const
    { return molarDensity_[phaseIdx]; };

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int scvfIdx, int timeIdx)
    {
        // reset all gradients to 0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFracGrad_[phaseIdx][compIdx] = Scalar(0);
            }
        }

        // retrieve a reference to the structure for the sub control
        // volume face
        const auto &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];
        
        // calculate gradients
        Vector tmp;
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            // the fluid state at the current vertex
            const auto &fluidState = elemCtx.volVars(scvIdx, timeIdx).fluidState();

            // FE gradient at vertex scvIdx
            const Vector &feGrad = scvf.grad[scvIdx];

            // the concentration gradient of the components
            // component in the phases

            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                for (int compIdx = 0; compIdx < numComponents; compIdx++)
                {
                    tmp = feGrad;
                    tmp *= fluidState.moleFraction(phaseIdx, compIdx);
                    moleFracGrad_[phaseIdx][compIdx] += tmp;
                }
            }
        }
    }

    void calculateDiffCoeffPM_(const ElementContext &elemCtx,
                               int scvfIdx, int timeIdx)
    {
        const VolumeVariables &volVarsI = elemCtx.volVars(this->insideIdx(), timeIdx);
        const VolumeVariables &volVarsJ = elemCtx.volVars(this->outsideIdx(), timeIdx);

        const auto &fsI = volVarsI.fluidState();
        const auto &fsJ = volVarsJ.fluidState();

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            molarDensity_[phaseIdx] = (fsI.molarDensity(phaseIdx) + fsJ.molarDensity(phaseIdx))/2;

            // calculate tortuosity at the nodes i and j needed
            // for porous media diffusion coefficient
            Scalar tauI =
                std::pow(volVarsI.porosity() * fsI.saturation(phaseIdx), 7.0/3)
                / (volVarsI.porosity() * volVarsI.porosity());
            Scalar tauJ =
                std::pow(volVarsJ.porosity() * fsJ.saturation(phaseIdx), 7.0/3)
                / (volVarsJ.porosity() * volVarsJ.porosity());

            // Diffusion coefficient in the porous medium
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar Di = 
                    volVarsI.porosity()
                    * fsI.saturation(phaseIdx)
                    * tauI
                    * volVarsI.diffusionCoefficient(phaseIdx, compIdx);
                Scalar Dj =
                    volVarsJ.porosity()
                    * fsJ.saturation(phaseIdx)
                    * tauJ
                    * volVarsJ.diffusionCoefficient(phaseIdx, compIdx);

                porousDiffCoeff_[phaseIdx][compIdx] = harmonicMean(Di, Dj);
            }
        }
    }

    // molar densities
    Scalar molarDensity_[numPhases];

    // gradients
    Vector moleFracGrad_[numPhases][numComponents];

    // the diffusivity matrix for the porous medium
    Scalar porousDiffCoeff_[numPhases][numComponents];
};

} // end namepace

#endif
