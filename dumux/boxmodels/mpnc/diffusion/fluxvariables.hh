// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief This file contains the diffusion module for the flux data of
 *        the fully coupled two-phase N-component model
 */
#ifndef DUMUX_MPNC_DIFFUSION_FLUX_VARIABLES_HH
#define DUMUX_MPNC_DIFFUSION_FLUX_VARIABLES_HH

#include <dune/common/fvector.hh>

#include "../mpncproperties.hh"

namespace Dumux {

template<class TypeTag, bool enableDiffusion>
class MPNCFluxVariablesDiffusion
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        lPhaseIdx = FluidSystem::lPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx
    };

    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;



public:
    MPNCFluxVariablesDiffusion()
    {}

    void update(const ElementContext &elemCtx, int scvfIdx)
    {
        int i = elemCtx.fvElemGeom().subContVolFace[scvfIdx].i;
        int j = elemCtx.fvElemGeom().subContVolFace[scvfIdx].j;

        const VolumeVariables &volVarsI = elemCtx.volVars(i);
        const VolumeVariables &volVarsJ = elemCtx.volVars(j);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFrac_[phaseIdx][compIdx]  = volVarsI.fluidState().moleFraction(phaseIdx, compIdx);
                moleFrac_[phaseIdx][compIdx] += volVarsJ.fluidState().moleFraction(phaseIdx, compIdx);
                moleFrac_[phaseIdx][compIdx] /= 2;
            }
        }

        // initialize the diffusion coefficients to zero
        for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
            porousDiffCoeffL_[compIIdx] = 0.0;
            for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
                porousDiffCoeffG_[compIIdx][compJIdx] = 0.0;
        }
        
        // update the concentration gradients using two-point
        // gradients
        const GlobalPosition &normal = elemCtx.fvElemGeom().subContVolFace[scvfIdx].normal;

        GlobalPosition tmp = elemCtx.pos(j);
        tmp -= elemCtx.pos(i);
        Scalar dist = tmp.two_norm();
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // concentration gradients
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                moleFracGrad_[phaseIdx][compIdx] = normal;
                moleFracGrad_[phaseIdx][compIdx]
                    *=
                    (volVarsJ.fluidState().moleFraction(phaseIdx, compIdx) -
                     volVarsI.fluidState().moleFraction(phaseIdx, compIdx))
                    / (dist * normal.two_norm());
            }
        }
    

        // calculate the diffusion coefficients at the integration
        // point in the porous medium
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // make sure to only calculate diffusion coefficents
            // for phases which exist in both finite volumes
            if (volVarsI.fluidState().saturation(phaseIdx) <= 1e-4 ||
                volVarsJ.fluidState().saturation(phaseIdx) <= 1e-4)
            {
                continue;
            }

            // reduction factor for the diffusion coefficients in the
            // porous medium in nodes i and j. this is the tortuosity
            // times porosity times phase saturation at the nodes i
            // and j
            //
            // TODO (?): move this calculation to the soil (possibly
            // that's a bad idea, though)
            Scalar red_i =
                volVarsI.fluidState().saturation(phaseIdx)/volVarsI.porosity() *
                std::pow(volVarsI.porosity() * volVarsI.fluidState().saturation(phaseIdx), 7.0/3);
            Scalar red_j =
                volVarsJ.fluidState().saturation(phaseIdx)/volVarsJ.porosity() *
                std::pow(volVarsJ.porosity() * volVarsJ.fluidState().saturation(phaseIdx), 7.0/3);

            if (phaseIdx == FluidSystem::lPhaseIdx) {
                // Liquid phase diffusion coefficients in the porous medium
                for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
                    // -> arithmetic mean
                    porousDiffCoeffL_[compIIdx]
                        = 1./2*(red_i * volVarsI.diffCoeff(lPhaseIdx, 0, compIIdx) +
                                red_j * volVarsJ.diffCoeff(lPhaseIdx, 0, compIIdx));
                }
            }
            else {
                // Gas phase diffusion coefficients in the porous medium
                for (int compIIdx = 0; compIIdx < numComponents; ++compIIdx) {
                    for (int compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
                        // -> arithmetic mean
                        porousDiffCoeffG_[compIIdx][compJIdx]
                            = 1./2*(red_i * volVarsI.diffCoeff(gPhaseIdx, compIIdx, compJIdx) +
                                    red_j * volVarsJ.diffCoeff(gPhaseIdx, compIIdx, compJIdx));
                    }
                }
            }
        }
    };

    Scalar porousDiffCoeffL(int compIdx) const
    {
        // TODO: tensorial diffusion coefficients
        return porousDiffCoeffL_[compIdx];
    };

    Scalar porousDiffCoeffG(int compIIdx, int compJIdx) const
    {
        // TODO: tensorial diffusion coefficients
        return porousDiffCoeffG_[compIIdx][compJIdx];
    };

    Scalar moleFraction(int phaseIdx,
                    int compIdx) const
    {
        return moleFrac_[phaseIdx][compIdx];
    };

    const GlobalPosition &moleFracGrad(int phaseIdx,
                                            int compIdx) const
    {
        return moleFracGrad_[phaseIdx][compIdx];
    };

protected:
    // the diffusion coefficients for the porous medium for the
    // liquid phase
    Scalar porousDiffCoeffL_[numComponents];

    // the diffusion coefficients for the porous medium for the
    // gas phase
    Scalar porousDiffCoeffG_[numComponents][numComponents];

    // the concentration gradients of all components in all phases
    GlobalPosition moleFracGrad_[numPhases][numComponents];

    // the mole fractions of each component at the integration point
    Scalar moleFrac_[numPhases][numComponents];
};


template<class TypeTag>
class MPNCFluxVariablesDiffusion<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    MPNCFluxVariablesDiffusion()
    {}

    void update(const ElementContext &elemCtx, int scvfIdx)
    {
    };
};

}

#endif // DUMUX_MPNC_DIFFUSION_FLUX_VARIABLES_HH
