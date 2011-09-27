// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008,2009 by Andreas Lauser                               *
 *   Copyright (C) 2008,2009 by Melanie Darcis                               *
 *                                                                           *
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
 * \brief Contains the quantities to calculate the energy flux in the
 *        MpNc box model.
 */
#ifndef DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH
#define DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/boxmodels/mpnc/mpncproperties.hh>
#include <dumux/common/spline.hh>

namespace Dumux
{

template <class TypeTag, bool enableEnergy/*=false*/, bool kineticEnergyTransfer/*=false*/>
class MPNCFluxVariablesEnergy
{
    static_assert(!(kineticEnergyTransfer && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!kineticEnergyTransfer,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

public:
    MPNCFluxVariablesEnergy()
    {
    }

    void update(const ElementContext &elemCtx, int scvfIdx)
    {};
};

template <class TypeTag>
class MPNCFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true,  /*kineticEnergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        gPhaseIdx = FluidSystem::gPhaseIdx,
        lPhaseIdx = FluidSystem::lPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef Dune::FieldVector<CoordScalar, dimWorld>  Vector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    MPNCFluxVariablesEnergy()
    {}

    void update(const ElementContext &elemCtx, int scvfIdx)
    {
        const FVElementGeometry &fvElemGeom = elemCtx.fvElemGeom();

        // calculate temperature gradient using finite element
        // gradients
        Vector tmp(0.0);
        Vector temperatureGradient(0.);
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            tmp = fvElemGeom.subContVolFace[scvfIdx].grad[scvIdx];
            tmp *= elemCtx.volVars(scvIdx).fluidState().temperature();
            temperatureGradient += tmp;
        }

        // project the heat flux vector on the face's normal vector
        temperatureGradientNormal_ = 
            temperatureGradient * fvElemGeom.subContVolFace[scvfIdx].normal;

        lambdaPm_ = lumpedLambdaPm_(elemCtx, scvfIdx) ;
    }

    /*!
     * \brief The lumped / average conductivity of solid plus phases \f$[W/mK]\f$.
     */
    Scalar lambdaPm() const
    { return lambdaPm_; }

    /*!
     * \brief The normal of the gradient of temperature .
     */
    Scalar temperatureGradientNormal() const
    { return temperatureGradientNormal_; }

private:
    Scalar lumpedLambdaPm_(const ElementContext &elemCtx, int scvfIdx)
    {
        const FVElementGeometry &fvElemGeom = elemCtx.fvElemGeom();

        // arithmetic mean of the liquid saturation and the porosity
        const int i = fvElemGeom.subContVolFace[scvfIdx].i;
        const int j = fvElemGeom.subContVolFace[scvfIdx].j;
        
        const Scalar Sli = elemCtx.volVars(i).fluidState().saturation(lPhaseIdx);
        const Scalar Slj = elemCtx.volVars(j).fluidState().saturation(lPhaseIdx);
        
        const Scalar Sl = std::max<Scalar>(0.0, 0.5*(Sli + Slj));

        const Scalar lambda_solid =
            (elemCtx.volVars(i).thermalConductivitySolid() 
             + elemCtx.volVars(j).thermalConductivitySolid())
            / 2;
        
        //        const Scalar lambdaDry = 0.583; // W / (K m) // works, orig
        //        const Scalar lambdaWet = 1.13; // W / (K m) // works, orig
        
#warning "TODO: this is certainly not correct: it does not include porosity, and the mutable parameters probably need to be upwinded! Change lambda_dry and lambda_wet to spatial parameters? what about the case with more than 2 fluid phases, then?"
#warning "TODO/more: assumes thermal conductivity of the fluids to be independent of pressure, temperature and composition!"
        typename FluidSystem::MutableParameters mutParams; //dummy
        Scalar lambdaDry = 
            0.5 * (lambda_solid + FluidSystem::computeThermalConductivity(mutParams, gPhaseIdx) );
        Scalar lambdaWet =
            0.5 * (lambda_solid + FluidSystem::computeThermalConductivity(mutParams, lPhaseIdx));
        
        // the heat conductivity of the matrix. in general this is a
        // tensorial value, but we assume isotropic heat conductivity.
        // This is the Sommerton approach with lambdaDry =
        // lambdaSn100%.  Taken from: H. Class: "Theorie und
        // numerische Modellierung nichtisothermer Mehrphasenprozesse
        // in NAPL-kontaminierten poroesen Medien", PhD Thesis, University of
        // Stuttgart, Institute of Hydraulic Engineering, p. 57
        Scalar result;
        if (Sl < 0.1) {
            // regularization
            Dumux::Spline<Scalar> sp(0, 0.1, // x1, x2
                                    0, std::sqrt(0.1), // y1, y2
                                    5*0.5/std::sqrt(0.1), 0.5/std::sqrt(0.1)); // m1, m2
            result = lambdaDry + sp.eval(Sl)*(lambdaWet - lambdaDry);
        }
        else
            result = lambdaDry + std::sqrt(Sl)*(lambdaWet - lambdaDry);

        return result;
    }

    Scalar lambdaPm_;
    Scalar temperatureGradientNormal_;
};

} // end namepace

#endif
