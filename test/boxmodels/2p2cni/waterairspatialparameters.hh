// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
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
 * \brief Definition of the spatial parameters for the water-air problem.
 */
#ifndef DUMUX_WATER_AIR_SPATIAL_PARAMETERS_HH
#define DUMUX_WATER_AIR_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/Mp/2padapter.hh>

#include <dumux/material/heatconduction/somerton.hh>

#include <dumux/boxmodels/2p2cni/2p2cnimodel.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class WaterAirSpatialParameters;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(WaterAirSpatialParameters);

// Set the spatial parameters
SET_TYPE_PROP(WaterAirSpatialParameters, SpatialParameters, Dumux::WaterAirSpatialParameters<TypeTag>);

// Set the material Law
SET_PROP(WaterAirSpatialParameters, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> EffMaterialLaw;

    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffMaterialLaw> TwoPMaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { lPhaseIdx = FluidSystem::lPhaseIdx };

public:
    typedef TwoPAdapter<lPhaseIdx, TwoPMaterialLaw> type;
};

// Set the heat conduction law
SET_PROP(WaterAirSpatialParameters, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Dumux::Somerton<FluidSystem::lPhaseIdx, Scalar> type;
};
}

/*!
 * \ingroup TwoPTwoCNIModel
 * \ingroup BoxTestProblems
 * \brief Definition of the spatial parameters for the water-air problem
 */
template<class TypeTag>
class WaterAirSpatialParameters : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;
    enum {
        lPhaseIdx = FluidSystem::lPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> Vector;

    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;

public:
    /*!
     * \brief The constructor
     *
     * \param gv The grid view
     */
    WaterAirSpatialParameters(const GridView &gv)
        : ParentType(gv)
    {
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = 1e-13;
        coarseK_ = 1e-12;

        // porosities
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setSwr(0.2);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.2);
        coarseMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        fineMaterialParams_.setPe(1e4);
        coarseMaterialParams_.setPe(1e4);
        fineMaterialParams_.setLambda(2.0);
        coarseMaterialParams_.setLambda(2.0);

        // parameters for the somerton law of heat conduction
        computeHeatCondParams_(fineHeatCondParams_, finePorosity_);
        computeHeatCondParams_(coarseHeatCondParams_, coarsePorosity_);
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    const Scalar intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }


    /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
    * \param element The current finite element
    * \param fvElemGeom The current finite volume geometry of the element
    * \param scvIdx The index of the sub-control volume
    */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context, int spaceIdx, int timeIdx) const
    {
        return
            790 // specific heat capacity of granite [J / (kg K)]
            * 2700; // density of granite [kg/m^3]
    }

    /*!
     * \brief Return the parameter object for the heat conductivty law
     *        for a given position
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConducionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineHeatCondParams_;
        return coarseHeatCondParams_;
    }

private:
    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1-poro)) * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1-poro));

        params.setFullySaturatedLambda(gPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(lPhaseIdx, lambdaWet);
    }

    bool isFineMaterial_(const GlobalPosition &pos) const
    { return pos[dim-1] > layerBottom_; };

    Scalar fineK_;
    Scalar coarseK_;
    Scalar layerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams fineHeatCondParams_;
    HeatConductionLawParams coarseHeatCondParams_;
};

}

#endif
