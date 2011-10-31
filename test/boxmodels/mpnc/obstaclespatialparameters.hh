// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#ifndef DUMUX_OBSTACLE_SPATIAL_PARAMS_HH
#define DUMUX_OBSTACLE_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedlinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/heatconduction/somerton.hh>

#include <dumux/boxmodels/mpnc/mpncmodel.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mplinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class ObstacleSpatialParameters;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(ObstacleSpatialParameters);

// Set the spatial parameters
SET_TYPE_PROP(ObstacleSpatialParameters, SpatialParameters, Dumux::ObstacleSpatialParameters<TypeTag>);

// Set the material Law
SET_PROP(ObstacleSpatialParameters, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum {
        lPhaseIdx = FluidSystem::lPhaseIdx
    };
    // define the material law
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    //    typedef RegularizedBrooksCorey<Scalar> EffMaterialLaw;
        typedef RegularizedLinearMaterial<Scalar> EffMaterialLaw;
        typedef EffToAbsLaw<EffMaterialLaw> TwoPMaterialLaw;
public:
    typedef TwoPAdapter<lPhaseIdx, TwoPMaterialLaw> type;
};

// Set the heat conduction law
SET_PROP(ObstacleSpatialParameters, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Dumux::Somerton<FluidSystem::lPhaseIdx, Scalar> type;
};
}

/**
 * \ingroup MPNCModel
 * \ingroup BoxTestProblems
 * \brief Definition of the spatial params properties for the obstacle problem
 *
 */
template<class TypeTag>
class ObstacleSpatialParameters : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    enum {
        lPhaseIdx = FluidSystem::lPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename HeatConductionLaw::Params HeatConductionLawParams;

    ObstacleSpatialParameters(const GridView &gv)
        : ParentType(gv)
    {
        // intrinsic permeabilities
        coarseK_ = 1e-12;
        fineK_ = 1e-15;

        // the porosity
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setSwr(0.0);
        fineMaterialParams_.setSnr(0.0);
        coarseMaterialParams_.setSwr(0.0);
        coarseMaterialParams_.setSnr(0.0);

        // parameters for the linear law, i.e. minimum and maximum
        // pressures
        fineMaterialParams_.setEntryPC(0.0);
        coarseMaterialParams_.setEntryPC(0.0);
        fineMaterialParams_.setMaxPC(0.0);
        coarseMaterialParams_.setMaxPC(0.0);

        /*
        // entry pressures for Brooks-Corey
        fineMaterialParams_.setPe(5e3);
        coarseMaterialParams_.setPe(1e3);

        // Brooks-Corey shape parameters
        fineMaterialParams_.setLambda(2);
        coarseMaterialParams_.setLambda(2);
        */

        // parameters for the somerton law of heat conduction
        computeHeatCondParams_(fineHeatCondParams_, finePorosity_);
        computeHeatCondParams_(coarseHeatCondParams_, coarsePorosity_);

    }

    ~ObstacleSpatialParameters()
    {}

    /*!
     * \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     *
     * \param globalSol The current solution vector
     */
    void update(const SolutionVector &globalSol)
    {
    };

    /*!
     * \brief Returns the intrinsic permeability tensor.
     *
     * \param element       The current finite element
     * \param fvElemGeom    The current finite volume geometry of the element
     * \param scvIdx        The index sub-control volume where the
     *                      intrinsic permeability is given.
     */
    template <class Context>
    Scalar intrinsicPermeability(const Context &context, int localIdx) const
    {
        if (isFineMaterial_(context.pos(localIdx)))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param element     The finite element
     * \param fvElemGeom  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    template <class Context>
    Scalar porosity(const Context &context, int localIdx) const
    {
        const GlobalPosition &pos = context.pos(localIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     *
     * \param pos The global position of the sub-control volume.
     * \return the material parameters object
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context, int localIdx) const
    {
        const GlobalPosition &pos = context.pos(localIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \brief Returns the volumetric heat capacity \f$[J/m^3 K]\f$ of
     *        the rock matrix.
     *
     * Porosity is _not_ taken into account by this method. This is
     * only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context, int localIdx) const
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
    heatConducionParams(const Context &context, int localIdx) const
    {
        const GlobalPosition &pos = context.pos(localIdx);
        if (isFineMaterial_(pos))
            return fineHeatCondParams_;
        return coarseHeatCondParams_;
    }

private:
    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = pow(lambdaGranite, (1-poro)) * pow(lambdaWater, poro);
        Scalar lambdaDry = pow(lambdaGranite, (1-poro));

        params.setFullySaturatedLambda(gPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(lPhaseIdx, lambdaWet);
    }

    /*!
     * \brief Returns whether a given global position is in the
     *        fine-permeability region or not.
     */
    static bool isFineMaterial_(const GlobalPosition &pos)
    {
        return
            10 <= pos[0] && pos[0] <= 20 &&
            0 <= pos[1] && pos[1] <= 35;
    };

    Scalar coarseK_;
    Scalar fineK_;

    Scalar coarsePorosity_;
    Scalar finePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams fineHeatCondParams_;
    HeatConductionLawParams coarseHeatCondParams_;

};

}

#endif
