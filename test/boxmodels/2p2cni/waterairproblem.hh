// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Klaus Mosthaf                                     *
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
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */
#ifndef DUMUX_WATER_AIR_PROBLEM_HH
#define DUMUX_WATER_AIR_PROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>

#include <dumux/boxmodels/2p2cni/2p2cnimodel.hh>

#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>

#include <dumux/material/heatconduction/somerton.hh>

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class WaterAirProblem;

namespace Properties
{
#if !ISOTHERMAL
NEW_TYPE_TAG(WaterAirProblem, INHERITS_FROM(BoxTwoPTwoCNI));
#else
NEW_TYPE_TAG(WaterAirProblem, INHERITS_FROM(BoxTwoPTwoC));
#endif

// Set the grid type
SET_PROP(WaterAirProblem, Grid)
{
    typedef Dune::YaspGrid<2> type;
};

// Set the problem property
SET_PROP(WaterAirProblem, Problem)
{
    typedef Dumux::WaterAirProblem<TypeTag> type;
};
// Set the material Law

SET_PROP(WaterAirProblem, MaterialLaw)
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
SET_PROP(WaterAirProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Dumux::Somerton<FluidSystem, Scalar> type;
};


// Set the wetting phase
SET_TYPE_PROP(WaterAirProblem, FluidSystem, 
              Dumux::FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar),
                                         /*complexRelations=*/true>);

// Enable gravity
SET_BOOL_PROP(WaterAirProblem, EnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(WaterAirProblem, NumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(WaterAirProblem, NewtonWriteConvergence, false);
}


/*!
 * \ingroup TwoPTwoCNIModel
 * \ingroup BoxTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium. During
 *        buoyancy driven upward migration the gas passes a high
 *        temperature area.
 *
 * The domain is sized 40 m times 40 m. The rectangular area with the
 * increased temperature (380 K) starts at (20 m, 5 m) and ends at (30
 * m, 35 m).
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top and on the bottom of the domain, while dirichlet conditions
 * apply on the left and the right boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the bottom boundary from 15 m to 25 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * This problem uses the \ref TwoPTwoCNIModel.
 *
 * This problem should typically be simulated for 300000 s.
 * A good choice for the initial time step size is 1000 s.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2p2cni -parameterFile test_2p2cni.input</tt>
 *  */
template <class TypeTag >
class WaterAirProblem 
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,
        numComponents = FluidSystem::numComponents,

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif
       
        // component indices
        N2Idx = FluidSystem::N2Idx,
        H2OIdx = FluidSystem::H2OIdx,

        // phase indices
        lPhaseIdx = FluidSystem::lPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        // equation indices
        contiN2EqIdx = Indices::conti0EqIdx + N2Idx,
        contiH2OEqIdx = Indices::conti0EqIdx + H2OIdx,

        // Phase State
        lPhaseOnly = Indices::lPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    WaterAirProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        maxDepth_ = 1000.0; // [m]
        eps_ = 1e-6;

        FluidSystem::init();

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
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "waterair"; }

#if ISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index (SCV index)
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    template <class Context>
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return 273.15 + 10; // -> 10°C
    };
#endif

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
     * \brief Returns the specific heat capacity \f$[J/(m^3 K)]\f$ of the rock matrix.
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
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineHeatCondParams_;
        return coarseHeatCondParams_;
    }


    template <class Context>
    void source(RateVector &values,
                const Context &context, int spaceIdx, int timeIdx) const
    {
        values = 0;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex for which the boundary type is set
     */
    template <class Context>
    void boundaryTypes(BoundaryTypes &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        if(globalPos[0] > 40 - eps_ || globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

#if !ISOTHERMAL
        values.setDirichlet(temperatureIdx, energyEqIdx);
#endif
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex for which the boundary type is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    template <class Context>
    void dirichlet(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    template <class Context>
    void neumann(RateVector &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        values = 0;

        // negative values for injection
        if (globalPos[0] > 15 && globalPos[0] < 25 &&
            globalPos[1] < eps_)
        {
            values[contiN2EqIdx] = -1e-3; // [kg/(s m^2)]
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
           const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        initial_(values, globalPos);

#if !ISOTHERMAL
            if (globalPos[0] > 20 && globalPos[0] < 30 && globalPos[1] < 30)
               values[temperatureIdx] = 380;
#endif
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;

        values.setPhasePresence(lPhaseOnly);

        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        values[switchIdx] = 0.0;
#if !ISOTHERMAL
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[1])*0.03;
#endif
    }

private:
    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {            
        Scalar lambdaGranite = 2.8; // [W / (K m)]

        // create a Fluid state which has all phases present
        Dumux::ImmiscibleFluidState<Scalar, FluidSystem> fs;
        fs.setTemperature(293.15);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fs.setPressure(phaseIdx, 1.0135e5);
        }

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar rho = FluidSystem::density(fs, paramCache, phaseIdx);
            fs.setDensity(phaseIdx, rho);
        }

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            Scalar lambdaSaturated;
            if (FluidSystem::isLiquid(phaseIdx)) {
                Scalar lambdaFluid =
                    FluidSystem::thermalConductivity(fs, paramCache, phaseIdx);
                lambdaSaturated = std::pow(lambdaGranite, (1-poro)) * std::pow(lambdaFluid, poro);
            }
            else
                lambdaSaturated = std::pow(lambdaGranite, (1-poro));
            
            params.setFullySaturatedLambda(phaseIdx, lambdaSaturated);
            if (!FluidSystem::isLiquid(phaseIdx))
                params.setVacuumLambda(lambdaSaturated);
        }
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

    Scalar maxDepth_;
    Scalar eps_;
};
} //end namespace

#endif
