// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Melanie Darcis                                    *
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
 * \brief Nonisothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 */
#ifndef DUMUX_INJECTION_PROBLEM_2PNI_HH
#define DUMUX_INJECTION_PROBLEM_2PNI_HH

#include <dumux/boxmodels/2pni/2pnimodel.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>
#include <dumux/material/heatconduction/somerton.hh>
#include <dumux/common/cubegridcreator.hh>

#include <dune/common/fvector.hh>

namespace Dumux {

template <class TypeTag>
class InjectionProblem2PNI;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(BoxTwoPNI));

// declare the properties specific for the non-isothermal immiscible
// injection problem
NEW_PROP_TAG(MaxDepth);

NEW_PROP_TAG(FluidSystemNumTemperature);
NEW_PROP_TAG(FluidSystemTemperatureLow);
NEW_PROP_TAG(FluidSystemTemperatureHigh);
NEW_PROP_TAG(FluidSystemNumPressure);
NEW_PROP_TAG(FluidSystemPressureLow);
NEW_PROP_TAG(FluidSystemPressureHigh);

// Set the grid type
SET_TYPE_PROP(InjectionProblem2PNI, 
              Grid,
              Dune::YaspGrid<2>);

// set the GridCreator property
SET_TYPE_PROP(InjectionProblem2PNI,
              GridCreator,
              CubeGridCreator<TypeTag>);

// Set the problem property
SET_TYPE_PROP(InjectionProblem2PNI,
              Problem,
              Dumux::InjectionProblem2PNI<TypeTag>);

#if 1
// Use the same fluid system as the 2p2c injection problem
SET_TYPE_PROP(InjectionProblem2PNI,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);
#else
// Set the wetting phase
SET_PROP(InjectionProblem2PNI, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(InjectionProblem2PNI, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::GasPhase<Scalar, Dumux::N2<Scalar> > type;
};
#endif

// Set the material Law
SET_PROP(InjectionProblem2PNI, MaterialLaw)
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
SET_PROP(InjectionProblem2PNI, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Dumux::Somerton<FluidSystem, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem2PNI, EnableGravity, true);

// write convergence behaviour to disk?
SET_BOOL_PROP(InjectionProblem2PNI, NewtonWriteConvergence, false);

//SET_TYPE_PROP(InjectionProblem2PNI, LinearSolver, SuperLUBackend<TypeTag>);

// define the properties specific for the non-isothermal immiscible
// injection problem
SET_SCALAR_PROP(InjectionProblem2PNI, GridSizeX, 60.0);
SET_SCALAR_PROP(InjectionProblem2PNI, GridSizeY, 40.0);
SET_SCALAR_PROP(InjectionProblem2PNI, GridSizeZ, 0.0);

SET_INT_PROP(InjectionProblem2PNI, GridCellsX, 24);
SET_INT_PROP(InjectionProblem2PNI, GridCellsY, 26);
SET_INT_PROP(InjectionProblem2PNI, GridCellsZ, 0);
SET_SCALAR_PROP(InjectionProblem2PNI, MaxDepth, 2700);

SET_INT_PROP(InjectionProblem2PNI, FluidSystemNumTemperature, 100);
SET_SCALAR_PROP(InjectionProblem2PNI, FluidSystemTemperatureLow, 290.0);
SET_SCALAR_PROP(InjectionProblem2PNI, FluidSystemTemperatureHigh, 500.0);
SET_INT_PROP(InjectionProblem2PNI, FluidSystemNumPressure, 200);
SET_SCALAR_PROP(InjectionProblem2PNI, FluidSystemPressureLow, 2e7);
SET_SCALAR_PROP(InjectionProblem2PNI, FluidSystemPressureHigh, 3e7);
}

/*!
 * \ingroup TwoPNIModel
 * \ingroup BoxTestProblems
 * \brief Nonisothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 *
 * The domain is sized 40 m times 40 m. The rectangular area with the increased temperature (380 K)
 * starts at (20 m, 5 m) and ends at (30 m, 35 m)
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
 * This problem uses the \ref TwoPNIModel.
 *
 * This problem should typically be simulated for 300000 seconds.
 * A good choice for the initial time step size is 1000 seconds.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2pni -parameterFile test_2pni.input</tt>
 */
template<class TypeTag>
class InjectionProblem2PNI 
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename HeatConductionLaw::Params HeatConductionLawParams;

    enum {
        numPhases = FluidSystem::numPhases,

        SnIdx = Indices::SnIdx,
        pwIdx = Indices::pwIdx,
        conti0EqIdx = Indices::conti0EqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        contiNEqIdx = Indices::conti0EqIdx + nPhaseIdx,

        temperatureIdx = Indices::temperatureIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;


public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem2PNI(TimeManager &timeManager)
        : ParentType(timeManager, 
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 1e-6;

        maxDepth_ = GET_PARAM(TypeTag, Scalar, MaxDepth); // [m]

        // initialize the tables of the fluid system
        FluidSystem::init(GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureLow),
                          GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureHigh),
                          GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, NumTemperature),
                          GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureLow),
                          GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureHigh),
                          GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, NumPressure));
        
        layerBottom_ = 22.0;

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(1e-13);
        coarseK_ = this->toDimMatrix_(1e-12);

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
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
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
        return coarsePorosity_;
    }

    /*!
     * \brief Return the parameter object for the Brooks-Corey material law which depends on the position
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
    { return "injection2pni"; }

    template <class Context>
    void source(RateVector &values,
                const Context &context, int spaceIdx, int timeIdx) const
    { values = 0; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Evaluate the boundary conditions.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        if (onInlet_(pos)) {
            RateVector massRate(0.0);
            
            // inject air. negative values mean injection
            massRate[contiNEqIdx] = -1e-3; // kg/(s m^2)
                       
            // calculate the rate of enthalpy inflow
            Dumux::ImmiscibleFluidState<Scalar, FluidSystem> fs;          
            initialFluidState_(fs, context, spaceIdx, timeIdx);

            values.setMassRate(massRate);
            values.setEnthalpyRate(massRate[contiNEqIdx]
                                   * fs.enthalpy(nPhaseIdx));
           
        }
        else if (onLeftBoundary_(pos)) {
            Dumux::ImmiscibleFluidState<Scalar, FluidSystem> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);

            // impose a free flow (Dirichlet) boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else
            values.setNoFlow();
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
        ImmiscibleFluidState<Scalar, FluidSystem> fs;
        
        initialFluidState_(fs, context, spaceIdx, timeIdx);
        values[pwIdx] = fs.pressure(wPhaseIdx);
        values[SnIdx] = fs.saturation(nPhaseIdx);
        values[temperatureIdx] = fs.temperature(/*phaseIdx=*/0);
    }
    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->bboxMax()[1] - eps_; }

    bool onInlet_(const GlobalPosition &pos) const
    { return onRightBoundary_(pos) && 5 < pos[1] && pos[1] < 15; }

    bool inHighTemperatureRegion_(const GlobalPosition &pos) const
    { return (pos[0] > 20) && (pos[0] < 30) && (pos[1] > 5) && (pos[1] < 35); }

    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1-poro)) + std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1-poro));

        params.setFullySaturatedLambda(nPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(wPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    template <class FluidState, class Context>
    void initialFluidState_(FluidState &fs, 
                            const Context &context,
                            int spaceIdx,
                            int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);
        
        Scalar densityW = 1000.0;
        
        // set temperature for _all_ phases
        fs.setTemperature(283.0 + (maxDepth_ - pos[1])*0.03);
        if (inHighTemperatureRegion_(pos))
            fs.setTemperature(380);
                    
        // set the wetting phase saturation and pressure
        fs.setPressure(wPhaseIdx, 1e5 + (maxDepth_ - pos[1])*densityW*9.81);
        fs.setSaturation(wPhaseIdx, 1.0);
        
        // set the non-wetting phase saturation
        fs.setSaturation(nPhaseIdx, 0.0);
        
        // set the non-wetting phase saturation
        Scalar pC[numPhases];
        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);
        fs.setPressure(nPhaseIdx, fs.pressure(wPhaseIdx) + (pC[nPhaseIdx] - pC[wPhaseIdx]));

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            fs.setDensity(phaseIdx, FluidSystem::density(fs, paramCache, phaseIdx));
            fs.setEnthalpy(phaseIdx, FluidSystem::enthalpy(fs, paramCache, phaseIdx));
        }
    }
        
    bool isFineMaterial_(const GlobalPosition &pos) const
    { return pos[dim-1] > layerBottom_; };

    DimMatrix fineK_;
    DimMatrix coarseK_;
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
