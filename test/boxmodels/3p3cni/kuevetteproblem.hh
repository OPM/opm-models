// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
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
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 */
#ifndef DUMUX_KUEVETTE_PROBLEM_HH
#define DUMUX_KUEVETTE_PROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkerVanGen3p.hh>
#include <dumux/material/fluidmatrixinteractions/mp/3padapter.hh>
#include <dumux/material/heatconduction/somerton.hh>
#include <dumux/boxmodels/3p3cni/3p3cnimodel.hh>

namespace Dumux
{
template <class TypeTag>
class KuevetteProblem;

namespace Properties
{
NEW_TYPE_TAG(KuevetteProblem, INHERITS_FROM(BoxThreePThreeCNI));

// Set the grid type
SET_TYPE_PROP(KuevetteProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(KuevetteProblem, Problem, Dumux::KuevetteProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(KuevetteProblem, 
              FluidSystem,
              Dumux::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(KuevetteProblem, EnableGravity, true);

// Use central differences (backward -1, forward +1)
SET_INT_PROP(KuevetteProblem, NumericDifferenceMethod, 0);

// Write newton convergence
//SET_BOOL_PROP(KuevetteProblem, NewtonWriteConvergence, true);

// Set the maximum time step
SET_SCALAR_PROP(KuevetteProblem, MaxTimeStepSize, 60.);

// set newton relative tolerance
SET_SCALAR_PROP(KuevetteProblem, NewtonRelTolerance, 1e-6);


// Set the material Law
SET_PROP(KuevetteProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

    // define the three-phase material law
    typedef Dumux::ParkerVanGen3P<Scalar> ThreePLaw;
    
public:
    // wrap the three-phase law in an adaptor to make use the generic
    // material law API
    typedef Dumux::ThreePAdapter<wPhaseIdx, nPhaseIdx, gPhaseIdx, ThreePLaw> type;
};

// Set the heat conduction law
SET_PROP(KuevetteProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Dumux::Somerton<FluidSystem, Scalar> type;
};
}


/*!
 * \ingroup ThreePThreeCNIBoxModel
 * \ingroup BoxTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 *
 * The domain is a quasi-two-dimensional container (kuevette). Its dimensions
 * are 1.5 m x 0.74 m. The top and bottom boundaries are closed, the right
 * boundary is a Dirichlet boundary allowing fluids to escape. From the left,
 * an injection of a hot water-air mixture is applied (Neumann boundary condition
 * for the mass components and the enthalpy), aimed at remediating an initial 
 * NAPL (Non-Aquoeus Phase Liquid) contamination in the heterogeneous domain.
 * The contamination is initially placed partly into the coarse sand
 * and partly into a fine sand lense.
 *
 * This simulation can be varied through assigning different boundary conditions
 * at the left boundary as described in Class (2001):
 * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in
 * NAPL-kontaminierten por"osen Medien, Dissertation, Eigenverlag des Instituts
 * f"ur Wasserbau
 *
 * This problem uses the \ref ThreePThreeCNIModel
 *
 * To see the basic effect and the differences to scenarios with pure steam or
 * pure air injection, it is sufficient to simulated for about 2-3 hours (10000 s).
 * Complete remediation of the domain requires much longer (about 10 days simulated time).
 * To adjust the simulation time it is necessary to edit the file test_3p3cni.input
 *
 * To run the simulation execute:
 *
 * <tt>./test_3p3cni -parameterFile test_3p3cni.input</tt>
 *  */
template <class TypeTag >
class KuevetteProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLawParams) HeatConductionLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, ThreePThreeCNIIndices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,

        pressure0Idx = Indices::pressure0Idx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,

        // Phase State
        threePhases = Indices::threePhases,
        wgPhaseOnly = Indices::wgPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    KuevetteProblem(TimeManager &timeManager)
        : ParentType(timeManager, 
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(1e-6)
    {
        FluidSystem::init();

        // intrinsic permeabilities
        fineK_ = 6.28e-12;
        coarseK_ = 9.14e-10;

        // porosities
        finePorosity_ = 0.42;
        coarsePorosity_ = 0.42;

        // residual saturations
        fineMaterialParams_.setSwr(0.12);
        fineMaterialParams_.setSwrx(0.12);
        fineMaterialParams_.setSnr(0.07);
        fineMaterialParams_.setSgr(0.01);
        coarseMaterialParams_.setSwr(0.12);
        coarseMaterialParams_.setSwrx(0.12);
        coarseMaterialParams_.setSnr(0.07);
        coarseMaterialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.005);
        fineMaterialParams_.setVgN(4.0);
        coarseMaterialParams_.setVgN(4.0);

        coarseMaterialParams_.setkrRegardsSnr(true);
        fineMaterialParams_.setkrRegardsSnr(true);

        // parameters for adsorption
        coarseMaterialParams_.setKdNAPL(0.);
        coarseMaterialParams_.setRhoBulk(1500.);
        fineMaterialParams_.setKdNAPL(0.);
        fineMaterialParams_.setRhoBulk(1500.);
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
    { return "kuevette3p3cni"; }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    Scalar intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
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
     * \brief return the parameter object for the material law which depends on the position
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
     * \brief Return the parameter object for the heat conductivty law
     *        for a given position
     */
    template <class Context>
    const HeatConductionLawParams&
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        return heatCondParams_;
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
            850 // specific heat capacity [J / (kg K)]
            * 2650; // density of sand [kg/m^3]
    }

    template <class Context>
    void source(RateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    { values = 0; }

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
    void boundaryTypes(BoundaryTypes &values, 
                       const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition globalPos = context.pos(spaceIdx, timeIdx);

        if(globalPos[0] > 1.5 - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

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
    { initial(values, context, spaceIdx, timeIdx); }

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
        if (globalPos[0] < eps_)
        {
            values[Indices::contiWEqIdx] = -0.1435; // 0.3435 [mol/(s m)] in total
            values[Indices::contiAEqIdx] = -0.2;
            values[Indices::contiCEqIdx] =  0.0;
            values[Indices::energyEqIdx] = -6929.;
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
        values[pressure0Idx] = 1e5 ;
        values[switch1Idx] = 0.12;
        values[switch2Idx] = 1.e-6;
        values[temperatureIdx] = 293.0;

        if((globalPos[0] >= 0.20) && (globalPos[0] <= 0.80) && (globalPos[1] >= 0.4) && (globalPos[1] <= 0.65))
        {
            values[switch2Idx] = 0.07;
            values.setPhasePresence(threePhases);
        }
        else
            values.setPhasePresence(wgPhaseOnly);
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
    {
        if (0.13 <= pos[0] && 1.20 >= pos[0] && 0.32 <= pos[1] && pos[1] <= 0.57)
            return true;
        else if (0.15 >= pos[1] && 1.20 <= pos[0])
            return true;
        else return false;
    };

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams heatCondParams_;

    const Scalar eps_;
};
} //end namespace

#endif
