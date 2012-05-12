// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
#ifndef DUMUX_CUVETTE_PROBLEM_HH
#define DUMUX_CUVETTE_PROBLEM_HH

#include <dumux/boxmodels/3p3cni/3p3cnimodel.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>
#include <dumux/material/fluidmatrixinteractions/3p/3pparkervangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/mp/3padapter.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mplinearmaterial.hh>
#include <dumux/material/heatconduction/somerton.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{
template <class TypeTag>
class CuvetteProblem;

namespace Properties
{
// create a new type tag for the cuvette steam injection problem
NEW_TYPE_TAG(CuvetteProblem, INHERITS_FROM(BoxThreePThreeCNI));

// Set the grid type
SET_TYPE_PROP(CuvetteProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(CuvetteProblem, Problem, Dumux::CuvetteProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(CuvetteProblem, 
              FluidSystem,
              Dumux::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(CuvetteProblem, EnableGravity, true);

// Set the maximum time step
SET_SCALAR_PROP(CuvetteProblem, MaxTimeStepSize, 600.);

// Set the material Law
SET_PROP(CuvetteProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;


    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { gPhaseIdx = FluidSystem::gPhaseIdx };

    // define the three-phase material law
    typedef Dumux::ThreePParkerVanGenuchten<Scalar> ThreePLaw;

    
public:
    // wrap the three-phase law in an adaptor to make use the generic
    // material law API
    typedef Dumux::ThreePAdapter<wPhaseIdx, nPhaseIdx, gPhaseIdx, ThreePLaw> type;

    //typedef Dumux::MpLinearMaterial<FluidSystem::numPhases, Scalar> type;
};

// Set the heat conduction law
SET_PROP(CuvetteProblem, HeatConductionLaw)
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
 * The domain is a quasi-two-dimensional container (cuvette). Its dimensions
 * are 1.5 m x 0.74 m. The top and bottom boundaries are closed, the right
 * boundary is a free-flow boundary allowing fluids to escape. From the left,
 * an injection of a hot water-air mixture is applied (Inflow boundary condition
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
 * To see the basic effect and the differences to scenarios with pure steam or
 * pure air injection, it is sufficient to simulated for about 2-3 hours (10000 s).
 * Complete remediation of the domain requires much longer (about 10 days simulated time).
 */
template <class TypeTag >
class CuvetteProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
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
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        numPhases = FluidSystem::numPhases,
        numComponents = FluidSystem::numComponents,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx,

        H2OIdx = FluidSystem::H2OIdx,
        airIdx = FluidSystem::airIdx,
        NAPLIdx = FluidSystem::NAPLIdx,

        conti0EqIdx = Indices::conti0EqIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    CuvetteProblem(TimeManager &timeManager)
        : ParentType(timeManager, 
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(1e-6)
    {
        FluidSystem::init(/*minT=*/283.15, /*maxT=*/500.0, /*nT=*/200,
                          /*minp=*/0.8e5, /*maxp=*/2e5, /*np=*/100);

        // intrinsic permeabilities
        fineK_ = this->toTensor_(6.28e-12);
        coarseK_ = this->toTensor_(9.14e-10);

        // porosities
        finePorosity_ = 0.42;
        coarsePorosity_ = 0.42;

        // parameters for the capillary pressure law
#if 1
        // three-phase van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.005);
        fineMaterialParams_.setVgN(4.0);
        coarseMaterialParams_.setVgN(4.0);

        coarseMaterialParams_.setkrRegardsSnr(true);
        fineMaterialParams_.setkrRegardsSnr(true);

        coarseMaterialParams_.setKdNAPL(0.);
        coarseMaterialParams_.setRhoBulk(1500.);
        fineMaterialParams_.setKdNAPL(0.);
        fineMaterialParams_.setRhoBulk(1500.);

        // residual saturations
        fineMaterialParams_.setSwr(0.1201);
        fineMaterialParams_.setSwrx(0.1201);
        fineMaterialParams_.setSnr(0.0701);
        fineMaterialParams_.setSgr(0.0101);
        coarseMaterialParams_.setSwr(0.1201);
        coarseMaterialParams_.setSwrx(0.1201);
        coarseMaterialParams_.setSnr(0.0701);
        coarseMaterialParams_.setSgr(0.0101);
#else
        // linear material law
        fineMaterialParams_.setPcMinSat(gPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(gPhaseIdx, 0);
        fineMaterialParams_.setPcMinSat(nPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(nPhaseIdx, -1000);
        fineMaterialParams_.setPcMinSat(wPhaseIdx, 0);
        fineMaterialParams_.setPcMaxSat(wPhaseIdx, -10000);

        coarseMaterialParams_.setPcMinSat(gPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(gPhaseIdx, 0);
        coarseMaterialParams_.setPcMinSat(nPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(nPhaseIdx, -100);
        coarseMaterialParams_.setPcMinSat(wPhaseIdx, 0);
        coarseMaterialParams_.setPcMaxSat(wPhaseIdx, -1000);

        // residual saturations
        fineMaterialParams_.setResidSat(wPhaseIdx, 0.1201);
        fineMaterialParams_.setResidSat(nPhaseIdx, 0.0701);
        fineMaterialParams_.setResidSat(gPhaseIdx, 0.0101);

        coarseMaterialParams_.setResidSat(wPhaseIdx, 0.1201);
        coarseMaterialParams_.setResidSat(nPhaseIdx, 0.0701);
        coarseMaterialParams_.setResidSat(gPhaseIdx, 0.0101);
#endif

        // initialize parameters for the heat conduction law
        computeHeatCondParams_(heatCondParams_, finePorosity_);

        initInjectFluidState_();
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
    { 
        static std::string tmp = std::string("cuvette")+this->model().name();
        return tmp.c_str();
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
    const Tensor &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 293.15; /* [K] */ }

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
    { return heatCondParams_; }

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
     * \brief Evaluate the boundary conditions.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onRightBoundary_(pos)) {
            CompositionalFluidState<Scalar, FluidSystem> fs;

            initialFluidState_(fs, context, spaceIdx, timeIdx);
            
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
            values.setNoFlow();
        }
        else  if (onLeftBoundary_(pos))
        {
            // injection
            RateVector molarRate;
            
            // inject with the same composition as the gas phase of
            // the injection fluid state
            Scalar molarInjectionRate = 0.3435; // [mol/(m^2 s)]
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                molarRate[conti0EqIdx + compIdx] =
                    - molarInjectionRate
                    * injectFluidState_.moleFraction(gPhaseIdx, compIdx);

            // calculate the total mass injection rate [kg / (m^2 s)
            Scalar massInjectionRate =
                molarInjectionRate
                * injectFluidState_.averageMolarMass(gPhaseIdx);

            // set the boundary rate vector
            values.setMolarRate(molarRate);
            values.setEnthalpyRate(- injectFluidState_.enthalpy(gPhaseIdx)
                                   * massInjectionRate); // [J / (m^2 s)]
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
        CompositionalFluidState<Scalar, FluidSystem> fs;

        initialFluidState_(fs, context, spaceIdx, timeIdx);

        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(fs, matParams, /*inEquilibrium=*/false);
    }

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->bboxMax()[1] - eps_; }

    bool isContaminated_(const GlobalPosition &pos) const
    {
        return 
            (0.20 <= pos[0]) && (pos[0] <= 0.80)
            && (0.4 <= pos[1]) && (pos[1] <= 0.65);
    };

    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        if (0.13 <= pos[0] && 1.20 >= pos[0] && 0.32 <= pos[1] && pos[1] <= 0.57)
            return true;
        else if (pos[1] <= 0.15 && 1.20 <= pos[0])
            return true;
        else return false;
    };

    template <class FluidState, class Context>
    void initialFluidState_(FluidState &fs,
                            const Context &context,
                            int spaceIdx,
                            int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);

        fs.setTemperature(293.0 /*[K]*/);

        Scalar pw = 1e5;

        if(isContaminated_(pos)) {
            fs.setSaturation(wPhaseIdx, 0.12);
            fs.setSaturation(nPhaseIdx, 0.07);
            fs.setSaturation(gPhaseIdx, 1 - 0.12 - 0.07);

            // set the capillary pressures
            const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
            Scalar pc[numPhases];
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, 
                               pw + (pc[phaseIdx] - pc[wPhaseIdx]));
            
            // compute the phase compositions
            typedef MiscibleMultiPhaseComposition<Scalar, FluidSystem> MMPC;
            typename FluidSystem::ParameterCache paramCache;
            MMPC::solve(fs, paramCache, /*setViscosity=*/false, /*setEnthalpy=*/false);
        }
        else {         
            fs.setSaturation(wPhaseIdx, 0.12);
            fs.setSaturation(gPhaseIdx, 1 - fs.saturation(wPhaseIdx));
            fs.setSaturation(nPhaseIdx, 0);

            // set the capillary pressures
            const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
            Scalar pc[numPhases];
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fs.setPressure(phaseIdx, 
                               pw + (pc[phaseIdx] - pc[wPhaseIdx]));

            // compute the phase compositions
            typedef MiscibleMultiPhaseComposition<Scalar, FluidSystem> MMPC;
            typename FluidSystem::ParameterCache paramCache;
            MMPC::solve(fs, paramCache, /*setViscosity=*/false, /*setEnthalpy=*/false);
            
            // set the contaminant mole fractions to zero. this is a
            // little bit hacky...
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                fs.setMoleFraction(phaseIdx, NAPLIdx, 0.0);

                if (phaseIdx == nPhaseIdx)
                    continue;

                Scalar sumx = 0;
                for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                    sumx += fs.moleFraction(phaseIdx, compIdx);

                for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                    fs.setMoleFraction(phaseIdx, 
                                       compIdx, 
                                       fs.moleFraction(phaseIdx, compIdx) / sumx);
            }
        }
    }

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
                lambdaSaturated = std::pow(lambdaGranite, (1-poro)) + std::pow(lambdaFluid, poro);
            }
            else
                lambdaSaturated = std::pow(lambdaGranite, (1-poro));
            
            params.setFullySaturatedLambda(phaseIdx, lambdaSaturated);
            if (!FluidSystem::isLiquid(phaseIdx))
                params.setVacuumLambda(lambdaSaturated);
        }
    }

    void initInjectFluidState_()
    {
        injectFluidState_.setTemperature(383.0); // [K]
        injectFluidState_.setPressure(gPhaseIdx, 1e5); // [Pa]
        injectFluidState_.setSaturation(gPhaseIdx, 1.0); // [-]

        Scalar xgH2O = 0.417;
        injectFluidState_.setMoleFraction(gPhaseIdx, 
                                          H2OIdx, 
                                          xgH2O); // [-]
        injectFluidState_.setMoleFraction(gPhaseIdx, 
                                          airIdx, 
                                          1 - xgH2O); // [-]
        injectFluidState_.setMoleFraction(gPhaseIdx, 
                                          NAPLIdx,
                                          0.0); // [-]
       
        // set the specific enthalpy of the gas phase
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(injectFluidState_, gPhaseIdx);

        Scalar h = FluidSystem::enthalpy(injectFluidState_, paramCache, gPhaseIdx);
        injectFluidState_.setEnthalpy(gPhaseIdx, h);
    }

    Tensor fineK_;
    Tensor coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams heatCondParams_;

    CompositionalFluidState<Scalar, FluidSystem> injectFluidState_;

    const Scalar eps_;
};
} //end namespace

#endif
