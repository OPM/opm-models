// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
/**
 * \file
 * \brief Problem where water and gas is injected by means of a
 *        dirchlet condition on the lower right of the domain which have to go
 *        around an obstacle with \f$10^3\f$ lower permeability.
 * \author Andreas Lauser, Klaus Mosthaf, Bernd Flemisch
 */
#ifndef DUMUX_OBSTACLEPROBLEM_HH
#define DUMUX_OBSTACLEPROBLEM_HH

#if OBSTACLE_USE_2P2C
#include <dumux/boxmodels/2p2c/2p2cmodel.hh>
#else
#include <dumux/boxmodels/mpnc/mpncmodel.hh>
#endif
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedlinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/mp/mplinearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>
#include <dumux/material/heatconduction/somerton.hh>

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

#include <iostream>

namespace Dumux
{

template <class TypeTag>
class ObstacleProblem;

namespace Properties
{
#if OBSTACLE_USE_2P2C
NEW_TYPE_TAG(ObstacleProblem, INHERITS_FROM(BoxTwoPTwoC));
#else
NEW_TYPE_TAG(ObstacleProblem, INHERITS_FROM(BoxMPNC));
#endif

// Set the grid type
SET_TYPE_PROP(ObstacleProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ObstacleProblem,
              Problem,
              Dumux::ObstacleProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(ObstacleProblem,
              FluidSystem,
              Dumux::FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the material Law
SET_PROP(ObstacleProblem, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum {
        lPhaseIdx = FluidSystem::lPhaseIdx,
        gPhaseIdx = FluidSystem::gPhaseIdx
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
SET_PROP(ObstacleProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Dumux::Somerton<FluidSystem, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(ObstacleProblem, EnableGravity, true);

// Write Newton convergence to disk?
SET_BOOL_PROP(ObstacleProblem, NewtonWriteConvergence, false);

// Use the line search strategy for the Newton update?
SET_BOOL_PROP(ObstacleProblem, NewtonUseLineSearch, false);

// Enable the re-use of the jacobian matrix whenever possible?
SET_BOOL_PROP(ObstacleProblem, EnableJacobianRecycling, false);

// Reassemble the jacobian matrix only where it changed?
SET_BOOL_PROP(ObstacleProblem, EnablePartialReassemble, false);

// use forward diffferences to approximate the partial derivatives
SET_INT_PROP(ObstacleProblem, NumericDifferenceMethod, +1);

#if OBSTACLE_USE_2P2C
/////////
// 2p-2c specific properties
/////////

// Verbosity of the 2p2c model (0=silent, 1=medium, 2=chatty)
SET_INT_PROP(ObstacleProblem, TwoPTwoCVerbosity, 1);

#else
/////////
// Mp-Nc specific properties
/////////

// Enable molecular diffusion of the components?
SET_BOOL_PROP(ObstacleProblem, EnableDiffusion, false);

// number of Newton iterations per time step where the update gets limited
SET_INT_PROP(ObstacleProblem, NewtonChoppedIterations, 0);
#endif
}


/*!
 * \ingroup MPNCModel
 * \ingroup BoxTestProblems
 * \brief Problem where liquid water is injected by means of a
 *        dirchlet condition on the lower right of the domain which have to go
 *        around an obstacle with \f$10^3\f$ lower permeability.
 *
 * The domain is sized 60m times 40m and consists of two media, a
 * moderately permeable soil (\f$ K_0=10e-12 m^2\f$) and an obstacle
 * at \f$[10; 20]m \times [0; 35]m \f$ with a lower permeablility of
 * \f$ K_1=K_0/1000\f$.
 *
 * Initially the whole domain is filled by nitrogen, the temperature
 * is \f$20\symbol{23}C\f$ at the whole domain. The gas pressure at
 * the left side of the domain is 1 bar, at the right side it is 2 bar
 * with a linear gradient in between.
 *
 * The boundary is no-flow except on the lower 10 meters of the left
 * and the right boundary which is a Dirichlet condition with the same
 * values as the initial condition.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_MpNc -parameterFile test_MpNc.input</tt>
 */
template <class TypeTag>
class ObstacleProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename HeatConductionLaw::Params HeatConductionLawParams;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),

        gPhaseIdx = FluidSystem::gPhaseIdx,
        lPhaseIdx = FluidSystem::lPhaseIdx,

        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx
    };

    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

public:
    ObstacleProblem(TimeManager &timeManager)
        : ParentType(timeManager, 
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 1e-6;
        temperature_ = 273.15 + 25; // -> 25°C

        // initialize the tables of the fluid system
        Scalar Tmin = temperature_ - 10.0;
        Scalar Tmax = temperature_ + 40.0;
        int nT = 3;

        Scalar pmin = 1.0e5 * 0.75;
        Scalar pmax = 2.0e5 * 1.25;
        int np = 1000;

        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

        // intrinsic permeabilities
        coarseK_ = this->toDimMatrix_(1e-12);
        fineK_ = this->toDimMatrix_(1e-15);

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

        initFluidStates_();
    }

    /*!
     * \brief Called directly after successful time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms of the individual phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            PrimaryVariables phaseStorage;
            this->model().globalPhaseStorage(phaseStorage, phaseIdx);

            if (this->gridView().comm().rank() == 0) {
                std::cout
                    <<"Storage in "
                    << FluidSystem::phaseName(phaseIdx)
                    << "Phase: ["
                    << phaseStorage
                    << "]"
                    << "\n";
            }
        }

        // Calculate total storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout
                <<"Storage total: [" << storage << "]"
                << "\n";
        }
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
    { return "obstacle"; }

    /*!
     * \brief Returns the temperature [K] within the domain.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    /*!
     * \brief Returns the intrinsic permeability tensor.
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        if (isFineMaterial_(context.pos(spaceIdx, timeIdx)))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
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
     * \brief Function for defining the parameters needed by constitutive relationships (kr-Sw, pc-Sw, etc.).
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
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

        if (onInlet_(pos))
            values.setFreeFlow(context, spaceIdx, timeIdx, inletFluidState_);
        else if (onOutlet_(pos))
            values.setFreeFlow(context, spaceIdx, timeIdx, outletFluidState_);
        else
            values.setNoFlow();
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of all conservation quantities that is generated or annihilated
     * per volume unit. Positive values mean that mass is created,
     * negative ones mean that it vanishes.
     */
    template <class Context>
    void source(RateVector &rate,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { rate = 0.0; }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &priVars, const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        priVars.assignMassConservative(outletFluidState_,
                                       matParams,
                                       /*inThermodynamicEquilibrium=*/false);
    }

    // \}

protected:
    /*!
     * \brief Returns whether a given global position is in the
     *        fine-permeability region or not.
     */
    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        return
            10 <= pos[0] && pos[0] <= 20 &&
            0 <= pos[1] && pos[1] <= 35;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x >= 60 - eps_ && y <= 10;
    }

    bool onOutlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x < eps_ && y <= 10;
    }

    
    void initFluidStates_()
    {
        initFluidState_(inletFluidState_, coarseMaterialParams_, /*isInlet=*/true);
        initFluidState_(outletFluidState_, coarseMaterialParams_, /*isInlet=*/false);
    }
    
    template <class FluidState>
    void initFluidState_(FluidState &fs, const MaterialLawParams &matParams, bool isInlet)
    {
        int refPhaseIdx;
        int otherPhaseIdx;

        // set the fluid temperatures
        fs.setTemperature(temperature_);

        if (isInlet) {
            // only liquid on inlet
            refPhaseIdx = lPhaseIdx;
            otherPhaseIdx = gPhaseIdx;

            // set liquid saturation
            fs.setSaturation(lPhaseIdx, 1.0);

            // set pressure of the liquid phase
            fs.setPressure(lPhaseIdx, 2e5);

            // set the liquid composition to pure water
            fs.setMoleFraction(lPhaseIdx, N2Idx, 0.0);
            fs.setMoleFraction(lPhaseIdx, H2OIdx, 1.0);
        }
        else {
            // elsewhere, only gas
            refPhaseIdx = gPhaseIdx;
            otherPhaseIdx = lPhaseIdx;

            // set gas saturation
            fs.setSaturation(gPhaseIdx, 1.0);

            // set pressure of the gas phase
            fs.setPressure(gPhaseIdx, 1e5);

            // set the gas composition to 99% nitrogen and 1% steam
            fs.setMoleFraction(gPhaseIdx, N2Idx, 0.99);
            fs.setMoleFraction(gPhaseIdx, H2OIdx, 0.01);
        };

        // set the other saturation
        fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

        // calulate the capillary pressure
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);
        fs.setPressure(otherPhaseIdx,
                       fs.pressure(refPhaseIdx)
                       + (pC[otherPhaseIdx] - pC[refPhaseIdx]));

        // make the fluid state consistent with local thermodynamic
        // equilibrium
        typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

        typename FluidSystem::ParameterCache paramCache;
        ComputeFromReferencePhase::solve(fs,
                                         paramCache,
                                         refPhaseIdx,
                                         /*setViscosity=*/false,
                                         /*setEnthalpy=*/false);
    }

    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1-poro)) * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1-poro));

        params.setFullySaturatedLambda(gPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(lPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    DimMatrix coarseK_;
    DimMatrix fineK_;

    Scalar coarsePorosity_;
    Scalar finePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams fineHeatCondParams_;
    HeatConductionLawParams coarseHeatCondParams_;

    CompositionalFluidState<Scalar, FluidSystem> inletFluidState_;
    CompositionalFluidState<Scalar, FluidSystem> outletFluidState_;
    
    Scalar temperature_;
    Scalar eps_;
};
} //end namespace

#endif
