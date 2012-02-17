// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                    *
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
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
#ifndef DUMUX_INFILTRATIONPROBLEM_HH
#define DUMUX_INFILTRATIONPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>

#include <dumux/boxmodels/3p3c/3p3cmodel.hh>

#include <dumux/material/fluidmatrixinteractions/3p/parkerVanGen3p.hh>
#include <dumux/material/fluidmatrixinteractions/mp/3padapter.hh>

#include <dumux/material/heatconduction/somerton.hh>

namespace Dumux
{
template <class TypeTag>
class InfiltrationProblem;

namespace Properties
{
NEW_TYPE_TAG(InfiltrationProblem, INHERITS_FROM(BoxThreePThreeC));

// Set the grid type
SET_TYPE_PROP(InfiltrationProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InfiltrationProblem, Problem, Dumux::InfiltrationProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(InfiltrationProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity?
SET_BOOL_PROP(InfiltrationProblem, EnableGravity, true);

// Write newton convergence?
SET_BOOL_PROP(InfiltrationProblem, NewtonWriteConvergence, true);

// Maximum tolerated relative error in the Newton method
SET_SCALAR_PROP(InfiltrationProblem, NewtonRelTolerance, 1e-8);

// -1 backward differences, 0: central differences, +1: forward differences
SET_INT_PROP(InfiltrationProblem, NumericDifferenceMethod, +1);

// Set the material Law
SET_PROP(InfiltrationProblem, MaterialLaw)
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
SET_PROP(InfiltrationProblem, HeatConductionLaw)
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
 * \ingroup ThreePThreeCBoxModel
 * \ingroup BoxTestProblems
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep, where
 * the lower part represents a slightly inclined groundwater table, and the
 * upper part is the vadose zone. 
 * A LNAPL (Non-Aqueous Phase Liquid which is lighter than water) infiltrates
 * (modelled with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and migrates
 * on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually trapped
 * immobile NAPL, which can in the following dissolve and evaporate slowly,
 * and eventually be transported by advection and diffusion.
 *
 * Left and right boundaries are constant head boundaries (Dirichlet),
 * Top and bottom are Neumann boundaries, all no-flow except for the small
 * infiltration zone in the upper left part.
 *
 * This problem uses the \ref ThreePThreeCModel.
 *
 * This problem should typically be simulated for 30 days.
 * A good choice for the initial time step size is 60 s.
 * To adjust the simulation time it is necessary to edit the file test_3p3cni.input
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_3p3c -parameterFile test_3p3c.input</tt>
 *  */
template <class TypeTag >
class InfiltrationProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef ThreePThreeCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, ThreePThreeCIndices) Indices;
    enum {
        pressure0Idx = Indices::pressure0Idx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        // Phase State
        wgPhaseOnly = Indices::wgPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InfiltrationProblem(TimeManager &timeManager)
        : ParentType(timeManager, 
                     GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(1e-6)
    {
        temperature_ = 273.15 + 10.0; // -> 10 degrees Celsius
        FluidSystem::init(/*tempMin=*/temperature_ - 1,
                          /*tempMax=*/temperature_ + 1,
                          /*nTemp=*/3,
                          /*pressMin=*/0.8*1e5,
                          /*pressMax=*/3*1e5,
                          /*nPress=*/200);

        // intrinsic permeabilities
        fineK_ = 1.e-11;
        coarseK_ = 1.e-11;

        // porosities
        porosity_ = 0.40;

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSwrx(0.12);
        materialParams_.setSnr(0.07);
        materialParams_.setSgr(0.03);

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(0.0005);
        materialParams_.setVgN(4.);
        materialParams_.setkrRegardsSnr(false);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);

        materialParams_.checkDefined();
    }

    /*!
     * \name Problem parameters
     */
    // \{
    
    bool shouldWriteRestartFile() const
    { return true; }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "infiltration"; }

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
        //const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        // if (isFineMaterial_(pos))
        //     return finePorosity_;
        // else
        //     return coarsePorosity_;
        return porosity_;
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
    { return materialParams_; }

#warning TODO
#if 0
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
#endif

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
            850. // specific heat capacity [J / (kg K)]
            * 2650.; // density of sand [kg/m^3]
    }

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
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

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
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        if(globalPos[0] > 500. - eps_)
            values.setAllDirichlet();
        else if(globalPos[0] < eps_)
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
    {
        const GlobalPosition globalPos = context.pos(spaceIdx, timeIdx);

        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar Sw;
        Scalar Swr=0.12;
        Scalar Sgr=0.03;

        values = 0.0;
        if (y >(-1.E-3*x+5)) {
            Scalar pc = 9.81 * 1000.0 * (y - (5 - 5e-4*x));
            if (pc < 0.0) pc = 0.0;

            Sw = invertPCGW_(pc, materialLawParams(context, spaceIdx, timeIdx));
            if (Sw < Swr) Sw = Swr;
            if (Sw > 1.-Sgr) Sw = 1.-Sgr;

            values[pressure0Idx] = 1e5 ;
            values[switch1Idx] = Sw;
            values[switch2Idx] = 1.e-6;

            values.setPhasePresent(FluidSystem::gPhaseIdx);
            values.setPhasePresent(FluidSystem::wPhaseIdx);

            Valgrind::CheckDefined(values);
        } 
        else {
            values[pressure0Idx] = 1e5 + 9.81 * 1000.0 * ((-5E-4*x+5) - y);
            values[switch1Idx] = 1.-Sgr;
            values[switch2Idx] = 1.e-6;

            values.setPhasePresent(FluidSystem::gPhaseIdx);
            values.setPhasePresent(FluidSystem::wPhaseIdx);

            Valgrind::CheckDefined(values);
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
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
        if ((globalPos[0] <= 75.+eps_) && (globalPos[0] >= 50.+eps_) && (globalPos[1] >= 10.-eps_))
        {
            values[Indices::contiWEqIdx] = -0.0;
            values[Indices::contiCEqIdx] = -0.001;
            values[Indices::contiAEqIdx] = -0.0;
            Valgrind::CheckDefined(values);
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

        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar Sw;
        Scalar Swr=0.12;
        Scalar Sgr=0.03;

        values = 0.0;
        if (y > -eps_*x + 5)
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5E-4*x+5));
            if (pc < 0.0) pc = 0.0;

            Sw = invertPCGW_(pc,
                             materialLawParams(context, spaceIdx, timeIdx));
            Valgrind::CheckDefined(Sw);
            if (Sw < Swr) Sw = Swr;
            if (Sw > 1.-Sgr) Sw = 1.-Sgr;

            values[pressure0Idx] = 1e5 ;
            values[switch1Idx] = Sw;
            values[switch2Idx] = 1.e-6;

            values.setPhasePresent(FluidSystem::gPhaseIdx);
            values.setPhasePresent(FluidSystem::wPhaseIdx);
            Valgrind::CheckDefined(values);
        } else {
            values[pressure0Idx] = 1e5 + 9.81 * 1000.0 * ((5 - x*5e-4) - y);
            Valgrind::CheckDefined(values[pressure0Idx]);
            Valgrind::CheckDefined(Sgr);
            values[switch1Idx] = 1.0 - Sgr;
            values[switch2Idx] = 1.0e-6;

            values.setPhasePresent(FluidSystem::gPhaseIdx);
            values.setPhasePresent(FluidSystem::wPhaseIdx);
            Valgrind::CheckDefined(values);
        }
    }

private:
    static Scalar invertPCGW_(Scalar pcIn, const MaterialLawParams &pcParams)
    {
        Scalar lower,upper;
        int k;
        int maxIt = 50;
        Scalar bisLimit = 1.;
        Scalar Sw, pcGW;
        lower=0.0; upper=1.0;
        for (k=1; k<=25; k++)
        {
            Sw = 0.5*(upper+lower);
            pcGW = MaterialLaw::pCGW(pcParams, Sw);
            Scalar delta = pcGW-pcIn;
            if (delta<0.) delta*=-1.;
            if (delta<bisLimit)
            {
                return(Sw);
            }
            if (k==maxIt) {
                return(Sw);
            }
            if (pcGW>pcIn) lower=Sw;
            else upper=Sw;
        }
        return(Sw);
    }

    bool isFineMaterial_(const GlobalPosition &pos) const
    {
        return
            70. <= pos[0] && pos[0] <= 85. &&
            7.0 <= pos[1] && pos[1] <= 7.50;
    };

    Scalar fineK_;
    Scalar coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;

    Scalar temperature_;
    Scalar eps_;
};
} //end namespace

#endif
