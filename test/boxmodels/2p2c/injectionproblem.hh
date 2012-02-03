// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
/*!
 * \file
 *
 * \brief Definition of a problem, where air is injected under a low permeable layer.
 */
#ifndef DUMUX_INJECTION_PROBLEM_HH
#define DUMUX_INJECTION_PROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfs.hh>

#include <dumux/boxmodels/2p2c/2p2cmodel.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>

#include "injectionspatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem, INHERITS_FROM(BoxTwoPTwoC, InjectionSpatialParameters));

// declare some injection problem specific property tags
NEW_PROP_TAG(FluidSystemPressureLow);
NEW_PROP_TAG(FluidSystemPressureHigh);
NEW_PROP_TAG(FluidSystemNumPressure);
NEW_PROP_TAG(FluidSystemTemperatureLow);
NEW_PROP_TAG(FluidSystemTemperatureHigh);
NEW_PROP_TAG(FluidSystemNumTemperature);

NEW_PROP_TAG(InitialConditionsMaxDepth);
NEW_PROP_TAG(InitialConditionsTemperature);
NEW_PROP_TAG(SimulationControlName);

// Set the grid type
SET_PROP(InjectionProblem, Grid)
{
    typedef Dune::SGrid<2,2> type;
};

// Set the problem property
SET_PROP(InjectionProblem, Problem)
{
    typedef Dumux::InjectionProblem<TypeTag> type;
};

// Set fluid configuration
SET_PROP(InjectionProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    static const bool useComplexRelations = false;
public:
    typedef Dumux::FluidSystems::H2ON2<Scalar, useComplexRelations> type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem, EnableGravity, true);

SET_BOOL_PROP(InjectionProblem, EnableJacobianRecycling, true);
SET_BOOL_PROP(InjectionProblem, EnableVelocityOutput, false);

// set the defaults for some injection problem specific properties
SET_SCALAR_PROP(InjectionProblem, FluidSystemPressureLow, 1e6);
SET_SCALAR_PROP(InjectionProblem, FluidSystemPressureHigh, 3e7);
SET_INT_PROP(InjectionProblem, FluidSystemNumPressure, 100);
SET_SCALAR_PROP(InjectionProblem, FluidSystemTemperatureLow, 273.15);
SET_SCALAR_PROP(InjectionProblem, FluidSystemTemperatureHigh, 373.15);
SET_INT_PROP(InjectionProblem, FluidSystemNumTemperature, 100);

SET_SCALAR_PROP(InjectionProblem, InitialConditionsMaxDepth, 2500);
SET_SCALAR_PROP(InjectionProblem, InitialConditionsTemperature, 293.15);
SET_STRING_PROP(InjectionProblem, SimulationControlName, "injection");
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup BoxTestProblems
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable spatial parameters (\f$ K=10e-12\f$) for \f$ y>22m\f$ and one with a lower permeablility (\f$ K=10e-13\f$)
 * in the rest of the domain.
 *
 * Air enters a water-filled aquifer, which is situated 2700m below sea level, at the right boundary
 * (\f$ 5m<y<15m\f$) and migrates upwards due to buoyancy. It accumulates and
 * partially enters the lower permeable aquitard.
 * This problem uses the \ref TwoPTwoCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2p2c -parameterFile ./test_2p2c.input</tt>
 */
template <class TypeTag>
class InjectionProblem : public TwoPTwoCProblem<TypeTag>
{
    typedef TwoPTwoCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, TwoPTwoCIndices) Indices;
    enum {
        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,


        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx,

        conti0EqIdx = Indices::conti0EqIdx,
        contiN2EqIdx = conti0EqIdx + N2Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {      
        eps_ = 1e-6;

        temperatureLow_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureLow);
        temperatureHigh_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureHigh);
        nTemperature_ = GET_PARAM_FROM_GROUP(TypeTag, int, FluidSystem, NumTemperature);

        nPressure_ = GET_PARAM_FROM_GROUP(TypeTag, int, FluidSystem, NumPressure);
        pressureLow_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureLow);
        pressureHigh_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureHigh);
        
        temperature_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, InitialConditions, Temperature);
        maxDepth_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, InitialConditions, MaxDepth);
        name_ = GET_PARAM_FROM_GROUP(TypeTag, std::string, SimulationControl, Name);

        // initialize the tables of the fluid system
        //FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
    }

    /*!
     * \brief Called directly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, lPhaseIdx);
        this->model().globalPhaseStorage(storageG, gPhaseIdx);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: liquid=[" << storageL << "]"
                     << " gas=[" << storageG << "]\n";
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
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; };

    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
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
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        if (globalPos[0] < eps_)
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
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        values = 0;
        if (globalPos[1] < 15 && globalPos[1] > 5) {
            values[contiN2EqIdx] = -1e-3; // kg/(s*m^2)
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
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        initial_(values, globalPos);

    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vert The vertex
     * \param globalIdx The index of the global vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vert,
                             int &globalIdx,
                             const GlobalPosition &globalPos) const
    { return Indices::lPhaseOnly; }

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = FluidSystem::H2O::liquidDensity(temperature_, 1e5);

        Scalar pl = 1e5 - densityW*this->gravity()[1]*(maxDepth_ - globalPos[1]);
        Scalar moleFracLiquidN2 = pl*0.95/BinaryCoeff::H2O_N2::henry(temperature_);
        Scalar moleFracLiquidH2O = 1.0 - moleFracLiquidN2;

        Scalar meanM =
            FluidSystem::molarMass(H2OIdx)*moleFracLiquidH2O +
            FluidSystem::molarMass(N2Idx)*moleFracLiquidN2;

        Scalar massFracLiquidN2 = moleFracLiquidN2*FluidSystem::molarMass(N2Idx)/meanM;

        values[Indices::pressureIdx] = pl;
        values[Indices::switchIdx] = massFracLiquidN2;
    }

    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;

    int nTemperature_;
    int nPressure_;

    std::string name_ ;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
};
} //end namespace

#endif
