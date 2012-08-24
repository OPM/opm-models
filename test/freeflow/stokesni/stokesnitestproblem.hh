// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * @file
 * @brief  Definition of a simple Stokes problem
 * @author Klaus Mosthaf, Andreas Lauser, Bernd Flemisch
 */
#ifndef DUMUX_STOKES_NI_TEST_PROBLEM_HH
#define DUMUX_STOKES_NI_TEST_PROBLEM_HH

#include <dumux/freeflow/stokes/stokesmodel.hh>
#include <dumux/common/simplexgridcreator.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{

template <class TypeTag>
class StokesNITestProblem;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesNITestProblem, INHERITS_FROM(BoxStokes));

// Set the grid type
SET_TYPE_PROP(StokesNITestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesNITestProblem, Problem, StokesNITestProblem<TypeTag>);

//! Select the fluid system
SET_TYPE_PROP(StokesNITestProblem, 
              FluidSystem,
              Dumux::FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Select the phase to be considered
SET_INT_PROP(StokesNITestProblem,
             StokesPhaseIndex,
             GET_PROP_TYPE(TypeTag, FluidSystem)::gPhaseIdx);

// Enable gravity
SET_BOOL_PROP(StokesNITestProblem, EnableGravity, true);

// Enable the energy equation
SET_BOOL_PROP(StokesNITestProblem, EnableEnergy, true);

// Enable constraints
SET_BOOL_PROP(StokesNITestProblem, EnableConstraints, true);
}

/*!
 * \ingroup BoxStokesNIModel
 * \ingroup BoxTestProblems
 * \brief StokesNI problem with air (N2) flowing
 *        from the left to the right.
 *
 * The domain of this problem is 1m times 1m. The upper and the lower
 * boundaries are fixed to the initial condition by means of
 * constraints, the left and the right boundaries are no-slip
 * conditions.
 */
template <class TypeTag>
class StokesNITestProblem 
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { // Number of equations and grid dimension
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dimWorld = GridView::dimensionworld
    };
    enum {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        moleFrac1Idx = Indices::moleFrac1Idx,
        velocity0Idx = Indices::velocity0Idx,
        temperatureIdx = Indices::temperatureIdx,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        momentum0EqIdx = Indices::momentum0EqIdx,
        energyEqIdx = Indices::energyEqIdx
    };
    enum { numComponents = FluidSystem::numComponents };
    enum { H2OIdx = FluidSystem::H2OIdx };
    enum { AirIdx = FluidSystem::AirIdx };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    StokesNITestProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 1e-6;

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/280.0, /*Tmax=*/285, /*nT=*/10,
                          /*pmin=*/1e5, /*pmax=*/1e5 + 100, /*np=*/200);
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
    { return "stokesni"; }

    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 273.15 + 10; }

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

        if (onUpperBoundary_(pos))
            values.setOutFlow(context, spaceIdx, timeIdx);
        else if(onLowerBoundary_(pos)) {
            // lower boundary is constraint!
            values = 0.0;
        }
        else {
            // left and right
            values.setNoFlow(context, spaceIdx, timeIdx);
        }
    }

    // \}

    /*!
     * \name Constraints
     */
    // \{

    /*!
     * \brief Set the constraints of this problem.
     *
     * This method sets temperature constraints for the finite volumes
     * adacent to the inlet.
     */
    template <class Context>
    void constraints(Constraints &values,
                     const Context &context,
                     int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onLowerBoundary_(pos) || onUpperBoundary_(pos))
        {
            PrimaryVariables initCond;
            initial(initCond, context, spaceIdx, timeIdx);

            values.setConstraint(temperatureIdx, energyEqIdx, initCond[temperatureIdx]);;
            values.setConstraint(pressureIdx, conti0EqIdx, initCond[pressureIdx]);
            values.setConstraint(moleFrac1Idx, conti0EqIdx+1, initCond[moleFrac1Idx]);;
            for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
                values.setConstraint(velocity0Idx + axisIdx,
                                     momentum0EqIdx + axisIdx,
                                     initCond[momentum0EqIdx + axisIdx]);
        }
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
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template <class Context>
    void source(RateVector &values,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { values = Scalar(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        
        Scalar moleFrac[numComponents];

        moleFrac[H2OIdx] = 1e-4;
        Scalar temperature = 283.15;
        if (inLens_(pos)) {
            moleFrac[H2OIdx] = 0.9e-4;
            temperature = 284.15;
        };
        moleFrac[AirIdx] = 1 - moleFrac[H2OIdx];
        
        
        // parabolic velocity profile
        Scalar y = this->bboxMax()[1] - pos[1];
        Scalar x = pos[0] - this->bboxMin()[0];
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];

        // parabolic velocity profile
        const Scalar maxVelocity = 1.0;

        Scalar a = - 4*maxVelocity/(width*width);
        Scalar b = - a*width;
        Scalar c = 0;
        
        DimVector velocity(0.0);
        velocity[1] = a * x*x + b * x + c;

        // hydrostatic pressure
        Scalar rho = 1.189;
        Scalar pressure = 1e5 - rho*this->gravity()[1]*y;

        for (int axisIdx = 0; axisIdx < dimWorld; ++ axisIdx)
            values[velocity0Idx + axisIdx] = velocity[axisIdx];

        values[pressureIdx] = pressure;
        values[moleFrac1Idx] = moleFrac[1];
        values[temperatureIdx] = temperature;
    }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[1] > this->bboxMax()[1] - eps_; }

    bool onBoundary_(const GlobalPosition &pos) const
    {
        return onLeftBoundary_(pos)
            || onRightBoundary_(pos)
            || onLowerBoundary_(pos)
            || onUpperBoundary_(pos);
    }

    bool inLens_(const GlobalPosition &pos) const
    {
        return 
            pos[0]<0.75 && pos[0]>0.25 &&
            pos[1]<0.75 && pos[1]>0.25;
    }

    Scalar eps_;
};
} //end namespace

#endif
