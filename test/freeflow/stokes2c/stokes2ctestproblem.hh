// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Klaus Mosthaf                                     *
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
 * \brief A simple problem using the isothermal Stokes model and two components.
 */
#ifndef DUMUX_STOKES_2C_TEST_PROBLEM_HH
#define DUMUX_STOKES_2C_TEST_PROBLEM_HH

#include <dumux/freeflow/stokes/stokesmodel.hh>
#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cTestProblem;

//////////
// Specify the properties for the stokes2c problem
//////////
namespace Properties
{
NEW_TYPE_TAG(Stokes2cTestProblem, INHERITS_FROM(BoxStokes));

// Set the grid type
SET_TYPE_PROP(Stokes2cTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(Stokes2cTestProblem, Problem, Dumux::Stokes2cTestProblem<TypeTag>);

//! Select the fluid system
SET_TYPE_PROP(Stokes2cTestProblem, 
              FluidSystem,
              Dumux::FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Select the phase to be considered
SET_INT_PROP(Stokes2cTestProblem,
             StokesPhaseIndex,
             GET_PROP_TYPE(TypeTag, FluidSystem)::gPhaseIdx);

// Disable gravity
SET_BOOL_PROP(Stokes2cTestProblem, EnableGravity, false);

// Enable constraints
SET_BOOL_PROP(Stokes2cTestProblem, EnableConstraints, true);
}

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxTestProblems
 * \brief Stokes transport problem with air flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the momentum balances
 * are all set to Dirichlet. The mass balance receives
 * outflow bcs, which are replaced in the localresidual by the sum
 * of the two momentum balances. In the middle of the right boundary,
 * one vertex receives Dirichlet bcs, to set the pressure level.
 *
 * This problem uses the \ref BoxStokes2cModel.
 */
template <class TypeTag>
class Stokes2cTestProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { dimWorld = GridView::dimensionworld };
    enum { numComponents = FluidSystem::numComponents };
    enum {
        // copy some indices for convenience
        conti0EqIdx = Indices::conti0EqIdx,
        momentum0EqIdx = Indices::momentum0EqIdx,

        velocity0Idx = Indices::velocity0Idx,
        moleFrac1Idx = Indices::moleFrac1Idx,
        pressureIdx = Indices::pressureIdx,

        H2OIdx = FluidSystem::H2OIdx,
        AirIdx = FluidSystem::AirIdx
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    Stokes2cTestProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 1e-6;

        // initialize the tables of the fluid system
        FluidSystem::init();
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
    { return "stokes2c"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    template <class Context>   
    Scalar temperature(const Context &context,
                       int spaceIdx, int timeIdx) const
    { return 273.15 + 10; /* -> 10 deg C */ }

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

        if (onLowerBoundary_(pos))
            values.setOutFlow(context, spaceIdx, timeIdx);
        else if(onUpperBoundary_(pos)) {
            // upper boundary is constraint!
            values = 0.0;
        }
        else {
            // left and right boundaries
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

        if (onUpperBoundary_(pos)) {
            PrimaryVariables initCond;
            initial(initCond, context, spaceIdx, timeIdx);

            values.setConstraint(pressureIdx, conti0EqIdx, initCond[pressureIdx]);;
            values.setConstraint(moleFrac1Idx, conti0EqIdx + 1, initCond[moleFrac1Idx]);
            for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
                values.setConstraint(velocity0Idx + axisIdx,
                                     momentum0EqIdx + axisIdx,
                                     initCond[velocity0Idx + axisIdx]);
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
    {
        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0.0);
    }

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
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        values = 0.0;

        //parabolic profile
        const Scalar v1 = 1.0;
        values[velocity0Idx + 1] =
            - v1*(globalPos[0] - this->bboxMin()[0])*(this->bboxMax()[0] - globalPos[0])
            / (0.25*(this->bboxMax()[0] - this->bboxMin()[0])*(this->bboxMax()[0] - this->bboxMin()[0]));

        Scalar moleFrac[numComponents];
        if (onUpperBoundary_(globalPos))
            moleFrac[H2OIdx] = 0.005;
        else
            moleFrac[H2OIdx] = 0.007;
        moleFrac[AirIdx] = 1.0 - moleFrac[H2OIdx];

        values[pressureIdx] = 1e5;
        values[velocity0Idx + 0] = 0.0;
        values[moleFrac1Idx] = moleFrac[1];
    }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bboxMax()[1] - eps_; }

    Scalar eps_;
};
} //end namespace

#endif
