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
 * \brief  Definition of a simple Stokes problem
 */
#ifndef DUMUX_STOKESTESTPROBLEM_HH
#define DUMUX_STOKESTESTPROBLEM_HH

#include <dumux/freeflow/stokes/stokesmodel.hh>
#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{

template <class TypeTag>
class StokesTestProblem;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesTestProblem, INHERITS_FROM(BoxStokes));

// Set the grid type
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesTestProblem, Problem, Dumux::StokesTestProblem<TypeTag>);

SET_PROP(StokesTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::GasPhase<Scalar, Dumux::N2<Scalar> > type;
};

// Disable gravity
SET_BOOL_PROP(StokesTestProblem, EnableGravity, false);

// Enable constraints
SET_BOOL_PROP(StokesTestProblem, EnableConstraints, true);

// Default simulation end time [s]
SET_SCALAR_PROP(StokesTestProblem, EndTime, 6e3);

// Default initial time step size [s]
SET_SCALAR_PROP(StokesTestProblem, InitialTimeStepSize, 10.0);

// Default grid file to load
SET_STRING_PROP(StokesTestProblem, GridFile, "grids/test_stokes.dgf");
}

/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxTestProblems
 * \brief Stokes flow problem with nitrogen (N2) flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the momentum balances
 * are set to Dirichlet with outflow on the right boundary. The mass balance has
 * outflow bcs, which are replaced in the localresidual by the sum
 * of the two momentum balances. In the middle of the right boundary,
 * one vertex receives Dirichlet bcs to set the pressure level.
 *
 * This problem uses the \ref BoxStokesModel.
 */
template <class TypeTag>
class StokesTestProblem 
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;

    enum {
        // Number of equations and grid dimension
        dimWorld = GridView::dimensionworld,

        // equation indices
        conti0EqIdx = Indices::conti0EqIdx,
        momentum0EqIdx = Indices::momentum0EqIdx,

        // primary variable indices
        velocity0Idx = Indices::velocity0Idx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    StokesTestProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    { eps_ = 1e-6; }

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
    { return "stokes"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    template <class Context>   
    Scalar temperature(const Context &context,
                       int spaceIdx, int timeIdx) const
    { return 273.15 + 10; } // -> 10 deg C
    
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

        Scalar y = pos[1] - this->bboxMin()[1];
        Scalar height = this->bboxMax()[1] - this->bboxMin()[1];

        // parabolic velocity profile
        const Scalar maxVelocity = 1.0;

        Scalar a = - 4*maxVelocity/(height*height);
        Scalar b = - a*height;
        Scalar c = 0;
        
        DimVector velocity(0.0);
        velocity[0] = a * y*y + b * y + c;

        if (onRightBoundary_(pos))
            values.setOutFlow(context, spaceIdx, timeIdx);
        else if(onLeftBoundary_(pos)) {
            // left boundary is constraint!
            values = 0.0;
        }
        else {
            // top and bottom
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

        if (onLeftBoundary_(pos) || onRightBoundary_(pos)) {
            PrimaryVariables initCond;
            initial(initCond, context, spaceIdx, timeIdx);

            values.setConstraint(pressureIdx, conti0EqIdx, initCond[pressureIdx]);;
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
        const auto &pos = context.pos(spaceIdx, timeIdx);

        Scalar y = pos[1] - this->bboxMin()[1];
        Scalar height = this->bboxMax()[1] - this->bboxMin()[1];

        // parabolic velocity profile on boundaries
        const Scalar maxVelocity = 1.0;

        Scalar a = - 4*maxVelocity/(height*height);
        Scalar b = - a*height;
        Scalar c = 0;
        
        DimVector velocity(0.0);
        velocity[0] = a * y*y + b * y + c;
        
        for (int axisIdx = 0; axisIdx < dimWorld; ++axisIdx)
            values[velocity0Idx + axisIdx] = velocity[axisIdx];
        values[pressureIdx] = 1e5;
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
        return 
            onLeftBoundary_(pos) || onRightBoundary_(pos) ||
            onLowerBoundary_(pos) || onUpperBoundary_(pos);
    }

    Scalar eps_;
};

} //end namespace

#endif
