// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Copyright (C) 2012 by Katherina Baber                                   *
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Definition of a problem, for the 1p2c box problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_OUTFLOW_PROBLEM_HH
#define DUMUX_OUTFLOW_PROBLEM_HH

#include <dumux/material/fluidstates/compositionalfluidstate.hh>
#include <dumux/boxmodels/pvs/pvsproperties.hh>
#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{
template <class TypeTag>
class OutflowProblem;

namespace Properties
{
NEW_TYPE_TAG(OutflowBaseProblem);

// Set the grid type
SET_PROP(OutflowBaseProblem, Grid)
{
    typedef Dune::YaspGrid<2> type;
};

// Set the problem property
SET_TYPE_PROP(OutflowBaseProblem, Problem, Dumux::OutflowProblem<TypeTag>);

// Set fluid configuration
SET_PROP(OutflowBaseProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar> type;
};

// Disable gravity
SET_BOOL_PROP(OutflowBaseProblem, EnableGravity, false);

// Also write mass fractions to the output
SET_BOOL_PROP(OutflowBaseProblem, VtkWriteMassFractions, true);
}


/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxTestProblems
 *
 * \brief Definition of a problem, for the 1p2c box problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1m times 1m with a discretization length of 0.05m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where an outflow boundary condition is applied.
 * This problem uses the \ref OnePTwoCBoxModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1p2c -parameterFile ./test_1p2c.input</tt>
 */
template <class TypeTag>
class OutflowProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // copy some indices for convenience
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // component indices
        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx
    };


    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    OutflowProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(1e-6)
    {
        temperature_ = 273.15 + 20;
        FluidSystem::init(/*minT=*/temperature_ - 1, /*maxT=*/temperature_ + 2, /*numT=*/3,
                          /*minp=*/0.8e5, /*maxp=*/2.5e5, /*nump=*/500);
        
        // set parameters of porous medium
        perm_ = this->toDimMatrix_(1e-10);
        porosity_ = 0.4;
        tortuosity_ = 0.28;
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
    { return "outflow"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 36 degrees Celsius.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; } // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Evaluate the boundary conditions for a boundary segment.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        
        if (onLeftBoundary_(globalPos)) {
            Dumux::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);
            fs.setPressure(/*phaseIdx=*/0, fs.pressure(/*phaseIdx=*/0) + 1e5);

            fs.setMoleFraction(/*phaseIdx=*/0, N2Idx, 2e-5);
            fs.setMoleFraction(/*phaseIdx=*/0, H2OIdx, 1 - 2e-5);

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else if (onRightBoundary_(globalPos)) {
            Dumux::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;
            initialFluidState_(fs, context, spaceIdx, timeIdx);

            // impose an outflow boundary condition
            values.setOutFlow(context, spaceIdx, timeIdx, fs);
        }
        else
            // no flow on top and bottom
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
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
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
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        Dumux::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        values.assignNaive(fs);
    }

    // \}

    /*!
     * \brief Define the intrinsic permeability \f$[m^2]\f$.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    { return perm_; }

    /*!
     * \brief Define the porosity \f$[-]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    { return porosity_; }

    /*!
     * \brief Define the tortuosity \f$[?]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    template <class Context>
    Scalar tortuosity(const Context &context, int spaceIdx, int timeIdx) const
    { return tortuosity_; }

    /*!
     * \brief Define the dispersivity \f$[?]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    template <class Context>
    Scalar dispersivity(const Context &context,
                        int spaceIdx, int timeIdx) const
    { return 0; }

private:
    bool onLeftBoundary_(const GlobalPosition &pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition &pos) const
    { return pos[0] > this->bboxMax()[0] - eps_; }

    template <class FluidState, class Context>
    void initialFluidState_(FluidState &fs, 
                            const Context &context,
                            int spaceIdx, int timeIdx) const
    {
        Scalar T = temperature(context, spaceIdx, timeIdx);
        //Scalar rho = FluidSystem::H2O::liquidDensity(T, /*pressure=*/1.5e5);
        //Scalar z = context.pos(spaceIdx, timeIdx)[dim - 1] - this->bboxMax()[dim - 1];
        //Scalar z = context.pos(spaceIdx, timeIdx)[dim - 1] - this->bboxMax()[dim - 1];

        fs.setSaturation(/*phaseIdx=*/0, 1.0);
        fs.setPressure(/*phaseIdx=*/0, 1e5 /* + rho*z */);
        fs.setMoleFraction(/*phaseIdx=*/0, H2OIdx, 1.0);
        fs.setMoleFraction(/*phaseIdx=*/0, N2Idx, 0);
        fs.setTemperature(T);
    }

    const Scalar eps_;

    MaterialLawParams materialParams_;
    DimMatrix perm_;
    Scalar temperature_;
    Scalar porosity_;
    Scalar tortuosity_;
}; 
} //end namespace

#endif
