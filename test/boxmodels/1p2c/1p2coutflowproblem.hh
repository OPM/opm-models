// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Katherina Baber                                   *
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Definition of a problem, for the 1p2c box problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_OUTFLOW_PROBLEM_HH
#define DUMUX_1P2C_OUTFLOW_PROBLEM_HH

#include <dumux/boxmodels/1p2c/1p2cmodel.hh>
#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{
template <class TypeTag>
class OnePTwoCOutflowProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCOutflowProblem, INHERITS_FROM(BoxOnePTwoC));

// Set the grid type
SET_PROP(OnePTwoCOutflowProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_TYPE_PROP(OnePTwoCOutflowProblem, Problem, Dumux::OnePTwoCOutflowProblem<TypeTag>);

// Set fluid configuration
SET_PROP(OnePTwoCOutflowProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar> type;
};

// Disable gravity
SET_BOOL_PROP(OnePTwoCOutflowProblem, EnableGravity, false);

// Also write mass fractions to the output
SET_BOOL_PROP(OnePTwoCOutflowProblem, VtkWriteMassFractions, true);
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
class OnePTwoCOutflowProblem
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, OnePTwoCIndices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,

        // indices of the equations
        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePTwoCOutflowProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
        , eps_(1e-6)
    {
        // set parameters of porous medium
        perm_ = 2.142e-11;
        porosity_ = 0.13;
        tortuosity_ = 0.706;

        // calculate the injection volume
        totalInjectionVolume_ = 0;
        FVElementGeometry fvGeom;
        ElementIterator elemIt = this->gridView().template begin<0>();
        const ElementIterator endIt = this->gridView().template end<0>();
        for (; elemIt != endIt; ++ elemIt) {
            fvGeom.update(this->gridView(), *elemIt);
            for (int i = 0; i < fvGeom.numVertices; ++i) {
                const GlobalPosition &pos = fvGeom.subContVol[i].global;
                if (inInjectionVolume_(pos))
                    totalInjectionVolume_ += fvGeom.subContVol[i].volume;
            };
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
    { return "outflow"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 36 degrees Celsius.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 273.15 + 36; } // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    template <class Context>
    void boundaryTypes(BoundaryTypes &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        values.setAllDirichlet();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    template <class Context>
    void dirichlet(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    template <class Context>
    void neumann(RateVector &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    {
        values = 0;
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
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);

        values = Scalar(0.0);

        if (inInjectionVolume_(globalPos)) {
            // total volumetric injection rate in ml/h
            Scalar injRateVol = 0.1;
            // convert to m^3/s
            injRateVol *= 1e-6/3600;
            // total mass injection rate. assume a density of 1030kg/m^3
            Scalar injRateMass = injRateVol*1030.0;

            // trail concentration in injected fluid in [mol/ml]
            Scalar trailInjRate = 1e-5;
            // convert to mol/m^3
            trailInjRate *= 1e6;
            // convert to mol/s
            trailInjRate *= injRateVol;

            // source term of the total mass
            values[contiEqIdx] = injRateMass / totalInjectionVolume_; // [kg/(s*m^3)]
            values[transEqIdx] = trailInjRate / totalInjectionVolume_; // [mol/(s*m^3)]
        }
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos
            = context.pos(spaceIdx, timeIdx);

        initial_(values, globalPos);
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
    Scalar intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        //context.pos(spaceIdx, timeIdx);
        return perm_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        //const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        return porosity_;
    }

    /*!
     * \brief Define the tortuosity \f$[?]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    template <class Context>
    Scalar tortuosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        //const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        return tortuosity_;
    }

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
    bool inInjectionVolume_(const GlobalPosition &globalPos) const
    {
        return
            10e-3 < globalPos[0] && globalPos[0] < 12e-3 &&
            10e-3 < globalPos[1] && globalPos[1] < 12e-3;
    };


    // the internal method for the initial condition
    void initial_(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = 0.0; //initial condition for the pressure
        values[x1Idx] = 0.0; //initial condition for the trail molefraction
//        if(globalPos[0] > 0.4 && globalPos[0] < 0.6 && globalPos[1] > 0.4 && globalPos[1] < 0.6)
//            values[x1Idx] = 0.6;
    }

    const Scalar eps_;

    Scalar perm_;
    Scalar porosity_;
    Scalar tortuosity_;
    Scalar totalInjectionVolume_;
}; 
} //end namespace

#endif
