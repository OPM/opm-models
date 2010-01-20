// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_TEST_2P_PROBLEM_HH
#define DUNE_TEST_2P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/material/fluids/water.hh>
#include <dumux/material/fluids/oil.hh>

#include <dumux/new_decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/new_decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/new_decoupled/2p/transport/fv/fvsaturation2p.hh>

#include "test_2p_soilproperties.hh"

namespace Dune
{

template<class TypeTag>
class Test2PProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(TwoPTestProblem, INHERITS_FROM(DecoupledTwoP, Transport));

// Set the grid type
SET_PROP(TwoPTestProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

// Set the problem property
SET_PROP(TwoPTestProblem, Problem)
{
public:
    typedef Dune::Test2PProblem<TTAG(TwoPTestProblem)> type;
};

// Set the model properties
SET_PROP(TwoPTestProblem, SaturationModel)
{
    typedef Dune::FVSaturation2P<TTAG(TwoPTestProblem)> type;
};
SET_PROP(TwoPTestProblem, PressureModel)
{
    typedef Dune::FVVelocity2P<TTAG(TwoPTestProblem)> type;
};

SET_INT_PROP(TwoPTestProblem, VelocityFormulation, TwoPCommonIndices::velocityW);

// Set the wetting phase
SET_TYPE_PROP(TwoPTestProblem, WettingPhase, Dune::Water);

// Set the non-wetting phase
SET_TYPE_PROP(TwoPTestProblem, NonwettingPhase, Dune::Oil);

// Set the soil properties
SET_PROP(TwoPTestProblem, Soil)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dune::Test2PSoil<Grid, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(TwoPTestProblem, EnableGravity, false);

SET_SCALAR_PROP(TwoPTestProblem, CFLFactor, 0.9);
}

/*!
 * \ingroup TwoPBoxProblems
 * \brief Soil decontamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically simulated until \f$t_{\text{end}} =
 * 50\,000\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 1\,000\;s\f$. 
 * To run the simulation execute  the following line in shell:
 * <tt>./test_2p 50000 1000</tt>
 */
template <class TypeTag = TTAG(TwoPTestProblem) >
class Test2PProblem : public IMPESProblem2P<TypeTag,
Test2PProblem<TypeTag> >
{
typedef Test2PProblem<TypeTag> ThisType;
typedef IMPESProblem2P<TypeTag, ThisType> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    wetting = TwoPCommonIndices::wPhase, nonwetting = TwoPCommonIndices::nPhase
};

typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
Test2PProblem(const GridView &gridView, const GlobalPosition LowerLeft = 0,
        const GlobalPosition UpperRight = 0)
: ParentType(gridView), Left_(LowerLeft[0]),
Right_(UpperRight[0])
{ }

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
{   return "test2p";}

bool shouldWriteRestartFile() const
{
    return false;
}

/*!
 * \brief Returns the temperature within the domain.
 *
 * This problem assumes a temperature of 10 degrees Celsius.
 */
Scalar temperature(const GlobalPosition& globalPos,
        const Element& element,
                const LocalPosition& localPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}


std::vector<Scalar> source(const GlobalPosition& globalPos,
        const Element& element,
        const LocalPosition& localPos)
{
    return std::vector<Scalar>(2,0.0);
}

typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
        const Intersection& intersection) const
{
    if ((globalPos[0] < eps_))
    return BoundaryConditions::dirichlet;
    // all other boundaries
    return BoundaryConditions::neumann;
}

BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos,
        const Intersection& intersection) const
{
    //        if (globalPos[0] > (Right_ - eps_) || globalPos[0] < eps_)
    if (globalPos[0] < eps_)
    return Dune::BoundaryConditions::dirichlet;
    else
    return Dune::BoundaryConditions::neumann;
}

Scalar dirichletPress(const GlobalPosition& globalPos,
        const Intersection& intersection) const
{
    if (globalPos[0] < eps_)
    return 2e5;
    // all other boundaries
    return 2e5;
}

Scalar dirichletSat(const GlobalPosition& globalPos,
        const Intersection& intersection) const
{
    if (globalPos[0] < eps_)
    return 0.8;
    // all other boundaries
    return 0;
}

std::vector<Scalar> neumannPress(const GlobalPosition& globalPos,
        const Intersection& intersection) const
{
    std::vector<Scalar> neumannFlux(2,0.0);
    if (globalPos[0]> Right_ - eps_)
    {
        neumannFlux[nonwetting] = 3e-4;
    }
    return neumannFlux;
}

Scalar neumannSat(const GlobalPosition& globalPos,
        const Intersection& intersection, Scalar factor) const
{
    if (globalPos[0]> Right_ - eps_)
    return factor;
    return 0;
}

Scalar initSat (const GlobalPosition& globalPos, const Element& element,
        const LocalPosition& localPos) const
{
    return 0.2;
}

private:
Scalar Left_;
Scalar Right_;

static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
