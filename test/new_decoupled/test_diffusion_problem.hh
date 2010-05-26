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
#ifndef DUMUX_TEST_2P_PROBLEM_HH
#define DUMUX_TEST_2P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/new_material/components/water.hh>
#include <dumux/new_material/components/oil.hh>

#include <dumux/new_decoupled/2p/diffusion/diffusionproblem2p.hh>
#include <dumux/new_decoupled/2p/diffusion/fv/fvpressure2p.hh>
#include <dumux/new_decoupled/2p/diffusion/fvmpfa/fvmpfaopressure2p.hh>
#include <dumux/new_decoupled/2p/diffusion/mimetic/mimeticpressure2p.hh>


#include "test_2p_spatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestDiffusionProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(DiffusionTestProblem, INHERITS_FROM(DecoupledTwoP, MPFAProperties));

// Set the grid type
SET_PROP(DiffusionTestProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

// Set the problem property
SET_PROP(DiffusionTestProblem, Problem)
{
public:
    typedef Dumux::TestDiffusionProblem<TTAG(DiffusionTestProblem)> type;
};

SET_PROP(DiffusionTestProblem, Model)
{
    //    typedef Dumux::FVPressure2P<TTAG(DiffusionTestProblem)> type;
    typedef Dumux::MimeticPressure2P<TTAG(DiffusionTestProblem)> type;
    //    typedef Dumux::FVMPFAOPressure2P<TTAG(DiffusionTestProblem)> type;
};

//SET_INT_PROP(DiffusionTestProblem, VelocityFormulation,
//        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::velocityW);

SET_INT_PROP(DiffusionTestProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))::pressureGlobal);

// Set the wetting phase
SET_PROP(DiffusionTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Water<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(DiffusionTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Oil<Scalar> > type;
};

// Set the soil properties
SET_PROP(DiffusionTestProblem, SpatialParameters)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

public:
    typedef Dumux::Test2PSpatialParams<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(DiffusionTestProblem, EnableGravity, false);
}

/*!
* \ingroup DecoupledProblems
*/
template<class TypeTag = TTAG(DiffusionTestProblem)>
class TestDiffusionProblem: public DiffusionProblem2P<TypeTag, TestDiffusionProblem<TypeTag> >
{
    typedef TestDiffusionProblem<TypeTag> ThisType;
    typedef DiffusionProblem2P<TypeTag, ThisType> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    TestDiffusionProblem(const GridView &gridView, const GlobalPosition lowerLeft = 0, const GlobalPosition upperRight = 0) :
        ParentType(gridView), lowerLeft_(lowerLeft), upperRight_(upperRight)
    {
        this->variables().saturation()=1;
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
        return "testdiffusion";
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 10 degrees Celsius.
    */
    Scalar temperature(const GlobalPosition& globalPos, const Element& element) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}

    Scalar referencePressure(const GlobalPosition& globalPos, const Element& element) const
    {
        return 1e5; // -> 10°C
    }

    std::vector<Scalar> source(const GlobalPosition& globalPos, const Element& element)
        {
        return std::vector<Scalar>(2, 0.0);
        }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if ((globalPos[0] < lowerLeft_[0] + eps_) || globalPos[0] > upperRight_[0] - eps_)
            return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < lowerLeft_[0] + eps_ || globalPos[0] > upperRight_[0] - eps_)
            return Dumux::BoundaryConditions::dirichlet;
        else
            return Dumux::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        if (globalPos[0] < lowerLeft_[0] + eps_)
        {
            return 2e5;
        }
        // all other boundaries
        return 1.99e5;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return 1;
    }

    std::vector<Scalar> neumannPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        std::vector<Scalar> neumannFlux(2, 0.0);

        return neumannFlux;
    }

private:
    GlobalPosition lowerLeft_;
    GlobalPosition upperRight_;

    static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
