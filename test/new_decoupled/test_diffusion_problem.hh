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


#include "test_diffusion_spatialparams.hh"

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
SET_TYPE_PROP(DiffusionTestProblem, Problem, Dumux::TestDiffusionProblem<TTAG(DiffusionTestProblem)>);

SET_TYPE_PROP(DiffusionTestProblem, Model, Dumux::FVPressure2P<TTAG(DiffusionTestProblem)>);
//SET_TYPE_PROP(DiffusionTestProblem, Model, Dumux::MimeticPressure2P<TTAG(DiffusionTestProblem)>);
//SET_TYPE_PROP(DiffusionTestProblem, Model, Dumux::FVMPFAOPressure2P<TTAG(DiffusionTestProblem)>);

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
    typedef Dumux::TestDiffusionSpatialParams<TypeTag> type;
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
    TestDiffusionProblem(const GridView &gridView, const double delta = 1.0) :
        ParentType(gridView), delta_(delta)
    {
        this->variables().saturation() = 1.0;
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

    void init()
    {
        // set the initial condition of the model
        this->model().initial();

        // write the inital solution to disk
        this->writeCurrentResult_();
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
        double pi = 4.0*atan(1.0);
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];
        double ux = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        double uy = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);
        double kxx = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        double kxy = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        double kyy = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;
        double f0 = sin(pi*globalPos[0])*sin(pi*globalPos[1])*pi*pi*(1.0 + delta_)*(globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])
            + cos(pi*globalPos[0])*sin(pi*globalPos[1])*pi*(1.0 - 3.0*delta_)*globalPos[0]
            + cos(pi*globalPos[1])*sin(pi*globalPos[0])*pi*(1.0 - 3.0*delta_)*globalPos[1]
            + cos(pi*globalPos[1])*cos(pi*globalPos[0])*2.0*pi*pi*(1.0 - delta_)*globalPos[0]*globalPos[1];

        std::vector<double> result(2, 0.0);
        result[wPhaseIdx]=(f0 + 2.0*(globalPos[0]*(kxx*ux + kxy*uy) + globalPos[1]*(kxy*ux + kyy*uy)))/rt;

        return (result);
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return BoundaryConditions::dirichlet;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return (exact(globalPos));
    }

    typename BoundaryConditions::Flags bctypeSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return BoundaryConditions::dirichlet;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        return 1.0;
    }

    std::vector<Scalar> neumannPress(const GlobalPosition& globalPos, const Intersection& intersection) const
    {
        std::vector<Scalar> neumannFlux(2, 0.0);

        return neumannFlux;
    }

    Scalar exact (const GlobalPosition& globalPos) const
    {
        double pi = 4.0*atan(1.0);

        return (sin(pi*globalPos[0])*sin(pi*globalPos[1]));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar,dim> grad(0);
        double pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        grad[1] = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);

        return grad;
    }

private:
    double delta_;
};
} //end namespace

#endif
