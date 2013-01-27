// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
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
 * \copydoc Ewoms::MPFATwoPTestProblem
 */

#ifndef EWOMS_TEST_MPFA2P_PROBLEM_HH
#define EWOMS_TEST_MPFA2P_PROBLEM_HH

#include <dune/grid/alugrid/2d/alugrid.hh>

#include <ewoms/io/cubegridcreator.hh>

#include <ewoms/material/fluidsystems/liquidphase.hh>
#include <ewoms/material/components/simpleh2o.hh>
#include <ewoms/material/components/dnapl.hh>

#include <ewoms/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <ewoms/decoupled/2p/diffusion/fv/fvpressureproperties2padaptive.hh>
#include <ewoms/decoupled/2p/diffusion/fvmpfa/omethod/fvmpfaopressureproperties2p.hh>
#include <ewoms/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfalpressureproperties2p.hh>
#include <ewoms/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfalpressureproperties2padaptive.hh>
#include <ewoms/decoupled/2p/transport/fv/fvtransportproperties2p.hh>
#include <ewoms/decoupled/2p/impes/impesproblem2p.hh>

#include<ewoms/decoupled/2p/transport/fv/evalcflfluxcoats.hh>
#include<ewoms/decoupled/2p/impes/gridadaptionindicator2plocal.hh>

#include <dune/common/fvector.hh>

#include <sstream>
#include <string>

#include "test_mpfa2pspatialparams.hh"

namespace Ewoms
{

template<class TypeTag>
class MPFATwoPTestProblem;

//////////
// Specify the properties
//////////
namespace Properties
{

NEW_TYPE_TAG(MPFATwoPTestProblem, INHERITS_FROM(Test2PSpatialParams));

NEW_PROP_TAG(OutputfileName);
NEW_PROP_TAG(ModelType);

// Set the grid type
SET_PROP(MPFATwoPTestProblem, Grid)
{
    typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> type;
};

// set the GridCreator property
SET_TYPE_PROP(MPFATwoPTestProblem, GridCreator, CubeGridCreator<TypeTag>);


// Set the problem property
SET_PROP(MPFATwoPTestProblem, Problem)
{
public:
    typedef Ewoms::MPFATwoPTestProblem<TypeTag> type;
};

// Set the wetting phase
SET_PROP(MPFATwoPTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Ewoms::LiquidPhase<Scalar, Ewoms::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(MPFATwoPTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Ewoms::LiquidPhase<Scalar, Ewoms::DNAPL<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(MPFATwoPTestProblem, EnableGravity, true);

SET_TYPE_PROP(MPFATwoPTestProblem, EvalCflFluxFunction, Ewoms::EvalCflFluxCoats<TypeTag>);
SET_SCALAR_PROP(MPFATwoPTestProblem, ImpetCflFactor, 1.0);
SET_TYPE_PROP(MPFATwoPTestProblem, GridAdaptIndicator, Ewoms::GridAdaptionIndicator2PLocal<TypeTag>);

// grid adaption parameters
SET_BOOL_PROP(MPFATwoPTestProblem, GridAdaptEnableInitializationIndicator, true);
SET_BOOL_PROP(MPFATwoPTestProblem, GridAdaptRefineAtDirichletBC, false);
SET_BOOL_PROP(MPFATwoPTestProblem, GridAdaptRefineAtFluxBC, true);
SET_BOOL_PROP(MPFATwoPTestProblem, GridAdaptRefineAtSource, false);

SET_SCALAR_PROP(MPFATwoPTestProblem, DomainSizeX, 20.0);
SET_SCALAR_PROP(MPFATwoPTestProblem, DomainSizeY, 10.0);
SET_SCALAR_PROP(MPFATwoPTestProblem, DomainSizeZ, 1.0);

SET_INT_PROP(MPFATwoPTestProblem, CellsX, 10);
SET_INT_PROP(MPFATwoPTestProblem, CellsY, 10);
SET_INT_PROP(MPFATwoPTestProblem, CellsZ, 1);
SET_SCALAR_PROP(MPFATwoPTestProblem, EndTime, 5e4);
SET_SCALAR_PROP(MPFATwoPTestProblem, InletWidth, 2.0);
SET_SCALAR_PROP(MPFATwoPTestProblem, InjectionFlux, 0.1);
SET_SCALAR_PROP(MPFATwoPTestProblem, OutputInterval, 0.0);
SET_SCALAR_PROP(MPFATwoPTestProblem, OutputTimeInterval, 1e4);
SET_STRING_PROP(MPFATwoPTestProblem, OutputfileName, "test_mpfa2p");
SET_STRING_PROP(MPFATwoPTestProblem, ModelType, "MPFAL");

NEW_TYPE_TAG(FVTwoPTestProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(FVAdaptiveTwoPTestProblem, INHERITS_FROM(FVPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFAOTwoPTestProblem, INHERITS_FROM(FVMPFAOPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFALTwoPTestProblem, INHERITS_FROM(FVMPFALPressureTwoP, FVTransportTwoP, IMPESTwoP, MPFATwoPTestProblem));
NEW_TYPE_TAG(MPFALAdaptiveTwoPTestProblem, INHERITS_FROM(FVMPFALPressureTwoPAdaptive, FVTransportTwoP, IMPESTwoPAdaptive, MPFATwoPTestProblem));

}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for sequential 2p models
 *
 * A DNAPL is injected from the top into a rectangular 2D domain saturated by water.
 * The remaining upper and the lower boundary is closed (Neumann = 0). At the sides a hydrostatic pressure condition
 * and free outflow for saturation are set. The domain is heterogeneous with a backround material and three lenses.
 *
 * To run the simulation execute the following line in shell:
 *
 * <tt>./test_mpfa2p</tt>,
 *
 * Additionally, the numerical model can be switched by executing with the parameter "ModelType":
 *
 * <tt>./test_mpfa2p --ModelType=...</tt>,
 *
 * where ModelType can be:
 * - FV (standard finite volume)
 * - FVAdaptive (adaptive finite volume)
 * - MPFAO (MPFA o-method)
 * - MPFAL (MPFA l-method)
 * - MPFALAdaptive (adaptive MPFA l-method)
 */
template<class TypeTag>
class MPFATwoPTestProblem: public IMPESProblem2P<TypeTag>
{
typedef IMPESProblem2P<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;

typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

enum
{
    dim = GridView::dimension, dimWorld = GridView::dimensionworld
};

enum
{
    nPhaseIdx = Indices::nPhaseIdx,
    pWIdx = Indices::pwIdx,
    SwIdx = Indices::SwIdx,
    eqIdxPress = Indices::pressEqIdx,
    eqIdxSat = Indices::satEqIdx
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
MPFATwoPTestProblem(TimeManager &timeManager)
    : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
{
    this->setGrid(GridCreator::grid());

    Scalar inletWidth = GET_PARAM(TypeTag, Scalar, InletWidth);
    GlobalPosition inletCenter = this->bboxMax();
    inletCenter[0] *= 0.5;

    inletLeftCoord_ = inletCenter;
    inletLeftCoord_[0] -=0.5*inletWidth;
    inletRightCoord_ = inletCenter;
    inletRightCoord_[0] +=0.5*inletWidth;

    inFlux_ = GET_PARAM(TypeTag, Scalar, InjectionFlux);

    int outputInterval = GET_PARAM(TypeTag, int, OutputInterval);
    this->setOutputInterval(outputInterval);

    Scalar outputTimeInterval = GET_PARAM(TypeTag, Scalar, OutputTimeInterval);
    this->setOutputTimeInterval(outputTimeInterval);
}

static void registerParameters()
{
    ParentType::registerParameters();
    SpatialParams::registerParameters();

    REGISTER_PARAM(TypeTag, Scalar, InletWidth, "The width on which mass is injected [m]" );
    REGISTER_PARAM(TypeTag, Scalar, InjectionFlux, "The rate of the injected mass [kg/m^2]");
    REGISTER_PARAM(TypeTag, bool, OutputInterval, "The number of time steps between writing VTK output");
    REGISTER_PARAM(TypeTag, bool, OutputTimeInterval, "The minumum time [s] between writing VTK output");
    REGISTER_PARAM(TypeTag, std::string, OutputfileName, "The name of the VTK output file");
    REGISTER_PARAM(TypeTag, std::string, ModelType, "The PDE model to be used (FV, FVAdaptive, MPFAL, MPFAO, MPFALAdaptive");
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
{
    std::ostringstream oss;
    oss << GET_PARAM(TypeTag, std::string, OutputfileName);

    std::string modelType = GET_PARAM(TypeTag, std::string, ModelType);
    if (modelType == "FV")
        oss << "_tpfa";
    else if (modelType == "FVAdaptive")
        oss << "_tpfa_adaptive";
    else if (modelType == "MPFAO")
        oss << "_mpfao";
    else if (modelType == "MPFAL")
        oss << "_mpfal";
    else if (modelType == "MPFALAdaptive")
        oss << "_mpfal_adaptive";

    return oss.str();
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
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}

//! Returns the reference pressure for evaluation of constitutive relations
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e5; // -> 10°C
}

void source(PrimaryVariables &values,const Element& element) const
{
    values = 0;
}

/*!
* \brief Returns the type of boundary condition.
*
* BC for pressure equation can be dirichlet (pressure) or neumann (flux).
*
* BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
*/
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (isInlet(globalPos))
    {
        bcTypes.setNeumann(eqIdxPress);
        bcTypes.setDirichlet(eqIdxSat);
    }
    else if (isBottom(globalPos) || isTop(globalPos))
    {
        bcTypes.setAllNeumann();
    }
    else
    {
        bcTypes.setDirichlet(eqIdxPress);
        bcTypes.setOutflow(eqIdxSat);
    }
}

//! set dirichlet condition  (pressure [Pa], saturation [-])
void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    values[pWIdx] = (1e5 - (this->bboxMax()- globalPos) * this->gravity() * WettingPhase::density(temp, pRef));
    values[SwIdx] = 1.0;

    if (isInlet(globalPos))
    {
        values[SwIdx] = 0.0;
    }
}

//! set neumann condition for phases (flux, [kg/(m^2 s)])
void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
{
    values = 0;
    if (isInlet(globalPos))
    {
        values[nPhaseIdx] = -inFlux_;
    }

}

//! return initial solution -> only saturation values have to be given!
void initialAtPos(PrimaryVariables &values,
        const GlobalPosition& globalPos) const
{
    values[pWIdx] = 0;
    values[SwIdx] = 1.0;
}

private:

bool isInlet(const GlobalPosition& globalPos) const
{
        if (!isTop(globalPos))
            return false;

    for (int i = 0; i < dim; i++)
    {
        if (globalPos[i] < inletLeftCoord_[i] - eps_)
            return false;
        if (globalPos[i] > inletRightCoord_[i] + eps_)
            return false;
    }
    return true;
}

bool isTop(const GlobalPosition& globalPos) const
{
    if (dim == 2)
    {
        if (globalPos[1] > this->bboxMax()[1] - eps_)
            return true;
    }
    if (dim == 3)
    {
        if (globalPos[2] > this->bboxMax()[2] - eps_)
            return true;
    }
    return false;
}

bool isBottom(const GlobalPosition& globalPos) const
{
    if (dim == 2)
    {
        if (globalPos[1] < eps_)
            return true;
    }
    if (dim == 3)
    {
        if (globalPos[2] < eps_)
            return true;
    }
    return false;
}

Scalar inFlux_;
GlobalPosition inletLeftCoord_;
GlobalPosition inletRightCoord_;
static const Scalar eps_;
};

template <class TypeTag>
const typename MPFATwoPTestProblem<TypeTag>::Scalar
MPFATwoPTestProblem<TypeTag>::eps_ = 1e-6;

} //end namespace

#endif
