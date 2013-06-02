// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Benjamin Faigle                              *
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
 * \copydoc Ewoms::TestMultTwoPTwoCProblem
 */
#ifndef EWOMS_TEST_2P2C_PROBLEM_HH
#define EWOMS_TEST_2P2C_PROBLEM_HH

#include "test_dec2p2c_spatialparams.hh"

#include <ewoms/decoupled/2p2c/2p2cproblem.hh>
#include <ewoms/decoupled/2p2c/fvpressure2p2cmultiphysics.hh>
#include <ewoms/decoupled/2p2c/fvtransport2p2cmultiphysics.hh>
#include <ewoms/decoupled/2p2c/cellData2p2cmultiphysics.hh>
#include <ewoms/material/fluidsystems/h2oairfluidsystem.hh>
#include <ewoms/io/cubegridcreator.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/common/fvector.hh>

namespace Ewoms
{

template<class TypeTag>
class TestMultTwoPTwoCProblem;

// Specify the properties
namespace Properties
{
NEW_TYPE_TAG(TestMultTwoPTwoCProblem, INHERITS_FROM(DecoupledTwoPTwoC, Test2P2CSpatialParams));

SET_TYPE_PROP(TestMultTwoPTwoCProblem, CellData, Ewoms::CellData2P2Cmultiphysics<TypeTag>);

SET_TYPE_PROP(TestMultTwoPTwoCProblem, GridCreator, Ewoms::CubeGridCreator<TypeTag>);

// Set the grid type
SET_PROP(TestMultTwoPTwoCProblem, Grid)
{
    //    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<3, 3> type;
};

// Set the problem property
SET_PROP(TestMultTwoPTwoCProblem, Problem)
{
    typedef Ewoms::TestMultTwoPTwoCProblem<TypeTag> type;
};

// Set the model properties
SET_PROP(TestMultTwoPTwoCProblem, TransportModel)
{
    typedef Ewoms::FVTransport2P2CMultiPhysics<TypeTag> type;
};

SET_PROP(TestMultTwoPTwoCProblem, PressureModel)
{
    typedef Ewoms::FVPressure2P2CMultiPhysics<TypeTag> type;
};

SET_INT_PROP(TestMultTwoPTwoCProblem, PressureFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::pressureNw);


//// Select fluid system
SET_PROP(TestMultTwoPTwoCProblem,
         FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Ewoms::FluidSystems::H2OAir<Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(TestMultTwoPTwoCProblem, EnableGravity, true);
SET_BOOL_PROP(TestMultTwoPTwoCProblem, EnableCapillarity, true);
SET_INT_PROP(TestMultTwoPTwoCProblem,
             BoundaryMobility,
             GET_PROP_TYPE(TypeTag, Indices)::satDependent);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, ImpetCFLFactor, 0.8);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, EndTime, 3000);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, InitialTimeStepSize, 200);
SET_INT_PROP(TestMultTwoPTwoCProblem, CellsX, 10);
SET_INT_PROP(TestMultTwoPTwoCProblem, CellsY, 10);
SET_INT_PROP(TestMultTwoPTwoCProblem, CellsZ, 10);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, DomainSizeX, 10.0);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, DomainSizeY, 10.0);
SET_SCALAR_PROP(TestMultTwoPTwoCProblem, DomainSizeZ, 10.0);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the sequential 2p2c model with multiphysics
 *
 * The domain is box shaped (3D). All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet). A Gas (Nitrogen)
 * is injected over a vertical well in the center of the domain.
 *
 * A multiphysics approach is used to adapt model complexity (see
 * description in the pressure module)
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_multiphyiscs2p2c</tt>
 * Optionally, simulation endtime and first timestep size can be
 * specified by programm arguments.
 */
template<class TypeTag>
class TestMultTwoPTwoCProblem: public IMPETProblem2P2C<TypeTag>
{
typedef IMPETProblem2P2C<TypeTag> ParentType;
typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

// boundary typedefs
typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

enum
{
    dim = GridView::dimension,
    dimWorld = GridView::dimensionworld
};

typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
typedef typename GridView::Intersection Intersection;
typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
TestMultTwoPTwoCProblem(TimeManager &timeManager) :
ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView()), eps_(1e-6), depthBOR_(1000.0)
{
    // Specifies how many time-steps are done before output will be written.
//    this->setOutputInterval(20);

    // initialize the tables of the fluid system
    FluidSystem::init(/*tempMin=*/280,
                      /*tempMax=*/290,
                      /*numTemp=*/10,
                      /*pMin=*/190000,
                      /*pMax=*/280000,
                      /*numP=*/400);
}

/*!
 * \name Problem parameters
 */
// \{

//! The problem name.
/*! This is used as a prefix for files generated by the simulation.
*/
const char *name() const
{
    return "test_multiphysics2p2c";
}
//!  Returns true if a restart file should be written.
/* The default behaviour is to write no restart file.
 */
bool shouldWriteRestartFile() const
{
    return false;
}

//! Returns the temperature within the domain.
/*! This problem assumes a temperature of 10 degrees Celsius.
 * \param globalPos The global Position
 */
Scalar temperatureAtPos(const GlobalPosition& globalPos) const
{
    return 273.15 + 10; // -> 10°C
}

// \}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::referencePressureAtPos()
 */
Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
{
    return 1e6;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::boundaryTypesAtPos()
 */
void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
{
    if (globalPos[0] > 10-1E-6 || globalPos[0] < 1e-6)
        bcTypes.setAllDirichlet();
    else
        // all other boundaries
        bcTypes.setAllNeumann();
}

/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::boundaryFormulation()
 */
const void boundaryFormulation(typename Indices::BoundaryFormulation &bcFormulation, const Intersection& intersection) const
{
    bcFormulation = Indices::concentration;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::dirichletAtPos()
 */
void dirichletAtPos(PrimaryVariables &bcValues ,const GlobalPosition& globalPos) const
{
    Scalar pRef = referencePressureAtPos(globalPos);
    Scalar temp = temperatureAtPos(globalPos);

    // Dirichlet for pressure equation
    bcValues[Indices::pressureEqIdx] = (globalPos[0] < 1e-6) ? (2.5e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1])
            : (2e5 - FluidSystem::H2O::liquidDensity(temp, pRef) * this->gravity()[dim-1]);

    // Dirichlet values for transport equations
    bcValues[Indices::contiWEqIdx] = 1.;
    bcValues[Indices::contiNEqIdx] = 1.- bcValues[Indices::contiWEqIdx];

}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::neumannAtPos()
 */
void neumannAtPos(PrimaryVariables &neumannValues, const GlobalPosition& globalPos) const
{
    neumannValues[Indices::contiNEqIdx] = 0.;
    neumannValues[Indices::contiWEqIdx] = 0.;
//    if (globalPos[1] < 15 && globalPos[1]> 5)
//    {
//        neumannValues[Indices::contiNEqIdx] = -0.015;
//    }
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::sourceAtPos()
 */
void sourceAtPos(PrimaryVariables &sourceValues, const GlobalPosition& globalPos) const
{
    this->setZero(sourceValues);
    if (std::abs(globalPos[0] - 4.8) < 0.5 && std::abs(globalPos[1] - 4.8) < 0.5)
        sourceValues[Indices::contiNEqIdx] = 0.0001;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::initialFormulation()
 */
const void initialFormulation(typename Indices::BoundaryFormulation &initialFormulation, const Element& element) const
{
    initialFormulation = Indices::concentration;
}
/*!
 * \copydoc Ewoms::TestDecTwoPTwoCProblem::initConcentrationAtPos()
 */
Scalar initConcentrationAtPos(const GlobalPosition& globalPos) const
{
    return 1;
}

private:
GlobalPosition lowerLeft_;
GlobalPosition upperRight_;

const Scalar eps_;
const Scalar depthBOR_;
};
} //end namespace

#endif
