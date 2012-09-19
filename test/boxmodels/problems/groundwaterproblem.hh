// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2009 by Anneli Schoeniger                                 *
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
 *
 * \brief Test for the immisicible box model with only a single phase
 *
 * This problem is inspired by groundwater flow. Don't expect it to be
 * realistic, though.
 */
#ifndef DUMUX_GROUND_WATER_PROBLEM_HH
#define DUMUX_GROUND_WATER_PROBLEM_HH

#include <dumux/boxmodels/immiscible/immiscibleproperties.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{
template <class TypeTag>
class GroundWaterProblem;

namespace Properties
{
NEW_TYPE_TAG(GroundWaterBaseProblem);

NEW_PROP_TAG(LensLowerLeftX);
NEW_PROP_TAG(LensLowerLeftY);
NEW_PROP_TAG(LensLowerLeftZ);
NEW_PROP_TAG(LensUpperRightX);
NEW_PROP_TAG(LensUpperRightY);
NEW_PROP_TAG(LensUpperRightZ);
NEW_PROP_TAG(Permeability);
NEW_PROP_TAG(PermeabilityLens);

SET_PROP(GroundWaterBaseProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(GroundWaterBaseProblem, Grid, Dune::YaspGrid<2>);
//SET_TYPE_PROP(GroundWaterBaseProblem, Grid, Dune::SGrid<2, 2>);

SET_TYPE_PROP(GroundWaterBaseProblem, Problem, Dumux::GroundWaterProblem<TypeTag>);

SET_SCALAR_PROP(GroundWaterBaseProblem, LensLowerLeftX, 0.25);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensLowerLeftY, 0.25);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensLowerLeftZ, 0.25);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensUpperRightX, 0.75);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensUpperRightY, 0.75);
SET_SCALAR_PROP(GroundWaterBaseProblem, LensUpperRightZ, 0.75);
SET_SCALAR_PROP(GroundWaterBaseProblem, Permeability, 1e-10);
SET_SCALAR_PROP(GroundWaterBaseProblem, PermeabilityLens, 1e-12);
// Linear solver settings
SET_TYPE_PROP(GroundWaterBaseProblem, LinearSolverWrapper, Dumux::Linear::SolverWrapperCG<TypeTag> );
SET_TYPE_PROP(GroundWaterBaseProblem, PreconditionerWrapper, Dumux::Linear::PreconditionerWrapperILU<TypeTag> );
SET_INT_PROP(GroundWaterBaseProblem, LinearSolverVerbosity, 0);

// Enable gravity
SET_BOOL_PROP(GroundWaterBaseProblem, EnableGravity, true);
}

/*!
 * \ingroup BoxTestProblems
 * \brief Test for the immisicible box model with only a single phase
 *
 * This problem is inspired by groundwater flow. Don't expect it to be
 * realistic, though.
 */
template <class TypeTag>
class GroundWaterProblem 
    : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressure0Idx = Indices::pressure0Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    GroundWaterProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    {
        eps_ = 1.0e-3;

        lensLowerLeft_[0] = GET_PARAM(TypeTag, Scalar, LensLowerLeftX);
        if (dim > 1)
            lensLowerLeft_[1] = GET_PARAM(TypeTag, Scalar, LensLowerLeftY);
        if (dim > 2)
            lensLowerLeft_[2] = GET_PARAM(TypeTag, Scalar, LensLowerLeftY);

        lensUpperRight_[0] = GET_PARAM(TypeTag, Scalar, LensUpperRightX);
        if (dim > 1)
            lensUpperRight_[1] = GET_PARAM(TypeTag, Scalar, LensUpperRightY);
        if (dim > 2)
            lensUpperRight_[2] = GET_PARAM(TypeTag, Scalar, LensUpperRightY);

        intrinsicPerm_ = this->toDimMatrix_(GET_PARAM(TypeTag, Scalar, Permeability));
        intrinsicPermLens_ = this->toDimMatrix_(GET_PARAM(TypeTag, Scalar, PermeabilityLens));
    }

    /*! \brief Define the porosity.
     *
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    { return 0.4; }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    { return isInLens_(context.pos(spaceIdx, timeIdx))?intrinsicPermLens_:intrinsicPerm_; }

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
    { return "1ptest"; }

    /*!
     * \brief Return the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return 273.15 + 10; } // 10C


    template <class Context>
    void source(RateVector &values,
                const Context &context,
                int spaceIdx, int timeIdx) const
    {
        values = 0;
    }
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

        if (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos)) {
            Scalar pressure;
            Scalar T = temperature(context, spaceIdx, timeIdx);
            if (onLowerBoundary_(globalPos))
                pressure = 2e5;
            else // on upper boundary
                pressure = 1e5;

            Dumux::ImmiscibleFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> fs;
            fs.setSaturation(/*phaseIdx=*/0, 1.0);
            fs.setPressure(/*phaseIdx=*/0, pressure);
            fs.setTemperature(T);

            // impose an freeflow boundary condition
            values.setFreeFlow(context, spaceIdx, timeIdx, fs);
        }
        else {
            // no flow boundary
            values.setNoFlow();
        }
            
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        //const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        values[pressure0Idx] = 1.0e+5;// + 9.81*1.23*(20-globalPos[dim-1]);
    }

    // \}

private:
    bool onLowerBoundary_(const GlobalPosition &pos) const
    { return pos[dim-1] < eps_; }

    bool onUpperBoundary_(const GlobalPosition &pos) const
    { return pos[dim-1] > this->bboxMax()[dim-1] - eps_; }

    bool isInLens_(const GlobalPosition &pos) const
    {
        return 
            lensLowerLeft_[0] <= pos[0] && pos[0] <= lensUpperRight_[0] &&
            lensLowerLeft_[1] <= pos[1] && pos[1] <= lensUpperRight_[1];
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;
    
    DimMatrix intrinsicPerm_;
    DimMatrix intrinsicPermLens_;

    Scalar eps_;
};
} //end namespace

#endif
