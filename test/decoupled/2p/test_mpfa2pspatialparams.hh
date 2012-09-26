// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Markus Wolff                                      *
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
#ifndef TEST_2P_SPATIALPARAMETERS_HH
#define TEST_2P_SPATIALPARAMETERS_HH

#include <dumux/decoupled/spatialparams/fvspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>


namespace Dumux
{

template<class TypeTag>
class Test2PSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(Test2PSpatialParams);

NEW_PROP_TAG(BackgroundEntryPressure);
NEW_PROP_TAG(LenseEntryPressure);
NEW_PROP_TAG(BackgroundPermeabilityXX);
NEW_PROP_TAG(BackgroundPermeabilityXY);
NEW_PROP_TAG(BackgroundPermeabilityYX);
NEW_PROP_TAG(BackgroundPermeabilityYY);
NEW_PROP_TAG(LensPermeabilityXX);
NEW_PROP_TAG(LensPermeabilityXY);
NEW_PROP_TAG(LensPermeabilityYX);
NEW_PROP_TAG(LensPermeabilityYY);
NEW_PROP_TAG(LensOneLowerLeftX);
NEW_PROP_TAG(LensOneUpperRightX);
NEW_PROP_TAG(LensTwoLowerLeftX);
NEW_PROP_TAG(LensTwoUpperRightX);
NEW_PROP_TAG(LensThreeLowerLeftX);
NEW_PROP_TAG(LensThreeUpperRightX);
NEW_PROP_TAG(LensOneLowerLeftY);
NEW_PROP_TAG(LensOneUpperRightY);
NEW_PROP_TAG(LensTwoLowerLeftY);
NEW_PROP_TAG(LensTwoUpperRightY);
NEW_PROP_TAG(LensThreeLowerLeftY);
NEW_PROP_TAG(LensThreeUpperRightY);
NEW_PROP_TAG(InletWidth);
NEW_PROP_TAG(InjectionFlux);
NEW_PROP_TAG(OutputInterval);
NEW_PROP_TAG(OutputTimeInterval);

// Set the spatial parameters
SET_TYPE_PROP(Test2PSpatialParams, SpatialParams, Dumux::Test2PSpatialParams<TypeTag>);

// Set the material law
SET_PROP(Test2PSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw> type;
};

SET_SCALAR_PROP(Test2PSpatialParams, BackgroundEntryPressure, 500.0);
SET_SCALAR_PROP(Test2PSpatialParams, LenseEntryPressure, 5000.0);
SET_SCALAR_PROP(Test2PSpatialParams, BackgroundPermeabilityXX, 1e-10);
SET_SCALAR_PROP(Test2PSpatialParams, BackgroundPermeabilityXY, 0.0);
SET_SCALAR_PROP(Test2PSpatialParams, BackgroundPermeabilityYX, 0.0);
SET_SCALAR_PROP(Test2PSpatialParams, BackgroundPermeabilityYY, 1e-10);
SET_SCALAR_PROP(Test2PSpatialParams, LensPermeabilityXX, 1e-14);
SET_SCALAR_PROP(Test2PSpatialParams, LensPermeabilityXY, 0.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensPermeabilityYX, 0.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensPermeabilityYY, 1e-14);
SET_SCALAR_PROP(Test2PSpatialParams, LensOneLowerLeftX, 7.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensOneUpperRightX, 13.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensTwoLowerLeftX, 2.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensTwoUpperRightX, 8.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensThreeLowerLeftX, 10.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensThreeUpperRightX, 3.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensOneLowerLeftY, 6.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensOneUpperRightY, 7.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensTwoLowerLeftY, 4.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensTwoUpperRightY, 5.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensThreeLowerLeftY, 2.0);
SET_SCALAR_PROP(Test2PSpatialParams, LensThreeUpperRightY, 3.0);
}

/** \todo Please doc me! */

template<class TypeTag>
class Test2PSpatialParams: public FVSpatialParams<TypeTag>
{
    typedef FVSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, Scalar, BackgroundEntryPressure, "The entry pressure of the background material [Pa]");
        REGISTER_PARAM(TypeTag, Scalar, LenseEntryPressure, "The entry pressure of the lenses [Pa]");
        REGISTER_PARAM(TypeTag, Scalar, BackgroundPermeabilityXX, "The xx-entry of the permebility tensor of the background medium [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, BackgroundPermeabilityXY, "The xy-entry of the permebility tensor of the background medium [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, BackgroundPermeabilityYX, "The yx-entry of the permebility tensor of the background medium [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, BackgroundPermeabilityYY, "The yy-entry of the permebility tensor of the background medium [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, LensPermeabilityXX, "The xx-entry of the permebility tensor of the lenses [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, LensPermeabilityXY, "The xy-entry of the permebility tensor of the lenses [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, LensPermeabilityYX, "The yx-entry of the permebility tensor of the lenses [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, LensPermeabilityYY, "The yy-entry of the permebility tensor of the lenses [m^2]");
        REGISTER_PARAM(TypeTag, Scalar, LensOneLowerLeftX, "The lower-left x-coordinate of the first lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensOneUpperRightX, "The upper-right x-coordinate of the first lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensTwoLowerLeftX, "The lower-left x-coordinate of the second lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensTwoUpperRightX, "The upper-right x-coordinate of the second lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensThreeLowerLeftX, "The lower-left x-coordinate of the third lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensThreeUpperRightX, "The upper-right x-coordinate of the third lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensOneLowerLeftY, "The lower-left y-coordinate of the first lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensOneUpperRightY, "The upper-right y-coordinate of the first lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensTwoLowerLeftY, "The lower-left y-coordinate of the second lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensTwoUpperRightY, "The upper-right y-coordinate of the second lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensThreeLowerLeftY, "The lower-left y-coordinate of the third lens [m]");
        REGISTER_PARAM(TypeTag, Scalar, LensThreeUpperRightY, "The upper-right y-coordinate of the third lens [m]");
    }

    const FieldMatrix& intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) const
    {
        if (isLensOne(globalPos))
            return permLenses_;
        else if (isLensTwo(globalPos))
            return permLenses_;
        else if (isLensThree(globalPos))
            return permLenses_;
        else
            return permBackground_;
    }

    Scalar porosity(const Element& element) const
    {
        return 0.4;
    }

    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (isLensOne(globalPos))
            return materialLawParamsLenses_;
        else if (isLensTwo(globalPos))
            return materialLawParamsLenses_;
        else if (isLensThree(globalPos))
            return materialLawParamsLenses_;
        else
            return materialLawParamsBackground_;
    }

    Test2PSpatialParams(const GridView& gridView) :
            ParentType(gridView), permBackground_(0), permLenses_(0),
            lensOneLowerLeft_(0), lensOneUpperRight_(0), lensTwoLowerLeft_(0), lensTwoUpperRight_(0), lensThreeLowerLeft_(0), lensThreeUpperRight_(0)
    {
        // residual saturations
        materialLawParamsBackground_.setSwr(0.);
        materialLawParamsBackground_.setSnr(0.);

        materialLawParamsLenses_.setSwr(0.);
        materialLawParamsLenses_.setSnr(0.);

        //parameters for Brooks-Corey law
        // entry pressures function
        materialLawParamsBackground_.setPe(GET_PARAM(TypeTag, Scalar, BackgroundEntryPressure));
        materialLawParamsLenses_.setPe(GET_PARAM(TypeTag, Scalar, LenseEntryPressure));

        // Brooks-Corey shape parameters
        materialLawParamsBackground_.setLambda(3);
        materialLawParamsLenses_.setLambda(2);

        permBackground_[0][0] = GET_PARAM(TypeTag, Scalar, BackgroundPermeabilityXX);
        permBackground_[0][1] = GET_PARAM(TypeTag, Scalar, BackgroundPermeabilityXY);
        permBackground_[1][0] = GET_PARAM(TypeTag, Scalar, BackgroundPermeabilityYX);
        permBackground_[1][1] = GET_PARAM(TypeTag, Scalar, BackgroundPermeabilityYY);
        permLenses_[0][0] = GET_PARAM(TypeTag, Scalar, LensPermeabilityXX);
        permLenses_[0][1] = GET_PARAM(TypeTag, Scalar, LensPermeabilityXY);
        permLenses_[1][0] = GET_PARAM(TypeTag, Scalar, LensPermeabilityYX);
        permLenses_[1][1] = GET_PARAM(TypeTag, Scalar, LensPermeabilityYY);
        lensOneLowerLeft_[0] = GET_PARAM(TypeTag, Scalar, LensOneLowerLeftX);
        lensOneUpperRight_[0] = GET_PARAM(TypeTag, Scalar, LensOneUpperRightX);
        lensTwoLowerLeft_[0] = GET_PARAM(TypeTag, Scalar, LensTwoLowerLeftX);
        lensTwoUpperRight_[0] = GET_PARAM(TypeTag, Scalar, LensTwoUpperRightX);
        lensThreeLowerLeft_[0] = GET_PARAM(TypeTag, Scalar, LensThreeLowerLeftX);
        lensThreeUpperRight_[0] = GET_PARAM(TypeTag, Scalar, LensThreeUpperRightX);
        lensOneLowerLeft_[1] = GET_PARAM(TypeTag, Scalar, LensOneLowerLeftY);
        lensOneUpperRight_[1] = GET_PARAM(TypeTag, Scalar, LensOneUpperRightY);
        lensTwoLowerLeft_[1] = GET_PARAM(TypeTag, Scalar, LensTwoLowerLeftY);
        lensTwoUpperRight_[1] = GET_PARAM(TypeTag, Scalar, LensTwoUpperRightY);
        lensThreeLowerLeft_[1] = GET_PARAM(TypeTag, Scalar, LensThreeLowerLeftY);
        lensThreeUpperRight_[1] = GET_PARAM(TypeTag, Scalar, LensThreeUpperRightY);
    }

private:

    bool isLensOne(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dim; i++)
        {
            if (globalPos[i] < lensOneLowerLeft_[i] || globalPos[i] > lensOneUpperRight_[i])
            {
                return false;
            }
        }
        return true;
    }
    bool isLensTwo(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dim; i++)
        {
            if (globalPos[i] < lensTwoLowerLeft_[i] || globalPos[i] > lensTwoUpperRight_[i])
            {
                return false;
            }
        }
        return true;
    }
    bool isLensThree(const GlobalPosition& globalPos) const
    {
        for (int i = 0; i < dim; i++)
        {
            if (globalPos[i] < lensThreeLowerLeft_[i] || globalPos[i] > lensThreeUpperRight_[i])
            {
                return false;
            }
        }
        return true;
    }

    MaterialLawParams materialLawParamsBackground_;
    MaterialLawParams materialLawParamsLenses_;
    FieldMatrix permBackground_;
    FieldMatrix permLenses_;
    GlobalPosition lensOneLowerLeft_;
    GlobalPosition lensOneUpperRight_;
    GlobalPosition lensTwoLowerLeft_;
    GlobalPosition lensTwoUpperRight_;
    GlobalPosition lensThreeLowerLeft_;
    GlobalPosition lensThreeUpperRight_;
};

} // end namespace
#endif
