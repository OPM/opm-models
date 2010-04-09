// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#ifndef TEST_2P_SPATIALPARAMETERS_HH
#define TEST_2P_SPATIALPARAMETERS_HH


//#include <dumux/new_material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/new_material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/new_material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

/** \todo Please doc me! */

template<class TypeTag>
class Test2PSpatialParams
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid))     Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename Grid::ctype                            CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    typedef RegularizedBrooksCorey<Scalar>                RawMaterialLaw;
//    typedef LinearMaterial<Scalar>                        RawMaterialLaw;
public:
    typedef EffToAbsLaw<RawMaterialLaw>               MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    void update (Scalar saturationW, const Element& element)
    {

    }

    const FieldMatrix& intrinsicPermeability  (const GlobalPosition& globalPos, const Element& element) const
    {
        if (globalPos[0] >= innerLowerLeft_[0] && globalPos[0] <= innerUpperRight_[0] && globalPos[1]
            >= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1])
            return innerPermeability_;
        else
            return outerPermeability_;
    }

    double porosity(const GlobalPosition& globalPos, const Element& element) const
    {
        return 0.2;
    }


    // return the brooks-corey context depending on the position
    const MaterialLawParams& materialLawParams(const GlobalPosition& globalPos, const Element &element) const
    {
            return materialLawParams_;
    }


    Test2PSpatialParams(const GridView& gridView)
    : constPermeability_(0)
    {
        // residual saturations
        materialLawParams_.setSwr(0);
        materialLawParams_.setSnr(0);

        // parameters for the Brooks-Corey Law
        // entry pressures
        materialLawParams_.setPe(10000);

        // Brooks-Corey shape parameters
        materialLawParams_.setAlpha(2);


        //define a lense
        innerLowerLeft_[0] = 30;
        innerLowerLeft_[1] = 30;
        innerUpperRight_[0] = 60;
        innerUpperRight_[1] = 45;

        outerPermeability_[0][0] = outerPermeability_[1][1] = 1e-10;
        outerPermeability_[0][1] = outerPermeability_[1][0] = 0;

        innerPermeability_[0][0] = innerPermeability_[1][1] = 1e-12;
        innerPermeability_[0][1] = innerPermeability_[1][0] = 0;
    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix constPermeability_;
    FieldMatrix outerPermeability_;
    FieldMatrix innerPermeability_;
    GlobalPosition innerLowerLeft_;
    GlobalPosition innerUpperRight_;

};

} // end namespace
#endif
