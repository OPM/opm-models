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
#include <dumux/new_material/fluidmatrixinteractions/2p/absolutesaturationslaw.hh>

namespace Dune
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

    typedef RegularizedBrooksCoreyParams<Scalar>                RawMaterialLawParams;
//    typedef LinearMaterialParams<Scalar>                        RawMaterialLawParams;
    typedef AbsoluteSaturationsLawParams<RawMaterialLawParams> MaterialLawParams;
    typedef RegularizedBrooksCorey<MaterialLawParams>           RawMaterialLaw;
//    typedef LinearMaterial<MaterialLawParams>                   RawMaterialLaw;
public:
    //typedef RawMaterialLaw                                       MaterialLaw;
    typedef AbsoluteSaturationsLaw<RawMaterialLaw>               MaterialLaw;

    void update (Scalar saturationW, const Element& element)
    {

    }

    const FieldMatrix& intrinsicPermeability  (const GlobalPosition& globalPos, const Element& element) const
    {
        return constPermeability_;
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
        materialLawParams_.setSwr(0.2);
        materialLawParams_.setSnr(0.2);

        // parameters for the Brooks-Corey Law
        // entry pressures
        materialLawParams_.setPe(0);

        // Brooks-Corey shape parameters
        materialLawParams_.setAlpha(2);

        for(int i = 0; i < dim; i++)
        {
            constPermeability_[i][i] = 1e-7;
        }

    }

private:
    MaterialLawParams materialLawParams_;
    FieldMatrix constPermeability_;

};

} // end namespace
#endif
