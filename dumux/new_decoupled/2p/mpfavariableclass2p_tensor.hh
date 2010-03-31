// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
#ifndef DUNE_MPFAVARIABLECLASS2P_HH
#define DUNE_MPFAVARIABLECLASS2P_HH

#include <dumux/new_decoupled/2p/variableclass2p.hh>
#include <dumux/new_decoupled/2p/2pfluidstate.hh>
#include "2pproperties.hh"

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations
 * @author Markus Wolff
 */

namespace Dumux
{
/*!
 * \ingroup DecoupledModel
 */
//! Class including the variables and data of discretized data of the constitutive relations.
/*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 */
template<class TypeTag>
class MPFAVariableClass2P: public VariableClass2P<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef VariableClass2P<TypeTag> ParentType;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW, pn = Indices::pressureNW, pglobal = Indices::pressureGlobal,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhase = GET_PROP_VALUE(TypeTag, PTAG(
                NumPhases))
    };

    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation));

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;
    typedef Dune::FieldVector<FieldMatrix, numPhase> FieldMatrixVector;
    typedef Dune::FieldVector<FieldMatrixVector, 2*dim> CornerCellsFieldMatrixVector;
    typedef Dune::FieldVector<CornerCellsFieldMatrixVector, 2*dim> InterfaceCornerCellsFieldMatrixVector;

public:
    typedef Dune::BlockVector<InterfaceCornerCellsFieldMatrixVector> UpwindMobilitiesVector;

private:
    UpwindMobilitiesVector upwindMobilities_;//store lambda for efficiency reasons

public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */

    MPFAVariableClass2P(const GridView& gridView, Scalar& initialSat = *(new Scalar(1)),
            Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        ParentType(gridView, initialSat, initialVel)
    {
        //resize to grid size
        upwindMobilities_.resize(this->gridSize());

        upwindMobilities_ = CornerCellsFieldMatrixVector(FieldMatrixVector(FieldMatrix(0)));
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param codim codimension of the entity of which data has to be strored
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    MPFAVariableClass2P(const GridView& gridView, int codim, Scalar& initialSat = *(new Scalar(1)), Dune::FieldVector<
            Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim>(0))) :
        ParentType(gridView, codim, initialSat, initialVel)
    {
        //resize to grid size
        upwindMobilities_.resize(this->gridSize());

        upwindMobilities_ = CornerCellsFieldMatrixVector(FieldMatrixVector(FieldMatrix(0)));
    }

public:

    //! Return vector of wetting phase mobilities
    FieldMatrix& upwindMobilitiesWetting(int globalIdx, int faceIdx, int cellIdx)
    {
        return upwindMobilities_[globalIdx][faceIdx][cellIdx][wPhaseIdx];
    }

    const FieldMatrix& upwindMobilitiesWetting(int globalIdx, int faceIdx, int cellIdx) const
    {
        return upwindMobilities_[globalIdx][faceIdx][cellIdx][wPhaseIdx];
    }

    //! Return vector of non-wetting phase mobilities
    FieldMatrix& upwindMobilitiesNonwetting(int globalIdx, int faceIdx, int cellIdx)
    {
        return upwindMobilities_[globalIdx][faceIdx][cellIdx][nPhaseIdx];
    }

    const FieldMatrix& upwindMobilitiesNonwetting(int globalIdx, int faceIdx, int cellIdx) const
    {
        return upwindMobilities_[globalIdx][faceIdx][cellIdx][nPhaseIdx];
    }
};
}
#endif
