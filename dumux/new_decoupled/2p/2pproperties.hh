// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
/*!
 * \file
 *
 * \brief Defines the properties required for the twophase BOX model.
 */

#ifndef DUMUX_2PPROPERTIES_HH
#define DUMUX_2PPROPERTIES_HH

#include <dumux/new_decoupled/common/impesproperties.hh>
#include <dumux/new_decoupled/2p/transport/transportproperties.hh>

namespace Dune
{
/*!
 * \ingroup fracflow
 * \addtogroup DecoupledModel
 */
// \{

////////////////////////////////
// forward declarations
////////////////////////////////


template<class TypeTag>
class VariableClass2P;

template<class TypeTag>
class FluidSystem2P;

template<class TypeTag>
class TwoPPhaseState;

template <class TypeTag>
struct TwoPCommonIndices;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

/*!
 * \ingroup fracflow
 * \addtogroup DecoupledModel
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the two-phase problems
NEW_TYPE_TAG(DecoupledTwoP, INHERITS_FROM(IMPES, Transport));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG ( TwoPIndices );
NEW_PROP_TAG( SpatialParameters ); //!< The type of the soil properties object
NEW_PROP_TAG( EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( PressureFormulation); //!< The formulation of the model
NEW_PROP_TAG( SaturationFormulation); //!< The formulation of the model
NEW_PROP_TAG( VelocityFormulation); //!< The formulation of the model
NEW_PROP_TAG( EnableCompressibility);// ! Returns whether compressibility is allowed
NEW_PROP_TAG( WettingPhase); //!< The wetting phase for two-phase models
NEW_PROP_TAG( NonwettingPhase); //!< The non-wetting phase for two-phase models
NEW_PROP_TAG( FluidSystem );
NEW_PROP_TAG( PhaseState );

//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_PROP(DecoupledTwoP, NumPhases) //!< The number of phases in the 2p model is 2
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p model!");
};

SET_PROP(DecoupledTwoP, TwoPIndices)
{
  typedef TwoPCommonIndices<TypeTag> type;
};

//! Set the default formulation
SET_INT_PROP(DecoupledTwoP,
        PressureFormulation,
        TwoPCommonIndices<TypeTag>::pressureW);

SET_INT_PROP(DecoupledTwoP,
        SaturationFormulation,
        TwoPCommonIndices<TypeTag>::saturationW);

SET_INT_PROP(DecoupledTwoP,
        VelocityFormulation,
        TwoPCommonIndices<TypeTag>::velocityTotal);

SET_BOOL_PROP(DecoupledTwoP, EnableCompressibility, false);

SET_TYPE_PROP(DecoupledTwoP, Variables, VariableClass2P<TypeTag>);

SET_TYPE_PROP(DecoupledTwoP, FluidSystem, FluidSystem2P<TypeTag>);

SET_TYPE_PROP(DecoupledTwoP, PhaseState, TwoPPhaseState<TypeTag>);

// \}
}

/*!
 * \brief The common indices for the two-phase model.
 */
template <class TypeTag>
struct TwoPCommonIndices
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    // Formulations
    //saturation flags
    static const int saturationW = 0;
    static const int saturationNW = 1;
    //pressure flags
    static const int pressureW = 0;
    static const int pressureNW = 1;
    static const int pressureGlobal = 2;
    //velocity flags
    static const int velocityW = 0;
    static const int velocityNW = 1;
    static const int velocityTotal = 2;

    // Phase indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx; //!< Index of the wetting phase in a phase vector
    static const int nPhaseIdx = FluidSystem::nPhaseIdx; //!< Index of the non-wetting phase in a phase vector
};

// \}

}

#endif
