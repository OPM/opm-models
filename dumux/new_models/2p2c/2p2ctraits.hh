/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: andreas.lauser _at_ iws.uni-stuttgart.de                         *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#ifndef DUMUX_2P2CTRAITS_HH
#define DUMUX_2P2CTRAITS_HH

#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include "2p2cnewtoncontroller.hh"

namespace Dune
{
////////////////////////////////
// forward declarations
////////////////////////////////
template<class TypeTag>
class TwoPTwoCBoxModel;

template<class TypeTag>
class TwoPTwoCBoxJacobian;

template <class TypeTag>
class TwoPTwoCPnSwTraits;

template <class TypeTag>
class TwoPTwoCPwSnTraits;

template <class TypeTag>
class TwoPTwoCVertexData;

template <class TypeTag>
class TwoPTwoCElementData;

template <class TypeTag>
class TwoPTwoCFluxData;

/*!
 * \brief The indices for the isothermal TwoPTwoC model.
 */
template <int eqOffset = 0>
class TwoPTwoCIndices
{
    // Primary variable indices
    static const int pressureIdx = 0;     //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx   = 1;     //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    // Phase state (-> 'pseudo' primary variable)
    static const int nPhaseOnly = 0; //!< Only the non-wetting phase is present
    static const int wPhaseOnly = 1; //!< Only the wetting phase is present
    static const int bothPhases = 2;  //!< Both phases are present
    
    // Formulations
    static const int pWsN = 0; //!< Pw and Sn as primary variables
    static const int pNsW = 1;  //!< Pn and Sw as primary variables

    // Phase indices
    static const int wPhase      = 0;    //!< Index of the wetting phase in a phase vector
    static const int nPhase      = 1;    //!< Index of the non-wetting phase in a phase vector

    // Component indices
    static const int wComp       = 0; //!< Index of the wetting component in a component vector
    static const int nComp       = 1; //!< Index of the non-wetting component in a compent vector
    
    /*!
     * \brief Map a component index to a mass index.
     *
     * (The mass index is the index of a component in the result
     * vector of primary variables in the storage or flux terms.)
     */
    static int comp2MassIdx(int compIdx) { return compIdx; }

    /*!
     * \brief Component index of the principal component in a phase.
     */
    static int mainCompIdx(int phaseIdx) { return phaseIdx; }
};

////////////////////////////////
// properties
////////////////////////////////

namespace Properties
{
//! set the number of equations to 2
SET_PROP_INT(BoxTwoPTwoC, NumEq, 2);

//! Use the pw-Sn formulation by default
SET_PROP(BoxTwoPTwoC, TwoPTwoCTraits)
{
    typedef TwoPTwoCPwSnTraits<TypeTag> type;
};

//! Use the 2p2c local jacobian operator for the 2p2c model
SET_PROP(BoxTwoPTwoC, LocalJacobian)
{
    typedef TwoPTwoCBoxJacobian<TypeTag> type;
};

//! Use the 2p2c specific newton controller for the 2p2c model
SET_PROP(BoxTwoPTwoC, NewtonController)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))  NewtonMethod;
  
public:
    typedef TwoPTwoCNewtonController<NewtonMethod> type;
};

//! the Model property
SET_PROP(BoxTwoPTwoC, Model)
{
    typedef TwoPTwoCBoxModel<TypeTag> type;
};

//! the VertexData property
SET_PROP(BoxTwoPTwoC, VertexData)
{
    typedef TwoPTwoCVertexData<TypeTag> type;
};

//! the ElementData property
SET_PROP(BoxTwoPTwoC, ElementData)
{
    typedef TwoPTwoCElementData<TypeTag> type;
};

//! the FluxData property
SET_PROP(BoxTwoPTwoC, FluxData)
{
    typedef TwoPTwoCFluxData<TypeTag> type;
};

//! the default upwind factor. default 1.0, i.e. fully upwind...
SET_PROP(BoxTwoPTwoC, UpwindAlpha)
{
private:    
    typedef typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    static const Scalar value = 1.0;
};

//! the upwind factor for the mobility. uses the value of UpwindAlpha
//! if the property is not overwritten elsewhere
SET_PROP(BoxTwoPTwoC, MobilityUpwindAlpha)
{
private:    
    typedef typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    static const Scalar value = GET_PROP_VALUE(TypeTag, UpwindAlpha);
};

//! The number of equations/primary variables in the 2p2c model is 2
SET_INT_PROP(BoxTwoPTwoC, NumEq, 2);

//! The number of phases in the 2p2c model is 2
SET_INT_PROP(BoxTwoPTwoC, NumPhases, 2);

//! The number of components in the 2p2c model is 2
SET_INT_PROP(BoxTwoPTwoC, NumComponents, 2);

//! Set the default formulation to pWsN
SET_INT_PROP(BoxTwoPTwoC, Formulation, GET_PROP(TypeTag, PTAG(TwoPTwoC) );

}

/*!
 * \brief Generic 2P-2C traits.
 *
 * This class specifies the exact behaviour of the two-phase
 * two-component model. By using a different traits class for the
 * model, the model can change its behaviour considerably.
 */
template <class TypeTag>
class TwoPTwoCBaseTraits
{
public:
    static const int numEq = 2;         //!< Number of primary variables / equations
    static const int numPhases = 2;     //!< Number of fluid phases
    static const int numComponents = 2; //!< Number of fluid components within a phase
    
    // Primary variable indices
    static const int pressureIdx = 0;     //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx   = 1;     //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    // present phases
    static const int nPhaseOnly = 0; //!< Only the non-wetting phase is present
    static const int wPhaseOnly = 1; //!< Only the wetting phase is present
    static const int bothPhases = 2;  //!< Both phases are present
    
    // formulation
    static const int pWsN = 0; //!< Pw and Sn as primary variables
    static const int pNsW = 1;  //!< Pn and Sw as primary variables
};

/*!
 * \brief The traits for the pw-Sn formulation of the 2p2c model.
 */
template <class Scalar>
class TwoPTwoCPwSnTraits : public TwoPTwoCBaseTraits<Scalar>
{
    typedef TwoPTwoCBaseTraits<Scalar>     ParentT;

public:
    static const int formulation = ParentT::pWsN; //!< Formulation to use
};

/*!
 * \brief The traits for the pn-Sw formulation of the 2p2c model.
 */
template <class Scalar>
class TwoPTwoCPnSwTraits : public TwoPTwoCBaseTraits<Scalar>
{
    typedef TwoPTwoCBaseTraits<Scalar>     ParentT;

public:
    static const int formulation = ParentT::pNsW; //!< Formulation to use
    static const int wPhase      = 1;    //!< Index of the wetting phase in a phase vector
    static const int nPhase      = 0;    //!< Index of the non-wetting phase in a phase vector

    static const int wComp       = 1; //!< Index of the wetting component in a solution vector
    static const int nComp       = 0; //!< Index of the non-wetting component in a solution vector
};

/*!
 * \brief The traits for the non-isothermal 2p2c model formulation.
 */
template <class Scalar,
          class BaseTraits = TwoPTwoCPwSnTraits<Scalar> >
class TwoPTwoCNITraits : public BaseTraits
{
public:
    static const int numEq = 3;  //!< Override the mumber of primary variables: We also have temperature
    // Primary variable indices
    static const int temperatureIdx = 2; //! The index for temperature in solution vectors.
};

}

#endif
