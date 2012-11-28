// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \ingroup RichardsModel
 *
 * \brief Contains the default definitions for the properties required
 *        by the Richards VCVF discretization.
 */
#ifndef EWOMS_RICHARDS_PROPERTY_DEFAULTS_HH
#define EWOMS_RICHARDS_PROPERTY_DEFAULTS_HH

#include "richardsmodel.hh"
#include "richardsindices.hh"
#include "richardsfluxvariables.hh"
#include "richardsratevector.hh"
#include "richardsboundaryratevector.hh"
#include "richardsprimaryvariables.hh"
#include "richardsvolumevariables.hh"
#include "richardsproperties.hh"
#include "richardsnewtonmethod.hh"

#include <ewoms/disc/vcfv/vcfvmultiphaseproblem.hh>
#include <ewoms/models/modules/velocity/vcfvvelocitymodules.hh>
#include <ewoms/material/components/nullcomponent.hh>
#include <ewoms/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <ewoms/material/heatconduction/dummyheatconductionlaw.hh>

namespace Ewoms {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////
//! Number of equations required by the model
SET_INT_PROP(VcfvRichards, NumEq, 1);
//! Number of fluid phases considered
SET_INT_PROP(VcfvRichards, NumPhases, 2);
//! Number of components considered
SET_INT_PROP(VcfvRichards, NumComponents, 2);

//! By default, assume that the first phase is the liquid one
SET_INT_PROP(VcfvRichards, LiquidPhaseIndex, 0);

//! The local residual operator
SET_TYPE_PROP(VcfvRichards,
              LocalResidual,
              RichardsLocalResidual<TypeTag>);

//! The global model used
SET_TYPE_PROP(VcfvRichards, Model, RichardsModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(VcfvRichards, BaseProblem, VcfvMultiPhaseProblem<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(VcfvRichards, RateVector, RichardsRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(VcfvRichards, BoundaryRateVector, RichardsBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(VcfvRichards, PrimaryVariables, RichardsPrimaryVariables<TypeTag>);

//! The class for the volume averaged quantities
SET_TYPE_PROP(VcfvRichards, VolumeVariables, RichardsVolumeVariables<TypeTag>);

//! The class for the quantities required for the flux calculation
SET_TYPE_PROP(VcfvRichards, FluxVariables, RichardsFluxVariables<TypeTag>);

//! The class of the Newton method
SET_TYPE_PROP(VcfvRichards, NewtonMethod, RichardsNewtonMethod<TypeTag>);

// disable the smooth upwinding method by default
SET_BOOL_PROP(VcfvRichards, EnableSmoothUpwinding, false);

//! The class with all index definitions for the model
SET_TYPE_PROP(VcfvRichards, Indices, Ewoms::RichardsIndices);

/*!
 * \brief Set type of the parameter objects for the material law
 *
 * By default this is just retrieved from the material law.
 */
SET_TYPE_PROP(VcfvRichards,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(VcfvRichards,
              HeatConductionLaw,
              Ewoms::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(VcfvRichards,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Use the Darcy relation by default
SET_TYPE_PROP(VcfvRichards, VelocityModule, Ewoms::VcfvDarcyVelocityModule<TypeTag>);

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. Please be aware that you
 * should be careful to use the Richards model in conjunction with
 * liquid non-wetting phases. This is only meaningful if the viscosity
 * of the liquid phase is _much_ lower than the viscosity of the
 * wetting phase.
 */
SET_PROP(VcfvRichards, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Ewoms::LiquidPhase<Scalar, Ewoms::NullComponent<Scalar> > type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. This doed not need to be
 * specified by the problem for the Richards model to work because the
 * Richards model does not conserve the non-wetting phase.
 */
SET_PROP(VcfvRichards, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Ewoms::GasPhase<Scalar, Ewoms::NullComponent<Scalar> > type;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the immiscible twophase fluid system. The
 * actual fluids used are specified using in the problem definition by
 * the WettingPhase and NonwettingPhase properties. Be aware that
 * using different fluid systems in conjunction with the Richards
 * model only makes very limited sense.
 */
SET_PROP(VcfvRichards, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Ewoms::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

} // end namepace Properties
} // end namespace Ewoms

#endif
