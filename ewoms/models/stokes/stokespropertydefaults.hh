// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser
  Copyright (C) 2012 by Klaus Mosthaf

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \ingroup StokesModel
 *
 * \brief Defines default values for the properties required by the
 *        Stokes model.
 */
#ifndef EWOMS_STOKES_PROPERTY_DEFAULTS_HH
#define EWOMS_STOKES_PROPERTY_DEFAULTS_HH

#include "stokesproperties.hh"
#include "stokesindices.hh"
#include "stokeslocalresidual.hh"
#include "stokesmodel.hh"
#include "stokesvolumevariables.hh"
#include "stokesfluxvariables.hh"
#include "stokesboundaryratevector.hh"

#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/components/NullComponent.hpp>
#include <opm/material/heatconduction/FluidConduction.hpp>
#include <opm/material/fluidsystems/1pFluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <ewoms/linear/superlubackend.hh>

namespace Opm {
namespace Properties {
SET_PROP(StokesModel, NumEq) //!< set the number of equations
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    static const int value = GridView::dimensionworld
                             + FluidSystem::numComponents;
};

//! the number of phases
SET_INT_PROP(StokesModel, NumPhases, 1);

//! the number of components
SET_INT_PROP(StokesModel, NumComponents, 1);

//! Use the Stokes local residual function for the Stokes model
SET_TYPE_PROP(StokesModel,
              LocalResidual,
              Ewoms::StokesLocalResidual<TypeTag>);

//! Use the Stokes local residual function for the Stokes model
SET_TYPE_PROP(StokesModel,
              BaseProblem,
              Ewoms::StokesProblem<TypeTag>);

//! Increase the relative tolerance of the newton method to 10^-7
SET_SCALAR_PROP(StokesModel, NewtonRelativeTolerance, 1e-7);

#if HAVE_SUPERLU
SET_TAG_PROP(StokesModel, LinearSolverSplice, SuperLULinearSolver);
#else
#warning "No SuperLU installed. SuperLU is the recommended linear " \
         "solver for the Stokes models."
#endif

//! the Model property
SET_TYPE_PROP(StokesModel, Model, Ewoms::StokesModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(StokesModel, VolumeVariables, Ewoms::StokesVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(StokesModel, FluxVariables, Ewoms::StokesFluxVariables<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(StokesModel, BoundaryRateVector, Ewoms::StokesBoundaryRateVector<TypeTag>);

//! The fluid system to use by default
SET_PROP(StokesModel, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::FluidSystems::OneP<Scalar, Fluid> type;
};

//! The fluid that is used in the single-phase fluidsystem.
SET_PROP(StokesModel, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

SET_TYPE_PROP(StokesModel, Indices, Ewoms::StokesIndices<TypeTag, /*PVOffset=*/0>);

//! Choose the type of the employed fluid state.
SET_PROP(StokesModel, FluidState)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Opm::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! set the heat conduction law to a dummy one by default
SET_TYPE_PROP(StokesModel,
              HeatConductionLaw,
              Opm::FluidHeatConduction<typename GET_PROP_TYPE(TypeTag, FluidSystem),
                                         typename GET_PROP_TYPE(TypeTag, Scalar),
                                         GET_PROP_VALUE(TypeTag, StokesPhaseIndex)>);

//! extract the type parameter objects for the heat conduction law
//! from the law itself
SET_TYPE_PROP(StokesModel,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

//! Set the phaseIndex per default to zero (important for two-phase fluidsystems).
SET_INT_PROP(StokesModel, StokesPhaseIndex, 0);

//! Disable the energy equation by default
SET_BOOL_PROP(StokesModel, EnableEnergy, false);

//! Disable the inertial term for the Stokes model by default
SET_BOOL_PROP(StokesModel, EnableNavierTerm, false);

//! The (Navier-)Stokes needs the gradients at the center of the SCVs, so
//! we enable them here.
SET_BOOL_PROP(StokesModel, RequireScvCenterGradients, true);

//! Enable the inertial term for the Navier-Stokes model
SET_BOOL_PROP(NavierStokesModel, EnableNavierTerm, true);
}} // namespace Opm,Properties

#endif
