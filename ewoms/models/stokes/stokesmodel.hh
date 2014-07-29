/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 *
 * \copydoc Ewoms::StokesModel
 */
#ifndef EWOMS_STOKES_MODEL_HH
#define EWOMS_STOKES_MODEL_HH

#include "stokesproperties.hh"
#include "stokeslocalresidual.hh"
#include "stokesproblem.hh"
#include "stokesindices.hh"
#include "stokeslocalresidual.hh"
#include "stokesmodel.hh"
#include "stokesintensivequantities.hh"
#include "stokesextensivequantities.hh"
#include "stokesboundaryratevector.hh"


#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/components/NullComponent.hpp>
#include <opm/material/heatconduction/FluidConduction.hpp>
#include <opm/material/fluidsystems/SinglePhaseFluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>

#include <ewoms/linear/superlubackend.hh>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class StokesModel;
}

namespace Opm {
namespace Properties {
//! The type tag for the problems using the Stokes equations
NEW_TYPE_TAG(StokesModel, INHERITS_FROM(VcfvDiscretization));

/*!
 * \brief The type tag for the problems using the Navier-Stokes equations.
 *
 * Basically this just takes everything from the Stokes model, but
 * sets the \c EnableNavierTerm property to true by default.
 */
NEW_TYPE_TAG(NavierStokesModel, INHERITS_FROM(StokesModel));

//! set the number of equations
SET_PROP(StokesModel, NumEq)
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

//! Increase the raw tolerance of the newton method to 10^-7
SET_SCALAR_PROP(StokesModel, NewtonRawTolerance, 1e-7);

#if HAVE_SUPERLU
SET_TAG_PROP(StokesModel, LinearSolverSplice, SuperLULinearSolver);
#else
#warning "No SuperLU installed. SuperLU is the recommended linear " \
         "solver for the Stokes models."
#endif

//! the Model property
SET_TYPE_PROP(StokesModel, Model, Ewoms::StokesModel<TypeTag>);

//! the IntensiveQuantities property
SET_TYPE_PROP(StokesModel, IntensiveQuantities, Ewoms::StokesIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(StokesModel, ExtensiveQuantities, Ewoms::StokesExtensiveQuantities<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(StokesModel, BoundaryRateVector, Ewoms::StokesBoundaryRateVector<TypeTag>);

//! The fluid system to use by default
SET_PROP(StokesModel, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::FluidSystems::SinglePhase<Scalar, Fluid> type;
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

namespace Ewoms {

/*!
 * \ingroup StokesModel
 *
 * \brief A model for the Navier-Stokes equations.
 *
 * This model implements Navier-Stokes flow of a single fluid. By
 * default, it solves the momentum balance of the time-dependent Stokes
 * equations, i.e.
 * \f[
 * \frac{\partial \left(\varrho\,\mathbf{v}\right)}
 *      {\partial t}
 * + \boldsymbol{\nabla} p
 * - \nabla \cdot
 * \left(
 * \mu \left(\boldsymbol{\nabla} \mathbf{v} + \boldsymbol{\nabla}
 *\mathbf{v}^T\right)
 * \right)
 * - \varrho\,\mathbf{g}
 * = 0\;,
 * \f]
 * and the mass balance equation
 * \f[
 * \frac{\partial \varrho}{\partial t}
 * + \nabla \cdot\left(\varrho\,\mathbf{v}\right)
 * - q
 * = 0 \;.
 * \f]
 *
 * If the property \c EnableNavierStokes is set to \c true, an
 * additional convective momentum flux term (Navier term) gets
 * included into the momentum conservation equations which allows to
 * capture turbolent flow regimes. This additional term is given by
 *  \f[
 * \varrho \left(\mathbf{v} \cdot \boldsymbol{\nabla} \right)
 * \mathbf{v} \;.
 * \f]
 *
 * These equations are discretized by a fully-coupled vertex-centered
 * finite volume scheme in space and using the implicit Euler method
 * in time. Be aware, that this discretization scheme is quite
 * unstable for the Navier-Stokes equations and quickly leads to
 * unphysical oscillations in the calculated solution. We intend to
 * use a more appropriate discretization scheme in the future, though.
 */
template <class TypeTag>
class StokesModel : public GET_PROP_TYPE(TypeTag, Discretization)
{
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dimWorld = GridView::dimensionworld };
    enum { phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex) };
    enum { numComponents = FluidSystem::numComponents };

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Ewoms::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    StokesModel(Simulator &simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Returns true iff a fluid phase is used by the model.
     *
     * \param phaseIdxQueried The index of the fluid phase in question
     */
    bool phaseIsConsidered(int phaseIdxQueried) const
    { return phaseIdxQueried == phaseIdx; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::ostringstream oss;
        if (pvIdx == Indices::pressureIdx)
            oss << "pressure";
        else if (Indices::moleFrac1Idx <= pvIdx && pvIdx < Indices::moleFrac1Idx + numComponents - 1)
            oss << "moleFraction^"
                << FluidSystem::componentName(pvIdx - Indices::moleFrac1Idx + 1);
        else if (Indices::velocity0Idx <= pvIdx && pvIdx < Indices::velocity0Idx + dimWorld)
            oss << "velocity_" << pvIdx - Indices::velocity0Idx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(int eqIdx) const
    {
        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "continuity^"
                << FluidSystem::componentName(eqIdx - Indices::conti0EqIdx);
        else if (Indices::momentum0EqIdx <= eqIdx &&
                 eqIdx < Indices::momentum0EqIdx + dimWorld)
            oss << "momentum_" << eqIdx - Indices::momentum0EqIdx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalDofIdx, int pvIdx) const
    {
        // for stokes flow the pressure gradients are often quite
        // small, so we need higher precision for pressure. TODO: find
        // a good weight for the pressure.
        if (Indices::pressureIdx == pvIdx) {
            return 1e-5;
        }

        return 1;
    }

    /*!
     * \copydoc FvBaseDiscretization::prepareOutputFields
     */
    void prepareOutputFields() const
    { }

    /*!
     * \copydoc FvBaseDiscretization::appendOutputFields
     */
    void appendOutputFields(BaseOutputWriter &writer) const
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&writer);
        if (!vtkWriter)
            return; // TODO (?): print warning

        typedef BaseOutputWriter::ScalarBuffer ScalarBuffer;
        typedef BaseOutputWriter::VectorBuffer VectorBuffer;

        // create the required scalar fields
        unsigned numVertices = this->gridView_.size(/*codim=*/dimWorld);
        ScalarBuffer &pressure = *vtkWriter->allocateManagedScalarBuffer(numVertices);
        ScalarBuffer &density = *vtkWriter->allocateManagedScalarBuffer(numVertices);
        ScalarBuffer &temperature = *vtkWriter->allocateManagedScalarBuffer(numVertices);
        ScalarBuffer &viscosity = *vtkWriter->allocateManagedScalarBuffer(numVertices);
        VectorBuffer &velocity = *vtkWriter->allocateManagedVectorBuffer(/*numOuter=*/numVertices,
                                                                         /*numInner=*/dimWorld);
        ScalarBuffer *moleFraction[numComponents];
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            moleFraction[compIdx] = vtkWriter->allocateManagedScalarBuffer(numVertices);

        // iterate over grid
        ElementContext elemCtx(this->simulator_);

        ElementIterator elemIt = this->gridView().template begin<0>();
        ElementIterator elemEndIt = this->gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updateAll(*elemIt);

            int numScv = elemCtx.numPrimaryDof(/*timeIdx=*/0);
            for (int dofIdx = 0; dofIdx < numScv; ++dofIdx)
            {
                int globalIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/dofIdx, /*timeIdx=*/0);
                const auto &intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/dofIdx, /*timeIdx=*/0);
                const auto &fluidState = intQuants.fluidState();

                pressure[globalIdx] = fluidState.pressure(phaseIdx);
                density[globalIdx] = fluidState.density(phaseIdx);
                temperature[globalIdx] = fluidState.temperature(phaseIdx);
                viscosity[globalIdx] = fluidState.viscosity(phaseIdx);

                // this seems to be a bug in Dune's dense vector
                // classes: It is possible to assign a DynamicVector
                // from a scalar and also to add a FieldVector to it,
                // but it is impossible to directly copy a FieldVector
                // into a DynamicVector...
                velocity[globalIdx] = 0;
                velocity[globalIdx] += intQuants.velocityCenter();

                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    (*moleFraction[compIdx])[globalIdx] = fluidState.moleFraction(phaseIdx, compIdx);
            };
        }

        std::ostringstream tmp;

        tmp.str("");
        tmp << "pressure_" << FluidSystem::phaseName(phaseIdx);
        vtkWriter->attachScalarVertexData(pressure, tmp.str());

        tmp.str("");
        tmp << "density_" << FluidSystem::phaseName(phaseIdx);
        vtkWriter->attachScalarVertexData(density, tmp.str());

        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            tmp.str("");
            tmp << "moleFraction_" << FluidSystem::phaseName(phaseIdx) << "^"
                << FluidSystem::componentName(compIdx);
            vtkWriter->attachScalarVertexData(*moleFraction[compIdx], tmp.str());
        }

        tmp.str("");
        tmp << "temperature_" << FluidSystem::phaseName(phaseIdx);
        vtkWriter->attachScalarVertexData(temperature, tmp.str());

        tmp.str("");
        tmp << "viscosity_" << FluidSystem::phaseName(phaseIdx);
        vtkWriter->attachScalarVertexData(viscosity, tmp.str());

        tmp.str("");
        tmp << "velocity_" << FluidSystem::phaseName(phaseIdx);
        vtkWriter->attachVectorVertexData(velocity, tmp.str());
    }
};
} // namespace Ewoms

#endif
