/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \copydoc Ewoms::ImmiscibleModel
 */
#ifndef EWOMS_IMMISCIBLE_MODEL_HH
#define EWOMS_IMMISCIBLE_MODEL_HH

#include <ewoms/parallel/mpihelper.hh>
#include "immiscibleproperties.hh"
#include "immiscibleindices.hh"
#include "immisciblefluxvariables.hh"
#include "immiscibleprimaryvariables.hh"
#include "immisciblevolumevariables.hh"
#include "immiscibleratevector.hh"
#include "immiscibleboundaryratevector.hh"
#include "immisciblelocalresidual.hh"

#include <ewoms/models/common/multiphasebasemodel.hh>
#include <ewoms/models/common/energymodule.hh>
#include <ewoms/io/vtkenergymodule.hh>
#include <opm/material/components/NullComponent.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/fluidsystems/1pFluidSystem.hpp>
#include <opm/material/fluidsystems/2pImmiscibleFluidSystem.hpp>

#include <dune/common/unused.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class ImmiscibleModel;
}

namespace Opm {
namespace Properties {
//! The generic type tag for problems using the immiscible multi-phase model
NEW_TYPE_TAG(ImmiscibleModel, INHERITS_FROM(MultiPhaseBaseModel, VtkEnergy));
//! The type tag for single-phase immiscible problems
NEW_TYPE_TAG(ImmiscibleOnePhaseModel, INHERITS_FROM(ImmiscibleModel));
//! The type tag for two-phase immiscible problems
NEW_TYPE_TAG(ImmiscibleTwoPhaseModel, INHERITS_FROM(ImmiscibleModel));

//! Use the immiscible multi-phase local jacobian operator for the immiscible multi-phase model
SET_TYPE_PROP(ImmiscibleModel, LocalResidual,
              Ewoms::ImmiscibleLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(ImmiscibleModel, Model, Ewoms::ImmiscibleModel<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(ImmiscibleModel, RateVector, Ewoms::ImmiscibleRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(ImmiscibleModel, BoundaryRateVector, Ewoms::ImmiscibleBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(ImmiscibleModel, PrimaryVariables, Ewoms::ImmisciblePrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(ImmiscibleModel, VolumeVariables, Ewoms::ImmiscibleVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(ImmiscibleModel, FluxVariables, Ewoms::ImmiscibleFluxVariables<TypeTag>);

//! The indices required by the isothermal immiscible multi-phase model
SET_TYPE_PROP(ImmiscibleModel, Indices, Ewoms::ImmiscibleIndices<TypeTag, /*PVOffset=*/0>);

//! Disable the energy equation by default
SET_BOOL_PROP(ImmiscibleModel, EnableEnergy, false);

/////////////////////
// set slightly different properties for the single-phase case
/////////////////////

//! The fluid system to use by default
SET_PROP(ImmiscibleOnePhaseModel, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
public:
    typedef Opm::FluidSystems::OneP<Scalar , Fluid> type;
};

SET_PROP(ImmiscibleOnePhaseModel, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

// disable output of a few quantities which make sense in a
// multi-phase but not in a single-phase context
SET_BOOL_PROP(ImmiscibleOnePhaseModel, VtkWriteSaturations, false);
SET_BOOL_PROP(ImmiscibleOnePhaseModel, VtkWriteMobilities, false);
SET_BOOL_PROP(ImmiscibleOnePhaseModel, VtkWriteRelativePermeabilities, false);

/////////////////////
// set slightly different properties for the two-phase case
/////////////////////
SET_PROP(ImmiscibleTwoPhaseModel, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

SET_PROP(ImmiscibleTwoPhaseModel, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

SET_PROP(ImmiscibleTwoPhaseModel, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Opm::FluidSystems::TwoPImmiscible<Scalar, WettingPhase,
                                              NonwettingPhase> type;
};

}} // namespace Properties, Opm

namespace Ewoms {
/*!
 * \ingroup ImmiscibleModel
 * \brief A fully-implicit multi-phase flow model which assumes
 *        immiscibility of the phases.
 *
 * This model implements multi-phase flow of \f$M > 0\f$ immiscible
 * fluids \f$\alpha\f$. By default, the standard multi-phase Darcy
 * approach is used to determine the velocity, i.e.
 * \f[
 * \mathbf{v}_\alpha =
 * - \frac{k_{r\alpha}}{\mu_\alpha}
 * \mathbf{K}\left(\mathbf{grad}\, p_\alpha -
 *                 \varrho_{\alpha} \mathbf{g} \right) \;,
 * \f]
 * although the actual approach which is used can be specified via the
 * \c VelocityModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, VelocityModule, Ewoms::ForchheimerVelocityModule<TypeTag>);
 * \endcode
 *
 * The core of the model is the conservation mass of each component by
 * means of the equation
 * \f[
 * \frac{\partial\;\phi S_\alpha \rho_\alpha }{\partial t}
 * - \mathrm{div} \left\{ \rho_\alpha \mathbf{v}_\alpha  \right\}
 * - q_\alpha = 0 \;.
 * \f]
 *
 * The model uses the following primary variables:
 * - The pressure \f$p_0\f$ in Pascal of the phase with the lowest index
 * - The saturations \f$S_\alpha\f$ of the \f$M - 1\f$ phases that
 *   exhibit the lowest indices
 * - The absolute temperature \f$T\f$ in Kelvin if energy is conserved
 *   via the energy equation
 */
template <class TypeTag>
class ImmiscibleModel
    : public Ewoms::MultiPhaseBaseModel<TypeTag>
{
    typedef Ewoms::MultiPhaseBaseModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { numComponents = FluidSystem::numComponents };

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    enum { dimWorld = GridView::dimensionworld };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    ImmiscibleModel(Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        if (enableEnergy)
            Ewoms::VtkEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "immiscible"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;

        if (pvIdx == Indices::pressure0Idx) {
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        }
        else if (Indices::saturation0Idx <= pvIdx
                 && pvIdx < Indices::saturation0Idx + numPhases - 1) {
            int phaseIdx = pvIdx - Indices::saturation0Idx;
            oss << "saturation_" << FluidSystem::phaseName(phaseIdx);
        }
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(int eqIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::eqName(eqIdx)).empty())
            return s;

        std::ostringstream oss;

        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx
                                                     + numComponents)
            oss << "conti_"
                << FluidSystem::phaseName(eqIdx - Indices::conti0EqIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        // find the a reference pressure. The first degree of freedom
        // might correspond to non-interior entities which would lead
        // to an undefined value, so we have to iterate...
        for (size_t dofIdx = 0; dofIdx < this->numDof(); ++ dofIdx) {
            if (this->dofTotalVolume(dofIdx) > 0) {
                referencePressure_ =
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressure0Idx];
                break;
            }
        }
    }

    /*!
     * \copydetails FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalDofIdx, int pvIdx) const
    {
        assert(referencePressure_ > 0);

        Scalar tmp = EnergyModule::primaryVarWeight(asImp_(), globalDofIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;
        if (Indices::pressure0Idx == pvIdx) {
            return 10 / referencePressure_;
        }
        return 1.0;
    }

    /*!
     * \copydetails FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(int globalDofIdx, int eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(asImp_(), globalDofIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;

#ifndef NDEBUG
        int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);
#endif

        // make all kg equal
        return 1.0;
    }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        if (enableEnergy)
            this->outputModules_.push_back(new Ewoms::VtkEnergyModule<TypeTag>(this->problem_));
    }

private:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    mutable Scalar referencePressure_;
};
} // namespace Ewoms

#endif
