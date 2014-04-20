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
 * \copydoc Ewoms::BlackOilModel
 */
#ifndef EWOMS_BLACK_OIL_MODEL_HH
#define EWOMS_BLACK_OIL_MODEL_HH

#include "blackoilproblem.hh"
#include "blackoilindices.hh"
#include "blackoilfluxvariables.hh"
#include "blackoilprimaryvariables.hh"
#include "blackoilvolumevariables.hh"
#include "blackoilratevector.hh"
#include "blackoilboundaryratevector.hh"
#include "blackoillocalresidual.hh"
#include "blackoilproperties.hh"

#include <ewoms/models/common/multiphasebasemodel.hh>
#include <ewoms/io/vtkcompositionmodule.hh>
#include <ewoms/io/vtkblackoilmodule.hh>

#include <ewoms/io/eclipseoutputblackoilmodule.hh>
#include <ewoms/io/eclipsewriter.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class BlackOilModel;
}

namespace Opm {
namespace Properties {
//! The type tag for the black-oil problems
NEW_TYPE_TAG(BlackOilModel, INHERITS_FROM(MultiPhaseBaseModel, VtkBlackOil, VtkComposition));

//! Set the local residual function
SET_TYPE_PROP(BlackOilModel, LocalResidual,
              Ewoms::BlackOilLocalResidual<TypeTag>);

//! The Model property
SET_TYPE_PROP(BlackOilModel, Model, Ewoms::BlackOilModel<TypeTag>);

//! The Problem property
SET_TYPE_PROP(BlackOilModel, BaseProblem, Ewoms::BlackOilProblem<TypeTag>);

//! The BlackOilFluidState property
SET_PROP(BlackOilModel, BlackOilFluidState)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef Opm::CompositionalFluidState<Scalar,
                                         FluidSystem,
                                         /*enableEnthalpy=*/false> type;
};

//! the RateVector property
SET_TYPE_PROP(BlackOilModel, RateVector, Ewoms::BlackOilRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(BlackOilModel, BoundaryRateVector, Ewoms::BlackOilBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(BlackOilModel, PrimaryVariables, Ewoms::BlackOilPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BlackOilModel, VolumeVariables, Ewoms::BlackOilVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BlackOilModel, FluxVariables, Ewoms::BlackOilFluxVariables<TypeTag>);

//! The indices required by the model
SET_TYPE_PROP(BlackOilModel, Indices, Ewoms::BlackOilIndices</*PVOffset=*/0>);

//! Set the fluid system to the black-oil fluid system by default
SET_TYPE_PROP(BlackOilModel, FluidSystem,
              Opm::FluidSystems::BlackOil<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Only produce Eclipse output if the necessary preconditions are fulfilled.
SET_PROP(BlackOilModel, EnableEclipseOutput)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

#if HAVE_DUNE_CORNERPOINT
    static const bool isCpGrid = std::is_same<Grid, Dune::CpGrid>::value;
#else
    static const bool isCpGrid = false;
#endif

public:
    // in addition to Dune::CpGrid, the ERT libraries must be
    // available...
    static const bool value = isCpGrid && HAVE_ERT;
};

}} // namespace Properties, Opm

namespace Ewoms {
/*!
 * \ingroup BlackOilModel
 * \brief A fully-implicit black-oil flow model.
 *
 * The black-oil model is a three-phase, three-component model widely
 * used for oil reservoir simulation.  The phases are denoted by lower
 * index \f$\alpha \in \{ w, g, o \}\f$ ("water", "gas" and "oil") and
 * the components by upper index \f$\kappa \in \{ W, G, O \}\f$
 * ("Water", "Gas" and "Oil"). The model assumes partial miscibility:
 *
 * - Water and the gas phases are immisicible and are assumed to be
 *   only composed of the water and gas components respectively-
 * - The oil phase is assumed to be a mixture of the gas and the oil
 *  components.
 *
 * The densities of the phases are determined by so-called
 * <i>formation volume factors</i>:
 *
 * \f[
 * B_\alpha := \frac{\varrho_\alpha(1\,\text{bar})}{\varrho_\alpha(p_\alpha)}
 * \f]
 *
 * Since the gas and water phases are assumed to be immiscible, this
 * is sufficint to calculate their density. For the formation volume
 * factor of the the oil phase \f$B_o\f$ determines the density of
 * *saturated* oil, i.e. the density of the oil phase if some gas
 * phase is present.
 *
 * The composition of the oil phase is given by the <i>gas formation
 * factor</i> \f$R_s\f$, which defined as the volume of gas at
 * atmospheric pressure that is dissolved in saturated oil at a given
 * pressure:
 *
 * \f[
 * R_s := \frac{x_o^G(p)\,\varrho_{mol,o}(p)}{\varrho_g(1\,\text{bar})}\;.
 * \f]
 *
 * This allows to calculate all quantities required for the
 * mass-conservation equations for each component, i.e.
 *
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha \right\}
 * - q^\kappa = 0 \;,
 * \f]
 * where \f$\mathrm{v}_\alpha\f$ is the filter velocity of the phase
 * \f$\alpha\f$.
 *
 * By default \f$\mathrm{v}_\alpha\f$ is determined by using the
 * standard multi-phase Darcy approach, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 *\left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
 * although the actual approach which is used can be specified via the
 * \c VelocityModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, VelocityModule, Ewoms::ForchheimerVelocityModule<TypeTag>);
 * \endcode
 *
 * The primary variables used by this model are:
 * - The pressure of the phase with the lowest index
 * - The two saturations of the phases with the lowest indices
 */
template<class TypeTag >
class BlackOilModel
    : public MultiPhaseBaseModel<TypeTag>
{
    typedef MultiPhaseBaseModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = FluidSystem::numComponents };

public:
    BlackOilModel(Simulator &simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkBlackOilModule<TypeTag>::registerParameters();
        Ewoms::VtkCompositionModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::init
     */
    void init()
    {
        ParentType::init();

        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "blackoil"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::ostringstream oss;

        if (pvIdx == Indices::pressure0Idx) {
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        }
        else if (Indices::saturation0Idx <= pvIdx
                 && pvIdx <= Indices::saturation0Idx + numPhases - 1)
        {
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
        std::ostringstream oss;

        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "conti_" << FluidSystem::phaseName(eqIdx - Indices::conti0EqIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalDofIdx, int pvIdx) const
    {
        if (Indices::pressure0Idx == pvIdx) {
            return 10/referencePressure_;
        }

        return 1;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(int globalDofIdx, int eqIdx) const
    {
        int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
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

// HACK: this should be made private and the BaseModel should be
// declared to be a friend. Since C++-2003 (and more relevantly GCC
// 4.4) don't support friend typedefs, we need to make this method
// public until the oldest supported compiler supports friend
// typedefs...

//protected:
//    friend typename GET_PROP_TYPE(TypeTag, Discretization);

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which make sense for the blackoil model
        this->outputModules_.push_back(new Ewoms::VtkBlackOilModule<TypeTag>(this->simulator_));
        this->outputModules_.push_back(new Ewoms::VtkCompositionModule<TypeTag>(this->simulator_));

        if (enableEclipseOutput_()) {
            // add the output module for the Eclipse binary output
            this->outputModules_.push_back(new Ewoms::EclipseOutputBlackOilModule<TypeTag>(this->simulator_));
        }
    }

private:
    static bool enableEclipseOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclipseOutput); }

    mutable Scalar referencePressure_;
};
} // namespace Ewoms

#endif
