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
#include "blackoilextensivequantities.hh"
#include "blackoilprimaryvariables.hh"
#include "blackoilintensivequantities.hh"
#include "blackoilratevector.hh"
#include "blackoilboundaryratevector.hh"
#include "blackoillocalresidual.hh"
#include "blackoilnewtonmethod.hh"
#include "blackoilproperties.hh"

#include <ewoms/models/common/multiphasebasemodel.hh>
#include <ewoms/io/vtkcompositionmodule.hh>
#include <ewoms/io/vtkblackoilmodule.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class BlackOilModel;

template <class TypeTag>
class EclGridManager;
}

namespace Opm {
namespace Properties {
//! The type tag for the black-oil problems
NEW_TYPE_TAG(BlackOilModel, INHERITS_FROM(MultiPhaseBaseModel,
                                          VtkBlackOil,
                                          VtkComposition));

//! Set the local residual function
SET_TYPE_PROP(BlackOilModel, LocalResidual,
              Ewoms::BlackOilLocalResidual<TypeTag>);

//! Use the black-oil specific newton method
SET_TYPE_PROP(BlackOilModel, NewtonMethod, Ewoms::BlackOilNewtonMethod<TypeTag>);

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

//! the IntensiveQuantities property
SET_TYPE_PROP(BlackOilModel, IntensiveQuantities, Ewoms::BlackOilIntensiveQuantities<TypeTag>);

//! the ExtensiveQuantities property
SET_TYPE_PROP(BlackOilModel, ExtensiveQuantities, Ewoms::BlackOilExtensiveQuantities<TypeTag>);

//! The indices required by the model
SET_TYPE_PROP(BlackOilModel, Indices, Ewoms::BlackOilIndices</*PVOffset=*/0>);

//! Set the fluid system to the black-oil fluid system by default
SET_TYPE_PROP(BlackOilModel, FluidSystem,
              Opm::FluidSystems::BlackOil<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Set the number of Newton-Raphson iterations for which the update should be chopped to
//! 4 by default
SET_INT_PROP(BlackOilModel, BlackoilNumChoppedIterations, 4);
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
 * The composition of the oil phase is given by the <i>gas dissolution factor</i>
 * \f$R_s\f$, which defined as the volume of gas at atmospheric pressure that is
 * dissolved in a given amount of oil at reservoir pressure:
 *
 * \f[
 * R_s := \frac{\varrho_{o}^G}{\varrho_o^O}\;.
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
 * \c FluxModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, FluxModule, Ewoms::ForchheimerFluxModule<TypeTag>);
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
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef MultiPhaseBaseModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = FluidSystem::numComponents };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

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
     * \copydoc FvBaseDiscretization::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

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

        if (pvIdx == Indices::gasPressureIdx)
            oss << "pressure_" << FluidSystem::phaseName(FluidSystem::gasPhaseIdx);
        else if (pvIdx == Indices::waterSaturationIdx)
            oss << "saturation_" << FluidSystem::phaseName(FluidSystem::waterPhaseIdx);
        else if (pvIdx == Indices::switchIdx)
            oss << "switching";
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
        // do not care about the auxiliary equations as they are supposed to scale
        // themselves
        if (globalDofIdx >= (int) this->numGridDof())
            return 1;

        if (Indices::gasPressureIdx == pvIdx)
            return 10/referencePressure_;

        return 1;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(int globalDofIdx, int eqIdx) const
    {
        // do not care about the auxiliary equations as they are supposed to scale
        // themselves
        if (globalDofIdx >= (int) this->numGridDof())
            return 1;

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
        for (size_t dofIdx = 0; dofIdx < this->numGridDof(); ++ dofIdx) {
            if (this->isLocalDof(dofIdx)) {
                referencePressure_ = this->solution(/*timeIdx=*/0)[dofIdx][Indices::gasPressureIdx];
                break;
            }
        }
    }

    /*!
     * \internal
     * \brief Do the primary variable switching after a Newton iteration.
     *
     * This is an internal method that needs to be public because it
     * gets called by the Newton method after an update.
     */
    void switchPrimaryVars_()
    {
        numSwitched_ = 0;

        int numDof = this->numGridDof();
        for (int globalDofIdx = 0; globalDofIdx < numDof; ++globalDofIdx) {
            auto &priVars = this->solution(/*timeIdx=*/0)[globalDofIdx];
            if (priVars.adaptSwitchingVariable()) {
                this->linearizer().markDofRed(globalDofIdx);
                ++numSwitched_;
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        numSwitched_ = this->gridView_.comm().sum(numSwitched_);

        this->simulator_.model().newtonMethod().endIterMsg()
            << ", num switched=" << numSwitched_;
    }

    /*!
     * \brief Write the current solution for a degree of freedom to a
     *        restart file.
     *
     * \param outstream The stream into which the vertex data should
     *                  be serialized to
     * \param dof The Dune entity which's data should be serialized
     */
    template <class DofEntity>
    void serializeEntity(std::ostream &outstream, const DofEntity &dof)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int dofIdx = asImp_().dofMapper().index(dof);
#else
        int dofIdx = asImp_().dofMapper().map(dof);
#endif

        // write phase state
        if (!outstream.good()) {
            OPM_THROW(std::runtime_error, "Could not serialize degree of freedom " << dofIdx);
        }

        // write the primary variables
        const auto& priVars = this->solution_[/*timeIdx=*/0][dofIdx];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            outstream << priVars[eqIdx] << " ";

        // write the pseudo primary variables
        outstream << priVars.switchingVarMeaning() << " ";
        outstream << priVars.pvtRegionIndex() << " ";
    }

    /*!
     * \brief Reads the current solution variables for a degree of
     *        freedom from a restart file.
     *
     * \param instream The stream from which the vertex data should
     *                  be deserialized from
     * \param dof The Dune entity which's data should be deserialized
     */
    template <class DofEntity>
    void deserializeEntity(std::istream &instream,
                           const DofEntity &dof)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int dofIdx = asImp_().dofMapper().index(dof);
#else
        int dofIdx = asImp_().dofMapper().map(dof);
#endif

        // read in the "real" primary variables of the DOF
        auto& priVars = this->solution_[/*timeIdx=*/0][dofIdx];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            if (!instream.good())
                OPM_THROW(std::runtime_error,
                          "Could not deserialize degree of freedom " << dofIdx);
            instream >> priVars[eqIdx];
        }

        // read the pseudo primary variables
        int switchingVarMeaning;
        instream >> switchingVarMeaning;

        int pvtRegionIdx;
        instream >> pvtRegionIdx;

        if (!instream.good())
            OPM_THROW(std::runtime_error,
                      "Could not deserialize degree of freedom " << dofIdx);

        typedef typename PrimaryVariables::SwitchingVarMeaning SVM;
        priVars.setSwitchingVarMeaning(static_cast<SVM>(switchingVarMeaning));
        priVars.setPvtRegionIndex(pvtRegionIdx);
    }

    /*!
     * \brief Deserializes the state of the model.
     *
     * \tparam Restarter The type of the serializer class
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        ParentType::deserialize(res);

        // set the PVT indices of the primary variables. This is also done by writing
        // them into the restart file and re-reading them, but it is better to calculate
        // them from scratch because the input could have been changed in this regard...
        ElementContext elemCtx(this->simulator_);
        auto elemIt = this->gridView().template begin</*codim=*/0>();
        auto elemEndIt = this->gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            elemCtx.updateStencil(*elemIt);
            for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timIdx=*/0); ++dofIdx) {
                int globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timIdx=*/0);
                updatePvtRegionIndex_(this->solution_[/*timeIdx=*/0][globalDofIdx],
                                      elemCtx, dofIdx, /*timeIdx=*/0);
            }
        }

        this->solution_[/*timeIdx=*/1] = this->solution_[/*timeIdx=*/0];
    }


// HACK: this should be made private and the BaseModel should be
// declared to be a friend. Since C++-2003 (and more relevantly GCC
// 4.4) don't support friend typedefs, we need to make this method
// public until the oldest supported compiler supports friend
// typedefs...

//protected:
//    friend typename GET_PROP_TYPE(TypeTag, Discretization);

    template <class Context>
    void supplementInitialSolution_(PrimaryVariables &priVars,
                                    const Context &context, int dofIdx, int timeIdx)
    { updatePvtRegionIndex_(priVars, context, dofIdx, timeIdx); }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();

        // add the VTK output modules which make sense for the blackoil model
        this->addOutputModule(new Ewoms::VtkBlackOilModule<TypeTag>(this->simulator_));
        this->addOutputModule(new Ewoms::VtkCompositionModule<TypeTag>(this->simulator_));
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    template <class Context>
    void updatePvtRegionIndex_(PrimaryVariables &priVars, const Context &context, int dofIdx, int timeIdx)
    {
        int regionIdx = context.problem().pvtRegionIndex(context, dofIdx, timeIdx);
        priVars.setPvtRegionIndex(regionIdx);
    }

    mutable Scalar referencePressure_;
    int numSwitched_;
};
} // namespace Ewoms

#endif
