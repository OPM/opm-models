/*
  Copyright (C) 2011-2013 by Andreas Lauser
  Copyright (C) 2012 by Philipp Nuske

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
 * \copydoc Ewoms::NcpModel
 */
#ifndef EWOMS_NCP_MODEL_HH
#define EWOMS_NCP_MODEL_HH

#include "ncpproperties.hh"
#include "ncplocalresidual.hh"
#include "ncpfluxvariables.hh"
#include "ncpprimaryvariables.hh"
#include "ncpboundaryratevector.hh"
#include "ncpratevector.hh"
#include "ncpvolumevariables.hh"
#include "ncpnewtonmethod.hh"
#include "ncpindices.hh"

#include <ewoms/models/common/multiphasebasemodel.hh>
#include <ewoms/models/common/energymodule.hh>
#include <ewoms/models/common/diffusionmodule.hh>
#include <ewoms/vtk/vtkcompositionmodule.hh>
#include <ewoms/vtk/vtkenergymodule.hh>
#include <ewoms/vtk/vtkdiffusionmodule.hh>
#include <opm/material/constraintsolvers/CompositionFromFugacities.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/unused.hh>

#include <sstream>
#include <string>
#include <vector>
#include <array>

namespace Ewoms {
template <class TypeTag>
class NcpModel;
}

namespace Opm {
namespace Properties {
/*!
 * \brief Define the type tag for the compositional NCP model.
 */
NEW_TYPE_TAG(NcpModel, INHERITS_FROM(MultiPhaseBaseModel,
                                     VtkComposition,
                                     VtkEnergy,
                                     VtkDiffusion));

/*!
 * \brief Set the themodynamic constraint solver which calculates the
 *        composition of any phase given all component fugacities.
 */
SET_PROP(NcpModel, NcpCompositionFromFugacitiesSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Opm::CompositionFromFugacities<Scalar, FluidSystem> type;
};

//! Use the Ncp local jacobian operator for the compositional NCP model
SET_TYPE_PROP(NcpModel,
              LocalResidual,
              Ewoms::NcpLocalResidual<TypeTag>);

//! Use the Ncp specific newton method for the compositional NCP model
SET_TYPE_PROP(NcpModel, NewtonMethod, Ewoms::NcpNewtonMethod<TypeTag>);

//! the Model property
SET_TYPE_PROP(NcpModel, Model, Ewoms::NcpModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(NcpModel, BaseProblem, Ewoms::MultiPhaseBaseProblem<TypeTag>);

//! Disable the energy equation by default
SET_BOOL_PROP(NcpModel, EnableEnergy, false);

//! disable diffusion by default
SET_BOOL_PROP(NcpModel, EnableDiffusion, false);

//! the RateVector property
SET_TYPE_PROP(NcpModel, RateVector, Ewoms::NcpRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(NcpModel, BoundaryRateVector, Ewoms::NcpBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(NcpModel, PrimaryVariables, Ewoms::NcpPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(NcpModel, VolumeVariables, Ewoms::NcpVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(NcpModel, FluxVariables, Ewoms::NcpFluxVariables<TypeTag>);

//! truncate the newton update for the first 3 iterations of a time step
SET_INT_PROP(NcpModel, NcpNewtonNumChoppedIterations, 3);

//! The indices required by the compositional NCP model
SET_TYPE_PROP(NcpModel, Indices, Ewoms::NcpIndices<TypeTag, 0>);

//! The unmodified weight for the pressure primary variable
SET_SCALAR_PROP(NcpModel, NcpPressureBaseWeight, 1.0);
//! The weight for the saturation primary variables
SET_SCALAR_PROP(NcpModel, NcpSaturationsBaseWeight, 1.0);
//! The unmodified weight for the fugacity primary variables
SET_SCALAR_PROP(NcpModel, NcpFugacitiesBaseWeight, 1.0);

}} // namespace Properties, Opm

namespace Ewoms {

/*!
 * \ingroup NcpModel
 *
 * \brief A compositional multi-phase model based on non-linear
 *        complementarity functions.
 *
 * This model implements a \f$M\f$-phase flow of a fluid mixture
 * composed of \f$N\f$ chemical species. The phases are denoted by
 * lower index \f$\alpha \in \{ 1, \dots, M \}\f$. All fluid phases
 * are mixtures of \f$N \geq M - 1\f$ chemical species which are
 * denoted by the upper index \f$\kappa \in \{ 1, \dots, N \} \f$.
 *
 *
 * By default, the standard multi-phase Darcy approach is used to determine
 * the velocity, i.e.
 * \f[
 *   \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 *   \left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;,
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
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha \right\}
 * - q^\kappa = 0 \;.
 * \f]
 *
 * For the missing \f$M\f$ model assumptions, the model uses
 * non-linear complementarity functions. These are based on the
 * observation that if a fluid phase is not present, the sum of the
 * mole fractions of this fluid phase is smaller than \f$1\f$, i.e.
 * \f[ \forall \alpha: S_\alpha = 0 \implies \sum_\kappa
 * x_\alpha^\kappa \leq 1 \f]
 *
 * Also, if a fluid phase may be present at a given spatial location
 * its saturation must be non-negative:
 * \f[ \forall \alpha: \sum_\kappa x_\alpha^\kappa = 1 \implies S_\alpha \geq 0
 *\f]
 *
 * Since at any given spatial location, a phase is always either
 * present or not present, one of the strict equalities on the
 * right hand side is always true, i.e.
 * \f[
 * \forall \alpha: S_\alpha \left( \sum_\kappa x_\alpha^\kappa - 1 \right) = 0
 * \f]
 * always holds.
 *
 * These three equations constitute a non-linear complementarity
 * problem, which can be solved using so-called non-linear
 * complementarity functions \f$\Phi(a, b)\f$. Such functions have the property
 * \f[\Phi(a,b) = 0 \iff a \geq0 \land b \geq0  \land a \cdot b = 0 \f]
 *
 * Several non-linear complementarity functions have been suggested,
 * e.g. the Fischer-Burmeister function
 * \f[ \Phi(a,b) = a + b - \sqrt{a^2 + b^2} \;. \f]
 * This model uses
 * \f[ \Phi(a,b) = \min \{a,  b \}\;, \f]
 * because of its piecewise linearity.
 *
 * These equations are then discretized using a fully-implicit vertex
 * centered finite volume scheme for spatial discretization and the
 * implicit Euler method as temporal discretization.
 *
 * The model assumes local thermodynamic equilibrium and uses the
 * following primary variables:
 * - The pressure of the first phase \f$p_1\f$
 * - The component fugacities \f$f^1, \dots, f^{N}\f$
 * - The saturations of the first \f$M-1\f$ phases \f$S_1, \dots, S_{M-1}\f$
 * - Temperature \f$T\f$ if the energy equation is enabled
 */
template <class TypeTag>
class NcpModel
    : public MultiPhaseBaseModel<TypeTag>
{
    typedef MultiPhaseBaseModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { ncp0EqIdx = Indices::ncp0EqIdx };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { dimWorld = GridView::dimensionworld };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    typedef Ewoms::EnergyModule<TypeTag, enableEnergy> EnergyModule;
    typedef Ewoms::DiffusionModule<TypeTag, enableDiffusion> DiffusionModule;

public:
    NcpModel(Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        DiffusionModule::registerParameters();
        EnergyModule::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkCompositionModule<TypeTag>::registerParameters();

        if (enableDiffusion)
            Ewoms::VtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Ewoms::VtkEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc FvBaseDiscretization::init
     */
    void init()
    {
        ParentType::init();
        minActivityCoeff_.resize(this->numDof());
        std::fill(minActivityCoeff_.begin(), minActivityCoeff_.end(), 1.0);
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
    { return "ncp"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;
        if (pvIdx == pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (saturation0Idx <= pvIdx && pvIdx < saturation0Idx
                                                    + (numPhases - 1))
            oss << "saturation_"
                << FluidSystem::phaseName(/*phaseIdx=*/pvIdx - saturation0Idx);
        else if (fugacity0Idx <= pvIdx && pvIdx < fugacity0Idx + numComponents)
            oss << "fugacity^"
                << FluidSystem::componentName(pvIdx - fugacity0Idx);
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
        if (conti0EqIdx <= eqIdx && eqIdx < conti0EqIdx + numComponents)
            oss << "continuity^"
                << FluidSystem::componentName(eqIdx - conti0EqIdx);
        else if (ncp0EqIdx <= eqIdx && eqIdx < ncp0EqIdx + numPhases)
            oss << "ncp_"
                << FluidSystem::phaseName(/*phaseIdx=*/eqIdx - ncp0EqIdx);
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
     * \copydoc FvBaseDiscretization::updatePVWeights
     */
    void updatePVWeights(const ElementContext &elemCtx) const
    {
        for (int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
            int globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                minActivityCoeff_[globalIdx][compIdx] = 1e100;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    const auto &fs = elemCtx.volVars(dofIdx, /*timeIdx=*/0).fluidState();

                    minActivityCoeff_[globalIdx][compIdx]
                        = std::min(minActivityCoeff_[globalIdx][compIdx],
                                   fs.fugacityCoefficient(phaseIdx, compIdx)
                                   * fs.pressure(phaseIdx));
                    Valgrind::CheckDefined(minActivityCoeff_[globalIdx][compIdx]);
                };
            };
        };
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalDofIdx, int pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalDofIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;
        else if (fugacity0Idx <= pvIdx && pvIdx < fugacity0Idx + numComponents) {
            // component fugacity
            int compIdx = pvIdx - fugacity0Idx;
            assert(0 <= compIdx && compIdx <= numComponents);

            Valgrind::CheckDefined(minActivityCoeff_[globalDofIdx][compIdx]);
            static const Scalar fugacityBaseWeight =
                GET_PROP_VALUE(TypeTag, NcpFugacitiesBaseWeight);
            return fugacityBaseWeight / minActivityCoeff_[globalDofIdx][compIdx];
        }
        else if (Indices::pressure0Idx == pvIdx) {
            static const Scalar pressureBaseWeight
                = GET_PROP_VALUE(TypeTag, NcpPressureBaseWeight);
            return pressureBaseWeight / referencePressure_;
        }

#ifndef NDEBUG
        int phaseIdx = pvIdx - saturation0Idx;
        assert(0 <= phaseIdx && phaseIdx < numPhases - 1);
#endif

        // saturation
        static const Scalar saturationsBaseWeight
            = GET_PROP_VALUE(TypeTag, NcpSaturationsBaseWeight);
        return saturationsBaseWeight;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(int globalDofIdx, int eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalDofIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;
        // an NCP
        else if (ncp0EqIdx <= eqIdx && eqIdx < Indices::ncp0EqIdx + numPhases)
            return 1.0;

        // mass conservation equation
        int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numComponents);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    /*!
     * \brief Returns the smallest activity coefficient of a component for the
     *        most current solution at a vertex.
     *
     * \param globalDofIdx The global index of the vertex (i.e. finite volume) of interest.
     * \param compIdx The index of the component of interest.
     */
    Scalar minActivityCoeff(int globalDofIdx, int compIdx) const
    { return minActivityCoeff_[globalDofIdx][compIdx]; }

    /*!
     * \internal
     */
    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        this->vtkOutputModules_.push_back(
            new Ewoms::VtkCompositionModule<TypeTag>(this->problem_));
        if (enableDiffusion)
            this->vtkOutputModules_.push_back(
                new Ewoms::VtkDiffusionModule<TypeTag>(this->problem_));
        if (enableEnergy)
            this->vtkOutputModules_.push_back(
                new Ewoms::VtkEnergyModule<TypeTag>(this->problem_));
    }

    mutable Scalar referencePressure_;
    mutable std::vector<ComponentVector> minActivityCoeff_;
};

} // namespace Ewoms

#endif
