// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Contains the classes required to extend the black-oil model by solvents.
 */
#ifndef EWOMS_BLACK_OIL_SOLVENT_MODULE_HH
#define EWOMS_BLACK_OIL_SOLVENT_MODULE_HH

#include "blackoilproperties.hh"
#include <ewoms/io/vtkblackoilsolventmodule.hh>
#include <ewoms/models/common/quantitycallbacks.hh>

#include <opm/material/fluidsystems/blackoilpvt/SolventPvt.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SsfnTable.hpp>
#endif

#include <opm/common/Valgrind.hpp>
#include <opm/common/Unused.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Ewoms {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by solvents.
 */
template <class TypeTag, bool enableSolventV = GET_PROP_VALUE(TypeTag, EnableSolvent)>
class BlackOilSolventModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Opm::SolventPvt<Scalar> SolventPvt;

    typedef typename Opm::Tabulated1DFunction<Scalar> TabulatedFunction;

    static constexpr unsigned solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr unsigned contiSolventEqIdx = Indices::contiSolventEqIdx;
    static constexpr unsigned enableSolvent = enableSolventV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize all internal data structures needed by the solvent module
     */
    static void initFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        // some sanity checks: if solvents are enabled, the SOLVENT keyword must be
        // present, if solvents are disabled the keyword must not be present.
        if (enableSolvent && !deck.hasKeyword("SOLVENT")) {
            OPM_THROW(std::runtime_error,
                      "Non-trivial solvent treatment requested at compile time, but "
                      "the deck does not contain the SOLVENT keyword");
        }
        else if (!enableSolvent && deck.hasKeyword("SOLVENT")) {
            OPM_THROW(std::runtime_error,
                      "Solvent treatment disabled at compile time, but the deck "
                      "contains the SOLVENT keyword");
        }

        if (!deck.hasKeyword("SOLVENT"))
            return; // solvent treatment is supposed to be disabled

        solventPvt_.initFromDeck(deck, eclState);

        // initialize the objects which deal with the SSFN keyword
        const auto& tableManager = eclState.getTableManager();
        const auto& ssfnTables = tableManager.getSsfnTables();
        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        setNumSatRegions(numSatRegions);
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++ satRegionIdx) {
            const auto& ssfnTable = ssfnTables.template getTable<Opm::SsfnTable>(satRegionIdx);
            ssfnKrg_[satRegionIdx].setXYContainers(ssfnTable.getSolventFractionColumn(),
                                                   ssfnTable.getGasRelPermMultiplierColumn(),
                                                   /*sortInput=*/true);
            ssfnKrs_[satRegionIdx].setXYContainers(ssfnTable.getSolventFractionColumn(),
                                                   ssfnTable.getSolventRelPermMultiplierColumn(),
                                                   /*sortInput=*/true);
        }
    }
#endif

    /*!
     * \brief Specify the number of satuation regions.
     *
     * This must be called before setting the SSFN of any region.
     */
    static void setNumSatRegions(unsigned numRegions)
    {
        ssfnKrg_.resize(numRegions);
        ssfnKrs_.resize(numRegions);
    }

    /*!
     * \brief Specify the solvent saturation functions of a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    static void setSsfn(unsigned satRegionIdx,
                        const TabulatedFunction& ssfnKrg,
                        const TabulatedFunction& ssfnKrs)
    {
        ssfnKrg_[satRegionIdx] = ssfnKrg;
        ssfnKrs_[satRegionIdx] = ssfnKrs;
    }

    /*!
     * \brief Specify the solvent PVT of a all PVT regions.
     */
    static void setSolventPvt(const SolventPvt& value)
    { solventPvt_ = value; }

    /*!
     * \brief Register all run-time parameters for the black-oil solvent module.
     */
    static void registerParameters()
    {
        if (!enableSolvent)
            // solvents have disabled at compile time
            return;

        Ewoms::VtkBlackOilSolventModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all solvent specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableSolvent)
            // solvents have disabled at compile time
            return;

        model.addOutputModule(new Ewoms::VtkBlackOilSolventModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enableSolvent)
            // solvents have disabled at compile time
            return false;

        return pvIdx == solventSaturationIdx;
    }

    static std::string primaryVarName(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        return "saturation_solvent";
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableSolvent)
            return false;

        return eqIdx == contiSolventEqIdx;
    }

    static std::string eqName(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return "conti^solvent";
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableSolvent)
            return;

        storage[contiSolventEqIdx] +=
            Toolbox::template decay<LhsEval>(intQuants.porosity())
            * Toolbox::template decay<LhsEval>(intQuants.solventSaturation())
            * Toolbox::template decay<LhsEval>(intQuants.solventDensity());
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enableSolvent)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        unsigned upIdx = extQuants.solventUpstreamIndex();
        unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);

        if (upIdx == inIdx)
            flux[contiSolventEqIdx] =
                extQuants.solventVolumeFlux()
                *up.solventDensity();
        else
            flux[contiSolventEqIdx] =
                extQuants.solventVolumeFlux()
                *Opm::decay<Scalar>(up.solventDensity());
    }

    /*!
     * \brief Assign the solvent specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar solventSaturation)
    {
        if (!enableSolvent)
            return;

        priVars[solventSaturationIdx] = solventSaturation;
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the solvents.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if (!enableSolvent)
            return;

        // do a plain unchopped Newton update
        newPv[solventSaturationIdx] = oldPv[solventSaturationIdx] - delta[solventSaturationIdx];
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the cange of solvent primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid)
    {
        // do not weight the residual of solvents when it comes to convergence
        return std::abs(Toolbox::scalarValue(resid[contiSolventEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableSolvent)
            return;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
        unsigned dofIdx = model.dofMapper().index(dof);
#else
        unsigned dofIdx = model.dofMapper().map(dof);
#endif

        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[solventSaturationIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableSolvent)
            return;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,4)
        unsigned dofIdx = model.dofMapper().index(dof);
#else
        unsigned dofIdx = model.dofMapper().map(dof);
#endif

        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[solventSaturationIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1 = priVars0[solventSaturationIdx];
    }

    static const SolventPvt& solventPvt()
    { return solventPvt_; }

    static const TabulatedFunction& ssfnKrg(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return ssfnKrg_[satnumRegionIdx];
    }

    static const TabulatedFunction& ssfnKrs(const ElementContext& elemCtx,
                                            unsigned scvIdx,
                                            unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return ssfnKrs_[satnumRegionIdx];
    }

private:
    static SolventPvt solventPvt_;

    static std::vector<TabulatedFunction> ssfnKrg_; // the krg(Fs) column of the SSFN table
    static std::vector<TabulatedFunction> ssfnKrs_; // the krs(Fs) column of the SSFN table
};

template <class TypeTag, bool enableSolventV>
typename BlackOilSolventModule<TypeTag, enableSolventV>::SolventPvt
BlackOilSolventModule<TypeTag, enableSolventV>::solventPvt_;

template <class TypeTag, bool enableSolventV>
std::vector<typename BlackOilSolventModule<TypeTag, enableSolventV>::TabulatedFunction>
BlackOilSolventModule<TypeTag, enableSolventV>::ssfnKrg_;

template <class TypeTag, bool enableSolventV>
std::vector<typename BlackOilSolventModule<TypeTag, enableSolventV>::TabulatedFunction>
BlackOilSolventModule<TypeTag, enableSolventV>::ssfnKrs_;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilSolventIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        solvents extension of the black-oil model.
 */
template <class TypeTag, bool enableSolventV = GET_PROP_VALUE(TypeTag, EnableSolvent)>
class BlackOilSolventIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilSolventModule<TypeTag> SolventModule;

    static constexpr int solventSaturationIdx = Indices::solventSaturationIdx;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;

public:
    /*!
     * \brief Called before the saturation functions are doing their magic
     *
     * At this point, the saturations of the fluid state correspond to those if the phases
     * were pure hydrocarbons.
     */
    void solventPreSatFuncUpdate_(const ElementContext& elemCtx,
                                  unsigned scvIdx,
                                  unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(scvIdx, timeIdx);
        auto& fs = asImp_().fluidState_;
        solventSaturation_ = priVars.makeEvaluation(solventSaturationIdx, timeIdx);
        hydrocarbonSaturation_ = fs.saturation(gasPhaseIdx);

        // make the saturation of the gas phase which is used by the saturation functions
        // the sum of the solvent "saturation" and the saturation the hydrocarbon gas.
        fs.setSaturation(gasPhaseIdx, hydrocarbonSaturation_ + solventSaturation_);
    }

    /*!
     * \brief Called after the saturation functions have been doing their magic
     *
     * At this point, the pressures of the fluid state are final. After this function,
     * all saturations and relative permeabilities must be final. (i.e., the "hydrocarbon
     * saturations".)
     */
    void solventPostSatFuncUpdate_(const ElementContext& elemCtx,
                                   unsigned scvIdx,
                                   unsigned timeIdx)
    {
        // revert the gas "saturation" of the fluid state back to the saturation of the
        // hydrocarbon gas.
        auto& fs = asImp_().fluidState_;
        fs.setSaturation(gasPhaseIdx, hydrocarbonSaturation_);

        const auto& ssfnKrg = SolventModule::ssfnKrg(elemCtx, scvIdx, timeIdx);
        const auto& ssfnKrs = SolventModule::ssfnKrs(elemCtx, scvIdx, timeIdx);

        // compute the mobility of the solvent "phase". this only covers the "immiscible"
        // case.
        solventMobility_ = 0.0;
        Evaluation Stot = hydrocarbonSaturation_ + solventSaturation_;
        if (Stot > 1e-12) { // apply a cut-off
            Evaluation Fhydgas = hydrocarbonSaturation_/Stot;
            Evaluation Fsolgas = solventSaturation_/Stot;

            Evaluation& krg = asImp_().mobility_[gasPhaseIdx];
            solventMobility_ = krg * ssfnKrs.eval(Fsolgas, /*extrapolate=*/true);
            krg *= ssfnKrg.eval(Fhydgas, /*extrapolate=*/true);
        }
    }

    /*!
     * \brief Update the intensive PVT properties needed to handle solvents from the
     *        primary variables.
     *
     * At this point the pressures and saturations of the fluid state are correct.
     */
    void solventPvtUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                           unsigned scvIdx OPM_UNUSED,
                           unsigned timeIdx OPM_UNUSED)
    {
        const auto& iq = asImp_();
        const auto& fs = iq.fluidState();
        const auto& solventPvt = SolventModule::solventPvt();

        unsigned pvtRegionIdx = iq.pvtRegionIndex();
        solventRefDensity_ = solventPvt.referenceDensity(pvtRegionIdx);
        const Evaluation& T = fs.temperature(gasPhaseIdx);
        const Evaluation& p = fs.pressure(gasPhaseIdx);
        solventInvFormationVolumeFactor_ = solventPvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p);

        // TODO: implement solvent "miscibility"
        solventDensity_ = solventInvFormationVolumeFactor_*solventRefDensity_;
        solventViscosity_ = solventPvt.viscosity(pvtRegionIdx, T, p);
        solventMobility_ /= solventViscosity_;
    }

    const Evaluation& solventSaturation() const
    { return solventSaturation_; }

    const Evaluation& solventDensity() const
    { return solventDensity_; }

    const Evaluation& solventViscosity() const
    { return solventViscosity_; }

    const Evaluation& solventMobility() const
    { return solventMobility_; }

    const Evaluation& solventInverseFormationVolumeFactor() const
    { return solventInvFormationVolumeFactor_; }

    // This could be stored pr pvtRegion instead
    const Scalar& solventRefDensity() const
    { return solventRefDensity_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation hydrocarbonSaturation_;
    Evaluation solventSaturation_;
    Evaluation solventDensity_;
    Evaluation solventViscosity_;
    Evaluation solventMobility_;
    Evaluation solventInvFormationVolumeFactor_;

    Scalar solventRefDensity_;
};

template <class TypeTag>
class BlackOilSolventIntensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;


public:
    void solventPreSatFuncUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                                  unsigned scvIdx OPM_UNUSED,
                                  unsigned timeIdx OPM_UNUSED)
    { }

    void solventPostSatFuncUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                                   unsigned scvIdx OPM_UNUSED,
                                   unsigned timeIdx OPM_UNUSED)
    { }

    void solventPvtUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                           unsigned scvIdx OPM_UNUSED,
                           unsigned timeIdx OPM_UNUSED)
    { }

    const Evaluation& solventSaturation() const
    { OPM_THROW(std::runtime_error, "solventSaturation() called but solvents are disabled"); }

    const Evaluation& solventDensity() const
    { OPM_THROW(std::runtime_error, "solventDensity() called but solvents are disabled"); }

    const Evaluation& solventViscosity() const
    { OPM_THROW(std::runtime_error, "solventViscosity() called but solvents are disabled"); }

    const Evaluation& solventMobility() const
    { OPM_THROW(std::runtime_error, "solventMobility() called but solvents are disabled"); }

    const Evaluation& solventInverseFormationVolumeFactor() const
     { OPM_THROW(std::runtime_error, "solventInverseFormationVolumeFactor() called but solvents are disabled"); }

    const Scalar& solventRefDensity() const
     { OPM_THROW(std::runtime_error, "solventRefDensity() called but solvents are disabled"); }
};

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilSolventExtensiveQuantities
 *
 * \brief Provides the solvent specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableSolventV = GET_PROP_VALUE(TypeTag, EnableSolvent)>
class BlackOilSolventExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

public:
    /*!
     * \brief Method which calculates the volume flux of the polymer "phase" using the
     *        pressure potential gradient of the gas phase and the intrinsic permeability
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateVolumeFluxPerm(const ElementContext& elemCtx,
                              unsigned scvfIdx,
                              unsigned timeIdx)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        Ewoms::PressureCallback<TypeTag> pressureCallback(elemCtx);

        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const auto& faceNormal = scvf.normal();

        unsigned i = scvf.interiorIndex();
        unsigned j = scvf.exteriorIndex();

        // calculate the "raw" pressure gradient
        DimEvalVector solventPGrad;
        pressureCallback.setPhaseIndex(gasPhaseIdx);
        gradCalc.calculateGradient(solventPGrad,
                                   elemCtx,
                                   scvfIdx,
                                   pressureCallback);
        Opm::Valgrind::CheckDefined(solventPGrad);

        // correct the pressure gradients by the gravitational acceleration
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            const auto& gIn = elemCtx.problem().gravity(elemCtx, i, timeIdx);
            const auto& gEx = elemCtx.problem().gravity(elemCtx, j, timeIdx);

            const auto& intQuantsIn = elemCtx.intensiveQuantities(i, timeIdx);
            const auto& intQuantsEx = elemCtx.intensiveQuantities(j, timeIdx);

            const auto& posIn = elemCtx.pos(i, timeIdx);
            const auto& posEx = elemCtx.pos(j, timeIdx);
            const auto& posFace = scvf.integrationPos();

            // the distance between the centers of the control volumes
            DimVector distVecIn(posIn);
            DimVector distVecEx(posEx);
            DimVector distVecTotal(posEx);

            distVecIn -= posFace;
            distVecEx -= posFace;
            distVecTotal -= posIn;
            Scalar absDistTotalSquared = distVecTotal.two_norm2();

            // calculate the hydrostatic pressure at the integration point of the face
            auto rhoIn = intQuantsIn.solventDensity();
            auto pStatIn = - rhoIn*(gIn*distVecIn);

            // the quantities on the exterior side of the face do not influence the
            // result for the TPFA scheme, so they can be treated as scalar values.
            Scalar rhoEx = Toolbox::value(intQuantsEx.solventDensity());
            Scalar pStatEx = - rhoEx*(gEx*distVecEx);

            // compute the hydrostatic gradient between the two control volumes (this
            // gradient exhibitis the same direction as the vector between the two
            // control volume centers and the length (pStaticExterior -
            // pStaticInterior)/distanceInteriorToExterior
            DimEvalVector f(distVecTotal);
            f *= (pStatEx - pStatIn)/absDistTotalSquared;

            // calculate the final potential gradient
            for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                solventPGrad[dimIdx] += f[dimIdx];

                if (!Opm::isfinite(solventPGrad[dimIdx])) {
                    OPM_THROW(Opm::NumericalProblem,
                              "Non-finite potential gradient for solvent 'phase'");
                }
            }
        }

        // determine the upstream and downstream DOFs
        Evaluation solventPGradNormal = 0.0;
        for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
            solventPGradNormal += solventPGrad[dimIdx]*faceNormal[dimIdx];

        if (solventPGradNormal > 0) {
            solventUpstreamDofIdx_ = j;
            solventDownstreamDofIdx_ = i;
        }
        else {
            solventUpstreamDofIdx_ = i;
            solventDownstreamDofIdx_ = j;
        }

        const auto& up = elemCtx.intensiveQuantities(solventUpstreamDofIdx_, timeIdx);

        // this is also slightly hacky because it assumes that the derivative of the
        // flux between two DOFs only depends on the primary variables in the
        // upstream direction. For non-TPFA flux approximation schemes, this is not
        // true...
        if (solventUpstreamDofIdx_ == i)
            solventVolumeFlux_ = solventPGradNormal*up.solventMobility();
        else
            solventVolumeFlux_ = solventPGradNormal*Opm::scalarValue(up.solventMobility());
    }

    /*!
     * \brief Method which calculates the volume flux of the polymer "phase" using the
     *        gas pressure potential difference between cells and transmissibilities
     */
    template <class Dummy = bool> // we need to make this method a template to avoid
                                  // compiler errors if it is not instantiated!
    void updateVolumeFluxTrans(const ElementContext& elemCtx,
                               unsigned scvfIdx,
                               unsigned timeIdx)
    {
        const ExtensiveQuantities& extQuants = asImp_();

        unsigned interiorDofIdx = extQuants.interiorIndex();
        unsigned exteriorDofIdx = extQuants.exteriorIndex();
        assert(interiorDofIdx != exteriorDofIdx);

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

        unsigned I = elemCtx.globalSpaceIndex(interiorDofIdx, timeIdx);
        unsigned J = elemCtx.globalSpaceIndex(exteriorDofIdx, timeIdx);

        Scalar thpres = elemCtx.problem().thresholdPressure(I, J);
        Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        Scalar zIn = elemCtx.problem().dofCenterDepth(elemCtx, interiorDofIdx, timeIdx);
        Scalar zEx = elemCtx.problem().dofCenterDepth(elemCtx, exteriorDofIdx, timeIdx);
        Scalar distZ = zIn - zEx;

        const Evaluation& rhoIn = intQuantsIn.solventDensity();
        Scalar rhoEx = Toolbox::value(intQuantsEx.solventDensity());
        const Evaluation& rhoAvg = rhoIn*0.5 + rhoEx*0.5;

        const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(gasPhaseIdx);
        Evaluation pressureExterior = Toolbox::value(intQuantsEx.fluidState().pressure(gasPhaseIdx));
        pressureExterior += distZ*g*rhoAvg;

        Evaluation pressureDiffSolvent = pressureExterior - pressureInterior;
        if (std::abs(Opm::scalarValue(pressureDiffSolvent)) > thpres) {
            if (pressureDiffSolvent < 0.0)
                pressureDiffSolvent += thpres;
            else
                pressureDiffSolvent -= thpres;
        }
        else
            pressureDiffSolvent = 0.0;

        if (pressureDiffSolvent > 0.0) {
            solventUpstreamDofIdx_ = exteriorDofIdx;
            solventDownstreamDofIdx_ = interiorDofIdx;
        }
        else if (pressureDiffSolvent < 0.0) {
            solventUpstreamDofIdx_ = interiorDofIdx;
            solventDownstreamDofIdx_ = exteriorDofIdx;
        }
        else {
            // pressure potential gradient is zero; force consistent upstream and
            // downstream indices over the intersection regardless of the side which it
            // is looked at.
            solventUpstreamDofIdx_ = std::min(interiorDofIdx, exteriorDofIdx);
            solventDownstreamDofIdx_ = std::max(interiorDofIdx, exteriorDofIdx);
            solventVolumeFlux_ = 0.0;
            return;
        }

        Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
        const IntensiveQuantities& up = elemCtx.intensiveQuantities(solventUpstreamDofIdx_, timeIdx);
        if (solventUpstreamDofIdx_ == interiorDofIdx)
            solventVolumeFlux_ =
                up.solventMobility()
                *(-trans/faceArea)
                *pressureDiffSolvent;
        else
            solventVolumeFlux_ =
                Opm::scalarValue(up.solventMobility())
                *(-trans/faceArea)
                *pressureDiffSolvent;
    }

    unsigned solventUpstreamIndex() const
    { return solventUpstreamDofIdx_; }

    unsigned solventDownstreamIndex() const
    { return solventDownstreamDofIdx_; }

    const Evaluation& solventVolumeFlux() const
    { return solventVolumeFlux_; }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation solventVolumeFlux_;
    unsigned solventUpstreamDofIdx_;
    unsigned solventDownstreamDofIdx_;
};

template <class TypeTag>
class BlackOilSolventExtensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

public:
    void updateVolumeFluxPerm(const ElementContext& elemCtx OPM_UNUSED,
                              unsigned scvfIdx OPM_UNUSED,
                              unsigned timeIdx OPM_UNUSED)
    { }

    void updateVolumeFluxTrans(const ElementContext& elemCtx OPM_UNUSED,
                              unsigned scvfIdx OPM_UNUSED,
                              unsigned timeIdx OPM_UNUSED)
    { }

    unsigned solventUpstreamIndex() const
    { OPM_THROW(std::runtime_error, "solventUpstreamIndex() called but solvents are disabled"); }

    unsigned solventDownstreamIndex() const
    { OPM_THROW(std::runtime_error, "solventDownstreamIndex() called but solvents are disabled"); }

    const Evaluation& solventVolumeFlux() const
    { OPM_THROW(std::runtime_error, "solventVolumeFlux() called but solvents are disabled"); }
};

} // namespace Ewoms

#endif
