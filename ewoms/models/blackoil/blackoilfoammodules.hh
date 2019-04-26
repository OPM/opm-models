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
 * \brief Contains the classes required to extend the black-oil model to include the effects of foam.
 */
#ifndef EWOMS_BLACK_OIL_FOAM_MODULE_HH
#define EWOMS_BLACK_OIL_FOAM_MODULE_HH

#include "blackoilproperties.hh"
//#include <ewoms/io/vtkblackoilfoammodule.hh>
#include <ewoms/models/common/quantitycallbacks.hh>

//#include <opm/material/common/Tabulated1DFunction.hpp>
//#include <opm/material/common/IntervalTabulated2DFunction.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>
#include <math.h>

namespace Ewoms {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model to include the effects of foam.
 */
template <class TypeTag, bool enableFoamV = GET_PROP_VALUE(TypeTag, EnableFoam)>
class BlackOilFoamModule
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

    //typedef typename Opm::Tabulated1DFunction<Scalar> TabulatedFunction;
    //typedef typename Opm::IntervalTabulated2DFunction<Scalar> TabulatedTwoDFunction;

    static constexpr unsigned foamConcentrationIdx = Indices::foamConcentrationIdx;
    static constexpr unsigned contiFoamEqIdx = Indices::contiFoamEqIdx;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;


    static constexpr unsigned enableFoam = enableFoamV;

    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    //enum AdsorptionBehaviour { Desorption = 1, NoDesorption = 2 };

    // a struct containing constants to calculate change to relative permeability,
    // based on model (1-9) in Table 1 of
    // Kun  Ma,  Guangwei  Ren,  Khalid  Mateen,  Danielle  Morel,  and  PhilippeCordelier.
    // Modeling techniques for foam flow in porous media.SPE Journal,20(03):453–470, jun 2015.
    // The constants are provided by the keyword ...
    struct FoamCoefficients {
        Scalar fm_mob = 2;

        Scalar fm_surf = 1;
        Scalar ep_surf = 1;

        Scalar fm_oil = 1;
        Scalar fl_oil = 0;
        Scalar ep_oil = 1;

        Scalar fm_dry = 1;
        Scalar ep_dry = 0;

        Scalar fm_cap = 1;
        Scalar ep_cap = 1;
    };

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the polymer module
     */
    static void initFromDeck(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        // some sanity checks: if foam is enabled, the FOAM keyword must be
        // present, if foam is disabled the keyword must not be present.
        if (enableFoam && !deck.hasKeyword("FOAM")) {
            throw std::runtime_error("Non-trivial foam treatment requested at compile time, but "
                                     "the deck does not contain the FOAM keyword");
        }
        else if (!enableFoam && deck.hasKeyword("FOAM")) {
            throw std::runtime_error("Foam treatment disabled at compile time, but the deck "
                                     "contains the FOAM keyword");
        }

        if (!deck.hasKeyword("FOAM"))
            return; // foam treatment is supposed to be disabled

        //unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        // TODO: do not hard code
        setNumSatRegions(1);

    }
#endif

    // a struct containing the constants to calculate foam viscosity
    // based on Mark-Houwink equation and Huggins equation, the constants are provided
    // by the keyword PLYVMH
    //struct PlyvmhCoefficients {
    //    Scalar k_mh;
    //    Scalar a_mh;
    //    Scalar gamma;
    //    Scalar kappa;
    //};

    /*!
     * \brief Specify the number of satuation regions.
     *
     * This must be called before setting the PLYROCK and PLYADS of any region.
     */
    static void setNumSatRegions(unsigned numRegions)
    {
        foamCoefficients_.resize(numRegions);
    //    plyrockDeadPoreVolume_.resize(numRegions);
    //    plyrockResidualResistanceFactor_.resize(numRegions);
    //    plyrockRockDensityFactor_.resize(numRegions);
    //    plyrockAdsorbtionIndex_.resize(numRegions);
    //    plyrockMaxAdsorbtion_.resize(numRegions);
    //    plyadsAdsorbedFoam_.resize(numRegions);
    }

    /*!
     * \brief Specify the foam properties a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    //static void setPlyrock(unsigned satRegionIdx,
    //                       const Scalar& plyrockDeadPoreVolume,
    //                       const Scalar& plyrockResidualResistanceFactor,
    //                       const Scalar& plyrockRockDensityFactor,
    //                       const Scalar& plyrockAdsorbtionIndex,
    //                       const Scalar& plyrockMaxAdsorbtion)
    //{
    //    plyrockDeadPoreVolume_[satRegionIdx] = plyrockDeadPoreVolume;
    //    plyrockResidualResistanceFactor_[satRegionIdx] = plyrockResidualResistanceFactor;
    //    plyrockRockDensityFactor_[satRegionIdx] = plyrockRockDensityFactor;
    //    plyrockAdsorbtionIndex_[satRegionIdx] = plyrockAdsorbtionIndex;
    //    plyrockMaxAdsorbtion_[satRegionIdx] = plyrockMaxAdsorbtion;
    //}

    /*!
     * \brief Specify the number of pvt regions.
     *
     * This must be called before setting the PLYVISC of any region.
     */
    //static void setNumPvtRegions(unsigned numRegions)
    //{
    //    plyviscViscosityMultiplierTable_.resize(numRegions);
    //}

    /*!
     * \brief Specify the foam viscosity a single region.
     *
     * The index of specified here must be in range [0, numSatRegions)
     */
    //static void setPlyvisc(unsigned satRegionIdx,
    //                       const TabulatedFunction& plyviscViscosityMultiplierTable)
    //{
    //    plyviscViscosityMultiplierTable_[satRegionIdx] = plyviscViscosityMultiplierTable;
    //}

    /*!
     * \brief Specify the number of mix regions.
     *
     * This must be called before setting the PLYMAC and PLMIXPAR of any region.
     */
    //static void setNumMixRegions(unsigned numRegions)
    //{
    //    plymaxMaxConcentration_.resize(numRegions);
    //    plymixparToddLongstaff_.resize(numRegions);

    //}

    /*!
     * \brief Specify the maximum foam concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    //static void setPlymax(unsigned mixRegionIdx,
    //                      const Scalar& plymaxMaxConcentration)
    //{
    //    plymaxMaxConcentration_[mixRegionIdx] = plymaxMaxConcentration;
    //}

    /*!
     * \brief Specify the maximum foam concentration a single region.
     *
     * The index of specified here must be in range [0, numMixRegionIdx)
     */
    //static void setPlmixpar(unsigned mixRegionIdx,
    //                        const Scalar& plymixparToddLongstaff)
    //{
    //    plymixparToddLongstaff_[mixRegionIdx] = plymixparToddLongstaff;
    //}

    /*!
     * \brief Register all run-time parameters for the black-oil foam module.
     */
    static void registerParameters()
    {
        if (!enableFoam)
            // foam has been disabled at compile time
            return;

        //Ewoms::VtkBlackOilFoamModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all foam specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableFoam)
            // foam have been disabled at compile time
            return;

        //model.addOutputModule(new Ewoms::VtkBlackOilFoamModule<TypeTag>(simulator));
    }

    //static bool primaryVarApplies(unsigned pvIdx)
    //{
    //    // foam has been disabled at compile time
    //    return enableFoam;
    //}

    //static std::string primaryVarName(unsigned pvIdx)
    //{
    //    assert(primaryVarApplies(pvIdx));

    //    if (pvIdx == polymerConcentrationIdx) {
    //        return "polymer_waterconcentration";
    //    }
    //    else {
    //        return "polymer_molecularweight";
    //    }
    //}

    //static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    //{
    //    assert(primaryVarApplies(pvIdx));

    //    // TODO: it may be beneficial to chose this differently.
    //    return static_cast<Scalar>(1.0);
    //}

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableFoam)
            return false;

        return eqIdx == contiFoamEqIdx;

    }

    static std::string eqName(unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        return "conti^foam";
    }

    //static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    //{
    //    assert(eqApplies(eqIdx));

    //    // TODO: it may be beneficial to chose this differently.
    //    return static_cast<Scalar>(1.0);
    //}

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableFoam)
            return;

        const auto& fs = intQuants.fluidState();

        LhsEval surfaceVolumeFreeGas =
            Toolbox::template decay<LhsEval>(fs.saturation(gasPhaseIdx))
            * Toolbox::template decay<LhsEval>(fs.invB(gasPhaseIdx))
            * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // Avoid singular matrix if no gas is present.
        surfaceVolumeFreeGas = Opm::max(surfaceVolumeFreeGas, 1e-10);

        // Foam/surfactant in gas phase.
        const LhsEval gasFoam = surfaceVolumeFreeGas
            * Toolbox::template decay<LhsEval>(br) // TODO
            * Toolbox::template decay<LhsEval>(intQuants.foamConcentration());

        // Adsorbed foam/surfactant.
        const LhsEval adsorbedFoam =
            Toolbox::template decay<LhsEval>(1.0 - intQuants.porosity())
            * Toolbox::template decay<LhsEval>(intQuants.foamRockDensity()) // TODO
            * Toolbox::template decay<LhsEval>(intQuants.foamAdsorption())  // TODO
            / Toolbox::template decay<LhsEval>(intQuants.porosity());

        LhsEval accumulationFoam = gasFoam + adsorbedFoam;
        storage[contiFoamEqIdx] += accumulationFoam;
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)

    {
        if (!enableFoam)
            return;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const unsigned upIdx = extQuants.upstreamIndex(FluidSystem::waterPhaseIdx);
        const unsigned inIdx = extQuants.interiorIndex();
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const unsigned contiWaterEqIdx = Indices::conti0EqIdx + Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);


        if (upIdx == inIdx) {
            flux[contiFoamEqIdx] = 0.;
                    //extQuants.volumeFlux(waterPhaseIdx)
                    //*up.fluidState().invB(waterPhaseIdx)
                    //*up.polymerViscosityCorrection()
                    ///extQuants.polymerShearFactor()
                    //*up.polymerConcentration();

            // modify water
            flux[contiWaterEqIdx] /= 1.;
                    //extQuants.waterShearFactor();
        }
        else {
            flux[contiFoamEqIdx] = 0.;
                    //extQuants.volumeFlux(waterPhaseIdx)
                    //*Opm::decay<Scalar>(up.fluidState().invB(waterPhaseIdx))
                    //*Opm::decay<Scalar>(up.polymerViscosityCorrection())
                    ///Opm::decay<Scalar>(extQuants.polymerShearFactor())
                    //*Opm::decay<Scalar>(up.polymerConcentration());

            // modify water
            flux[contiWaterEqIdx] /= 1.;
                    //Opm::decay<Scalar>(extQuants.waterShearFactor());
        }

    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider the change of foam primary variables for convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableFoam)
            return;

        //unsigned dofIdx = model.dofMapper().index(dof);
        //const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        //outstream << priVars[polymerConcentrationIdx];
        //outstream << priVars[polymerMoleWeightIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableFoam)
            return;

        //unsigned dofIdx = model.dofMapper().index(dof);
        //PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        //PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        //instream >> priVars0[polymerConcentrationIdx];
        //instream >> priVars0[polymerMoleWeightIdx];

        //// set the primary variables for the beginning of the current time step.
        //priVars1[polymerConcentrationIdx] = priVars0[polymerConcentrationIdx];
        //priVars1[polymerMoleWeightIdx] = priVars0[polymerMoleWeightIdx];
    }

    //static const Scalar plyrockDeadPoreVolume(const ElementContext& elemCtx,
    //                                          unsigned scvIdx,
    //                                          unsigned timeIdx)
    //{
    //    unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyrockDeadPoreVolume_[satnumRegionIdx];
    //}

    //static const Scalar plyrockResidualResistanceFactor(const ElementContext& elemCtx,
    //                                                    unsigned scvIdx,
    //                                                    unsigned timeIdx)
    //{
    //    unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyrockResidualResistanceFactor_[satnumRegionIdx];
    //}

    //static const Scalar plyrockRockDensityFactor(const ElementContext& elemCtx,
    //                                             unsigned scvIdx,
    //                                             unsigned timeIdx)
    //{
    //    unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyrockRockDensityFactor_[satnumRegionIdx];
    //}

    //static const Scalar plyrockAdsorbtionIndex(const ElementContext& elemCtx,
    //                                           unsigned scvIdx,
    //                                           unsigned timeIdx)
    //{
    //    unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyrockAdsorbtionIndex_[satnumRegionIdx];
    //}

    //static const Scalar plyrockMaxAdsorbtion(const ElementContext& elemCtx,
    //                                         unsigned scvIdx,
    //                                         unsigned timeIdx)
    //{
    //    unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyrockMaxAdsorbtion_[satnumRegionIdx];
    //}

    //static const TabulatedFunction& plyadsAdsorbedPolymer(const ElementContext& elemCtx,
    //                                                      unsigned scvIdx,
    //                                                      unsigned timeIdx)
    //{
    //    unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyadsAdsorbedPolymer_[satnumRegionIdx];
    //}

    //static const TabulatedFunction& plyviscViscosityMultiplierTable(const ElementContext& elemCtx,
    //                                                                unsigned scvIdx,
    //                                                                unsigned timeIdx)
    //{
    //    unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    //}

    //static const TabulatedFunction& plyviscViscosityMultiplierTable(unsigned pvtnumRegionIdx)
    //{
    //    return plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    //}

    //static const Scalar plymaxMaxConcentration(const ElementContext& elemCtx,
    //                                           unsigned scvIdx,
    //                                           unsigned timeIdx)
    //{
    //    unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plymaxMaxConcentration_[polymerMixRegionIdx];
    //}

    //static const Scalar plymixparToddLongstaff(const ElementContext& elemCtx,
    //                                           unsigned scvIdx,
    //                                           unsigned timeIdx)
    //{
    //    unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plymixparToddLongstaff_[polymerMixRegionIdx];
    //}

    static const FoamCoefficients& foamCoefficients(const ElementContext& elemCtx,
                                                        const unsigned scvIdx,
                                                        const unsigned timeIdx)
    {
        //const unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
        // TODO: Don't hard code!!!
        return foamCoefficients_[0];
    }

    //static const PlyvmhCoefficients& plyvmhCoefficients(const ElementContext& elemCtx,
    //                                                    const unsigned scvIdx,
    //                                                    const unsigned timeIdx)
    //{
    //    const unsigned polymerMixRegionIdx = elemCtx.problem().plmixnumRegionIndex(elemCtx, scvIdx, timeIdx);
    //    return plyvmhCoefficients_[polymerMixRegionIdx];
    //}

    ///*!
    // * \brief Computes the shear factor
    // *
    // * Input is polymer concentration and either the water velocity or the shrate if hasShrate_ is true.
    // * The pvtnumRegionIdx is needed to make sure the right table is used.
    // */
    //template <class Evaluation>
    //static Evaluation computeShearFactor(const Evaluation& polymerConcentration,
    //                                     unsigned pvtnumRegionIdx,
    //                                     const Evaluation& v0)
    //{
    //    typedef Opm::MathToolbox<Evaluation> ToolboxLocal;

    //    const auto& viscosityMultiplierTable = plyviscViscosityMultiplierTable_[pvtnumRegionIdx];
    //    Scalar viscosityMultiplier = viscosityMultiplierTable.eval(Opm::scalarValue(polymerConcentration), /*extrapolate=*/true);

    //    const Scalar eps = 1e-14;
    //    // return 1.0 if the polymer has no effect on the water.
    //    if (std::abs((viscosityMultiplier - 1.0)) < eps){
    //        return ToolboxLocal::createBlank(v0) + 1.;
    //    }

    //    const std::vector<Scalar>& shearEffectRefLogVelocity = plyshlogShearEffectRefLogVelocity_[pvtnumRegionIdx];
    //    auto v0AbsLog = Opm::log(Opm::abs(v0));
    //    // return 1.0 if the velocity /sharte is smaller than the first velocity entry.
    //    if (v0AbsLog < shearEffectRefLogVelocity[0])
    //        return ToolboxLocal::createBlank(v0) + 1.0;

    //    // compute shear factor from input
    //    // Z = (1 + (P - 1) * M(v)) / P
    //    // where M(v) is computed from user input
    //    // and P = viscosityMultiplier
    //    const std::vector<Scalar>& shearEffectRefMultiplier = plyshlogShearEffectRefMultiplier_[pvtnumRegionIdx];
    //    size_t numTableEntries = shearEffectRefLogVelocity.size();
    //    assert(shearEffectRefMultiplier.size() == numTableEntries);

    //    std::vector<Scalar> shearEffectMultiplier(numTableEntries, 1.0);
    //    for (size_t i = 0; i < numTableEntries; ++i) {
    //        shearEffectMultiplier[i] = (1.0 + (viscosityMultiplier - 1.0)*shearEffectRefMultiplier[i]) / viscosityMultiplier;
    //        shearEffectMultiplier[i] = Opm::log(shearEffectMultiplier[i]);
    //    }
    //    // store the logarithmic velocity and logarithmic multipliers in a table for easy look up and
    //    // linear interpolation in the logarithmic space.
    //    TabulatedFunction logShearEffectMultiplier = TabulatedFunction(numTableEntries, shearEffectRefLogVelocity, shearEffectMultiplier, /*bool sortInputs =*/ false);

    //    // Find sheared velocity (v) that satisfies
    //    // F = log(v) + log (Z) - log(v0) = 0;

    //    // Set up the function
    //    // u = log(v)
    //    auto F = [&logShearEffectMultiplier, &v0AbsLog](const Evaluation& u) {
    //        return u + logShearEffectMultiplier.eval(u, true) - v0AbsLog;
    //    };
    //    // and its derivative
    //    auto dF = [&logShearEffectMultiplier](const Evaluation& u) {
    //        return 1 + logShearEffectMultiplier.evalDerivative(u, true);
    //    };

    //    // Solve F = 0 using Newton
    //    // Use log(v0) as initial value for u
    //    auto u = v0AbsLog;
    //    bool converged = false;
    //    for (int i = 0; i < 20; ++i) {
    //        auto f = F(u);
    //        auto df = dF(u);
    //        u -= f/df;
    //        if (std::abs(Opm::scalarValue(f)) < 1e-12) {
    //            converged = true;
    //            break;
    //        }
    //    }
    //    if (!converged) {
    //        throw std::runtime_error("Not able to compute shear velocity. \n");
    //    }

    //    // return the shear factor
    //    return Opm::exp(logShearEffectMultiplier.eval(u, /*extrapolate=*/true));

    //}

    //const Scalar molarMass() const
    //{
    //    return 0.25; // kg/mol
    //}




private:
    //static std::vector<Scalar> plyrockDeadPoreVolume_;
    //static std::vector<Scalar> plyrockResidualResistanceFactor_;
    //static std::vector<Scalar> plyrockRockDensityFactor_;
    //static std::vector<Scalar> plyrockAdsorbtionIndex_;
    //static std::vector<Scalar> plyrockMaxAdsorbtion_;
    //static std::vector<TabulatedFunction> plyadsAdsorbedPolymer_;
    //static std::vector<TabulatedFunction> plyviscViscosityMultiplierTable_;
    //static std::vector<Scalar> plymaxMaxConcentration_;
    //static std::vector<Scalar> plymixparToddLongstaff_;
    //static std::vector<std::vector<Scalar>> plyshlogShearEffectRefMultiplier_;
    //static std::vector<std::vector<Scalar>> plyshlogShearEffectRefLogVelocity_;
    //static std::vector<Scalar> shrate_;
    //static bool hasShrate_;
    //static bool hasPlyshlog_;

    static std::vector<FoamCoefficients> foamCoefficients_;
    //static std::vector<PlyvmhCoefficients> plyvmhCoefficients_;
    //static std::map<int, TabulatedTwoDFunction> plymwinjTables_;
    //static std::map<int, TabulatedTwoDFunction> skprwatTables_;

    //static std::map<int, SkprpolyTable> skprpolyTables_;
};



//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockDeadPoreVolume_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockResidualResistanceFactor_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockRockDensityFactor_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockAdsorbtionIndex_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyrockMaxAdsorbtion_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedFunction>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyadsAdsorbedPolymer_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedFunction>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyviscViscosityMultiplierTable_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plymaxMaxConcentration_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plymixparToddLongstaff_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyshlogShearEffectRefMultiplier_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyshlogShearEffectRefLogVelocity_;
//
//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::Scalar>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::shrate_;
//
//template <class TypeTag, bool enablePolymerV>
//bool
//BlackOilPolymerModule<TypeTag, enablePolymerV>::hasShrate_;
//
//template <class TypeTag, bool enablePolymerV>
//bool
//BlackOilPolymerModule<TypeTag, enablePolymerV>::hasPlyshlog_;
//

template <class TypeTag, bool enableFoam>
std::vector<typename BlackOilFoamModule<TypeTag, enableFoam>::FoamCoefficients>
BlackOilFoamModule<TypeTag, enableFoam>::foamCoefficients_;

//template <class TypeTag, bool enablePolymerV>
//std::vector<typename BlackOilPolymerModule<TypeTag, enablePolymerV>::PlyvmhCoefficients>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plyvmhCoefficients_;
//
//template <class TypeTag, bool enablePolymerV>
//std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedTwoDFunction>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::plymwinjTables_;
//
//template <class TypeTag, bool enablePolymerV>
//std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::TabulatedTwoDFunction>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::skprwatTables_;
//
//template <class TypeTag, bool enablePolymerV>
//std::map<int, typename BlackOilPolymerModule<TypeTag, enablePolymerV>::SkprpolyTable>
//BlackOilPolymerModule<TypeTag, enablePolymerV>::skprpolyTables_;

/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilPolymerIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        polymers extension of the black-oil model.
 */
template <class TypeTag, bool enableFoam = GET_PROP_VALUE(TypeTag, EnableFoam)>
class BlackOilFoamIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilFoamModule<TypeTag> FoamModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    static constexpr int foamConcentrationIdx = Indices::foamConcentrationIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;


public:

    /*!
     * \brief Update the intensive properties needed to handle polymers from the
     *        primary variables
     *
     */
    void foamPropertiesUpdate_(const ElementContext& elemCtx,
                                  unsigned dofIdx,
                                  unsigned timeIdx)
    {
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        foamConcentration_ = priVars.makeEvaluation(foamConcentrationIdx, timeIdx);
        //const Scalar cmax = PolymerModule::plymaxMaxConcentration(elemCtx, dofIdx, timeIdx);

        //// permeability reduction due to polymer
        //const Scalar& maxAdsorbtion = PolymerModule::plyrockMaxAdsorbtion(elemCtx, dofIdx, timeIdx);
        //const auto& plyadsAdsorbedPolymer = PolymerModule::plyadsAdsorbedPolymer(elemCtx, dofIdx, timeIdx);
        //polymerAdsorption_ = plyadsAdsorbedPolymer.eval(polymerConcentration_, /*extrapolate=*/true);
        //if (PolymerModule::plyrockAdsorbtionIndex(elemCtx, dofIdx, timeIdx) == PolymerModule::NoDesorption) {
        //    const Scalar& maxPolymerAdsorption = elemCtx.problem().maxPolymerAdsorption(elemCtx, dofIdx, timeIdx);
        //    polymerAdsorption_ = std::max(Evaluation(maxPolymerAdsorption) , polymerAdsorption_);
        //}

        //// compute resitanceFactor
        //const Scalar& residualResistanceFactor = PolymerModule::plyrockResidualResistanceFactor(elemCtx, dofIdx, timeIdx);
        //const Evaluation resistanceFactor = 1.0 + (residualResistanceFactor - 1.0) * polymerAdsorption_ / maxAdsorbtion;

        //// compute effective viscosities
        //if (!enablePolymerMolarWeight) {
        //    const auto& fs = asImp_().fluidState_;
        //    const Evaluation& muWater = fs.viscosity(waterPhaseIdx);
        //    const auto& viscosityMultiplier = PolymerModule::plyviscViscosityMultiplierTable(elemCtx, dofIdx, timeIdx);
        //    const Evaluation viscosityMixture = viscosityMultiplier.eval(polymerConcentration_, /*extrapolate=*/true) * muWater;

        //    // Do the Todd-Longstaff mixing
        //    const Scalar plymixparToddLongstaff = PolymerModule::plymixparToddLongstaff(elemCtx, dofIdx, timeIdx);
        //    const Evaluation viscosityPolymer = viscosityMultiplier.eval(cmax, /*extrapolate=*/true) * muWater;
        //    const Evaluation viscosityPolymerEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(viscosityPolymer, 1.0 - plymixparToddLongstaff);
        //    const Evaluation viscosityWaterEffective = pow(viscosityMixture, plymixparToddLongstaff) * pow(muWater, 1.0 - plymixparToddLongstaff);

        //    const Evaluation cbar = polymerConcentration_ / cmax;
        //    // waterViscosity / effectiveWaterViscosity
        //    waterViscosityCorrection_ = muWater * ((1.0 - cbar) / viscosityWaterEffective + cbar / viscosityPolymerEffective);
        //    // effectiveWaterViscosity / effectivePolymerViscosity
        //    polymerViscosityCorrection_ =  (muWater / waterViscosityCorrection_) / viscosityPolymerEffective;
        //}
        //else { // based on PLYVMH
        //    const auto& plyvmhCoefficients = PolymerModule::plyvmhCoefficients(elemCtx, dofIdx, timeIdx);
        //    const Scalar k_mh = plyvmhCoefficients.k_mh;
        //    const Scalar a_mh = plyvmhCoefficients.a_mh;
        //    const Scalar gamma = plyvmhCoefficients.gamma;
        //    const Scalar kappa = plyvmhCoefficients.kappa;

        //    // viscosity model based on Mark-Houwink equation and Huggins equation
        //    // 1000 is a emperical constant, most likely related to unit conversion
        //    const Evaluation intrinsicViscosity = k_mh * pow(polymerMoleWeight_ * 1000., a_mh);
        //    const Evaluation x = polymerConcentration_ * intrinsicViscosity;
        //    waterViscosityCorrection_ = 1.0 / (1.0 + gamma * (x + kappa * x * x));
        //    polymerViscosityCorrection_ = 1.0;
        //}

        const auto& foamCoefficients = FoamModule::foamCoefficients(elemCtx, dofIdx, timeIdx);

        const Scalar fm_mob = foamCoefficients.fm_mob;

        const Scalar fm_surf = foamCoefficients.fm_surf;
        const Scalar ep_surf = foamCoefficients.ep_surf;

        const Scalar fm_oil = foamCoefficients.fm_oil;
        const Scalar fl_oil = foamCoefficients.fl_oil;
        const Scalar ep_oil = foamCoefficients.ep_oil;

        const Scalar fm_dry = foamCoefficients.fm_dry;
        const Scalar ep_dry = foamCoefficients.ep_dry;

        const Scalar fm_cap = foamCoefficients.fm_cap;
        const Scalar ep_cap = foamCoefficients.ep_cap;

        const auto& fs = asImp_().fluidState_;
        const Evaluation C_surf = foamConcentration_;
        const Evaluation Ca = 1e10; // TODO: replace with proper capillary number.
        const Evaluation S_o = fs.saturation(oilPhaseIdx);
        const Evaluation S_w = fs.saturation(waterPhaseIdx);

        Evaluation F1 = pow(C_surf/fm_surf, ep_surf);
        Evaluation F2 = pow((fm_oil-S_o)/(fm_oil-fl_oil), ep_oil);
        Evaluation F3 = pow(fm_cap/Ca, ep_cap);
        Evaluation F7 = 0.5 + atan(ep_dry*(S_w-fm_dry))/M_PI;

        Evaluation mobilityReductionFactor = 1./(1. + fm_mob*F1*F2*F3*F7);

        // adjust gas mobility
        asImp_().mobility_[gasPhaseIdx] *= mobilityReductionFactor;

        //// update rock properties
        //polymerDeadPoreVolume_ = PolymerModule::plyrockDeadPoreVolume(elemCtx, dofIdx, timeIdx);
        //polymerRockDensity_ = PolymerModule::plyrockRockDensityFactor(elemCtx, dofIdx, timeIdx);
    }

    //const Evaluation& polymerConcentration() const
    //{ return polymerConcentration_; }
    const Evaluation& foamConcentration() const
    { return foamConcentration_; }

    //const Evaluation& polymerMoleWeight() const
    //{
    //    if (!enablePolymerMolarWeight)
    //        throw std::logic_error("polymerMoleWeight() is called but polymer milecular weight is disabled");

    //    return polymerMoleWeight_;
    //}

    //const Scalar& polymerDeadPoreVolume() const
    //{ return polymerDeadPoreVolume_; }

    //const Evaluation& polymerAdsorption() const
    //{ return polymerAdsorption_; }

    //const Scalar& polymerRockDensity() const
    //{ return polymerRockDensity_; }

    //// effectiveWaterViscosity / effectivePolymerViscosity
    //const Evaluation& polymerViscosityCorrection() const
    //{ return polymerViscosityCorrection_; }

    //// waterViscosity / effectiveWaterViscosity
    //const Evaluation& waterViscosityCorrection() const
    //{ return waterViscosityCorrection_; }


protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation foamConcentration_;
    //// polymer molecular weight
    //Evaluation polymerMoleWeight_;
    //Scalar polymerDeadPoreVolume_;
    //Scalar polymerRockDensity_;
    //Evaluation polymerAdsorption_;
    //Evaluation polymerViscosityCorrection_;
    //Evaluation waterViscosityCorrection_;


};

template <class TypeTag>
class BlackOilFoamIntensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    void foamPropertiesUpdate_(const ElementContext& elemCtx OPM_UNUSED,
                                  unsigned scvIdx OPM_UNUSED,
                                  unsigned timeIdx OPM_UNUSED)
    { }

    //const Evaluation& polymerMoleWeight() const
    //{ throw std::logic_error("polymerMoleWeight() called but polymer molecular weight is disabled"); }

    //const Evaluation& polymerConcentration() const
    //{ throw std::runtime_error("polymerConcentration() called but polymers are disabled"); }

    const Evaluation& foamConcentration() const
    { throw std::runtime_error("foamConcentration() called but foam is disabled"); }

    //const Evaluation& polymerDeadPoreVolume() const
    //{ throw std::runtime_error("polymerDeadPoreVolume() called but polymers are disabled"); }

    //const Evaluation& polymerAdsorption() const
    //{ throw std::runtime_error("polymerAdsorption() called but polymers are disabled"); }

    //const Evaluation& polymerRockDensity() const
    //{ throw std::runtime_error("polymerRockDensity() called but polymers are disabled"); }

    //const Evaluation& polymerViscosityCorrection() const
    //{ throw std::runtime_error("polymerViscosityCorrection() called but polymers are disabled"); }

    //const Evaluation& waterViscosityCorrection() const
    //{ throw std::runtime_error("waterViscosityCorrection() called but polymers are disabled"); }
};


/*!
 * \ingroup BlackOil
 * \class Ewoms::BlackOilFoamExtensiveQuantities
 *
 * \brief Provides the foam specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableFoam = GET_PROP_VALUE(TypeTag, EnableFoam)>
class BlackOilFoamExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr unsigned waterPhaseIdx =  FluidSystem::waterPhaseIdx;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef BlackOilFoamModule<TypeTag> FoamModule;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

public:
    ///*!
    // * \brief Method which calculates the shear factor based on flow velocity
    // *
    // * This is the variant of the method which assumes that the problem is specified
    // * using permeabilities, i.e., *not* via transmissibilities.
    // */
    //template <class Dummy = bool> // we need to make this method a template to avoid
    //                              // compiler errors if it is not instantiated!
    //void updateShearMultipliersPerm(const ElementContext& elemCtx OPM_UNUSED,
    //                                unsigned scvfIdx OPM_UNUSED,
    //                                unsigned timeIdx OPM_UNUSED)
    //{
    //    throw std::runtime_error("The extension of the blackoil model for foam is not yet "
    //                             "implemented for problems specified using permeabilities.");
    //}

    ///*!
    // * \brief Method which calculates the shear factor based on flow velocity
    // *
    // * This is the variant of the method which assumes that the problem is specified
    // * using transmissibilities, i.e., *not* via permeabilities.
    // */
    //template <class Dummy = bool> // we need to make this method a template to avoid
    //// compiler errors if it is not instantiated!
    //void updateShearMultipliers(const ElementContext& elemCtx,
    //                            unsigned scvfIdx,
    //                            unsigned timeIdx)
    //{

    //    waterShearFactor_ = 1.0;
    //    polymerShearFactor_ = 1.0;

    //    if (!PolymerModule::hasPlyshlog())
    //        return;

    //    const ExtensiveQuantities& extQuants = asImp_();
    //    unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
    //    unsigned interiorDofIdx = extQuants.interiorIndex();
    //    unsigned exteriorDofIdx = extQuants.exteriorIndex();
    //    const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
    //    const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
    //    const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

    //    // compute water velocity from flux
    //    Evaluation poroAvg = intQuantsIn.porosity()*0.5 + intQuantsEx.porosity()*0.5;
    //    unsigned pvtnumRegionIdx = elemCtx.problem().pvtRegionIndex(elemCtx, scvfIdx, timeIdx);
    //    const Evaluation& Sw = up.fluidState().saturation(waterPhaseIdx);
    //    unsigned cellIdx = elemCtx.globalSpaceIndex(scvfIdx, timeIdx);
    //    const auto& materialLawManager = elemCtx.problem().materialLawManager();
    //    const auto& scaledDrainageInfo =
    //            materialLawManager->oilWaterScaledEpsInfoDrainage(cellIdx);
    //    const Scalar& Swcr = scaledDrainageInfo.Swcr;

    //    // guard against zero porosity and no mobile water
    //    Evaluation denom = Opm::max(poroAvg * (Sw - Swcr), 1e-12);
    //    Evaluation waterVolumeVelocity = extQuants.volumeFlux(waterPhaseIdx) / denom;

    //    // if shrate is specified. Compute shrate based on the water velocity
    //    if (PolymerModule::hasShrate()) {
    //        const Evaluation& relWater = up.relativePermeability(waterPhaseIdx);
    //        Scalar trans = elemCtx.problem().transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
    //        if (trans > 0.0) {
    //            Scalar faceArea = elemCtx.stencil(timeIdx).interiorFace(scvfIdx).area();
    //            auto dist = elemCtx.pos(interiorDofIdx, timeIdx) -  elemCtx.pos(exteriorDofIdx, timeIdx);
    //            // compute permeability from transmissibility.
    //            Scalar absPerm = trans / faceArea * dist.two_norm();
    //            waterVolumeVelocity *=
    //                PolymerModule::shrate(pvtnumRegionIdx)*Opm::sqrt(poroAvg*Sw / (relWater*absPerm));
    //            assert(Opm::isfinite(waterVolumeVelocity));
    //        }
    //    }

    //    // compute share factors for water and polymer
    //    waterShearFactor_ =
    //        PolymerModule::computeShearFactor(up.polymerConcentration(),
    //                                          pvtnumRegionIdx,
    //                                          waterVolumeVelocity);
    //    polymerShearFactor_ =
    //        PolymerModule::computeShearFactor(up.polymerConcentration(),
    //                                          pvtnumRegionIdx,
    //                                          waterVolumeVelocity*up.polymerViscosityCorrection());

    //}

    //const Evaluation& polymerShearFactor() const
    //{ return polymerShearFactor_; }

    //const Evaluation& waterShearFactor() const
    //{ return waterShearFactor_; }


private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    //Evaluation polymerShearFactor_;
    //Evaluation waterShearFactor_;

};

template <class TypeTag>
class BlackOilFoamExtensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

public:
    //void updateShearMultipliers(const ElementContext& elemCtx OPM_UNUSED,
    //                            unsigned scvfIdx OPM_UNUSED,
    //                            unsigned timeIdx OPM_UNUSED)
    //{ }

    //void updateShearMultipliersPerm(const ElementContext& elemCtx OPM_UNUSED,
    //                                unsigned scvfIdx OPM_UNUSED,
    //                                unsigned timeIdx OPM_UNUSED)
    //{ }

    //const Evaluation& polymerShearFactor() const
    //{ throw std::runtime_error("polymerShearFactor() called but polymers are disabled"); }

    //const Evaluation& waterShearFactor() const
    //{ throw std::runtime_error("waterShearFactor() called but polymers are disabled"); }
};


} // namespace Ewoms

#endif
