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
 * \copydoc Ewoms::EclSummaryWriter
 */
#ifndef EWOMS_ECL_SUMMARY_WRITER_HH
#define EWOMS_ECL_SUMMARY_WRITER_HH

#include "ertwrappers.hh"
#include "eclwellmanager.hh"
#include "ecldeckunits.hh"

#include <ewoms/common/propertysystem.hh>

#include <opm/common/Valgrind.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/Deck/Section.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>

#if HAVE_ERT
#include <ert/ecl/ecl_sum.h>
#endif

#include <boost/algorithm/string.hpp>

#include <string>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(EnableEclSummaryOutput);
}

template <class TypeTag>
class EclSummaryWriter;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Implements writing of ECL summary files.
 *
 * i.e., well rates, bottom hole pressures, etc.
 */
template <class TypeTag>
class EclSummaryWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Ewoms::EclWellManager<TypeTag> WellManager;
    typedef Ewoms::ErtSummary<TypeTag> ErtSummary;


    struct ErtWellInfo {
        smspec_node_type* wbhpErtHandle;
        smspec_node_type* wthpErtHandle;
        smspec_node_type* wgorErtHandle;
        smspec_node_type* wwirErtHandle;
        smspec_node_type* wgirErtHandle;
        smspec_node_type* woirErtHandle;
        smspec_node_type* wwprErtHandle;
        smspec_node_type* wgprErtHandle;
        smspec_node_type* woprErtHandle;
        smspec_node_type* wwitErtHandle;
        smspec_node_type* wgitErtHandle;
        smspec_node_type* woitErtHandle;
        smspec_node_type* wwptErtHandle;
        smspec_node_type* wgptErtHandle;
        smspec_node_type* woptErtHandle;
    };

    static const unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static const unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static const unsigned oilPhaseIdx = FluidSystem::oilPhaseIdx;

public:
    EclSummaryWriter(const Simulator& simulator)
        : simulator_(simulator)
#if HAVE_ERT
        , ertSummary_(simulator)
#endif
    {
        const auto& deck = simulator.gridManager().deck();

        // populate the set of quantities to write
        if (deck.hasKeyword("ALL"))
            handleAllKeyword__();
        else
            addPresentSummaryKeywords_(deck);

        addVariables_(simulator.gridManager().eclState(), simulator.gridManager().schedule());
    }

    ~EclSummaryWriter()
    { }

    /*!
     * \brief Adds an entry to the summary file.
     */
    void write(const WellManager& wellsManager, bool isInitial = false)
    {
        unsigned reportIdx = simulator_.episodeIndex();
        Scalar t = simulator_.time();
        if (!isInitial) {
            t += simulator_.timeStepSize();
            reportIdx += 1;
        }

        ErtSummaryTimeStep<TypeTag> ertSumTimeStep(ertSummary_, t, reportIdx);

        typedef EclDeckUnits<TypeTag> DeckUnits;
        const DeckUnits& deckUnits = simulator_.problem().deckUnits();

        // add the well quantities
        for (unsigned wellIdx = 0; wellIdx < wellsManager.numWells(); ++wellIdx) {
            const auto& well = wellsManager.well(wellIdx);
            const auto& summaryInfo = ertWellInfo_.at(well->name());

            if (writeWbhp_()) {
                Scalar bhpPascal = well->bottomHolePressure();
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wbhpErtHandle),
                                   deckUnits.siToDeck(bhpPascal, DeckUnits::pressure));
            }

            if (writeWthp_()) {
                Scalar thpPascal = well->tubingHeadPressure();
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wthpErtHandle),
                                   deckUnits.siToDeck(thpPascal, DeckUnits::pressure));
            }

            if (writeWgor_()) {
                // since I'm usure what the gas-to-oil ratio exactly expresses, I just
                // assume "volume of gas at standard conditions divided by volume of oil
                // at standard conditions". Mass-based measures would be drastically
                // different. (As will be if imperial units are used where the volume of
                // gas is MCF and the volume of oil is bbl)
                Scalar gasRate = std::abs(well->surfaceRate(gasPhaseIdx));
                Scalar oilRate = std::abs(well->surfaceRate(oilPhaseIdx));

                Scalar gasToOilRatio = 0;
                if (std::abs(oilRate) > 1e-3)
                    gasToOilRatio = gasRate/oilRate;

                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wgorErtHandle),
                                   deckUnits.siToDeck(gasToOilRatio, DeckUnits::gasOilRatio));
            }

            //////////
            // injection surface rates
            if (writeWwir_()) {
                Scalar ratePerSecond = std::max<Scalar>(0.0, well->surfaceRate(waterPhaseIdx));
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wwirErtHandle),
                                   deckUnits.siToDeck(ratePerSecond, DeckUnits::liquidRate));
            }

            if (writeWgir_()) {
                Scalar ratePerSecond = std::max<Scalar>(0.0, well->surfaceRate(gasPhaseIdx));
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wgirErtHandle),
                                   deckUnits.siToDeck(ratePerSecond, DeckUnits::gasRate));
            }

            if (writeWoir_()) {
                Scalar ratePerSecond = std::max<Scalar>(0.0, well->surfaceRate(oilPhaseIdx));
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.woirErtHandle),
                                   deckUnits.siToDeck(ratePerSecond, DeckUnits::liquidRate));
            }
            //////////

            //////////
            // total injected surface volume
            if (writeWwit_()) {
                Scalar totalVolume = wellsManager.totalInjectedVolume(well->name(), waterPhaseIdx);
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wwitErtHandle),
                                   deckUnits.siToDeck(totalVolume, DeckUnits::liquidSurfaceVolume));
            }

            if (writeWgit_()) {
                Scalar totalVolume = wellsManager.totalInjectedVolume(well->name(), gasPhaseIdx);
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wgitErtHandle),
                                   deckUnits.siToDeck(totalVolume, DeckUnits::gasSurfaceVolume));
            }

            if (writeWoit_()) {
                Scalar totalVolume = wellsManager.totalInjectedVolume(well->name(), oilPhaseIdx);
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.woitErtHandle),
                                   deckUnits.siToDeck(totalVolume, DeckUnits::liquidSurfaceVolume));
            }
            //////////

            //////////
            // production surface rates
            if (writeWwpr_()) {
                Scalar ratePerSecond = std::max<Scalar>(0.0, -well->surfaceRate(waterPhaseIdx));
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wwprErtHandle),
                                   deckUnits.siToDeck(ratePerSecond, DeckUnits::liquidRate));
            }

            if (writeWgpr_()) {
                Scalar ratePerSecond = std::max<Scalar>(0.0, -well->surfaceRate(gasPhaseIdx));
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wgprErtHandle),
                                   deckUnits.siToDeck(ratePerSecond, DeckUnits::gasRate));
            }

            if (writeWopr_()) {
                Scalar ratePerSecond = std::max<Scalar>(0.0, -well->surfaceRate(oilPhaseIdx));
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.woprErtHandle),
                                   deckUnits.siToDeck(ratePerSecond, DeckUnits::liquidRate));
            }
            //////////

            //////////
            // total producted surface volume
            if (writeWwpt_()) {
                Scalar totalVolume = wellsManager.totalProducedVolume(well->name(), waterPhaseIdx);
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wwptErtHandle),
                                   deckUnits.siToDeck(totalVolume, DeckUnits::liquidSurfaceVolume));
            }

            if (writeWgpt_()) {
                Scalar totalVolume = wellsManager.totalProducedVolume(well->name(), gasPhaseIdx);
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.wgptErtHandle),
                                   deckUnits.siToDeck(totalVolume, DeckUnits::gasSurfaceVolume));
            }

            if (writeWopt_()) {
                Scalar totalVolume = wellsManager.totalProducedVolume(well->name(), oilPhaseIdx);
                ecl_sum_tstep_iset(ertSumTimeStep.ertHandle(),
                                   smspec_node_get_params_index(summaryInfo.woptErtHandle),
                                   deckUnits.siToDeck(totalVolume, DeckUnits::liquidSurfaceVolume));
            }
            //////////
        }

        // write the _complete_ summary file!
        ecl_sum_fwrite(ertSummary_.ertHandle());
    }

private:
    static bool enableEclSummaryOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclSummaryOutput); }

    bool writeWbhp_() const
    { return summaryKeywords_.count("WBHP") > 0; }

    bool writeWthp_() const
    { return summaryKeywords_.count("WTHP") > 0; }

    // gas-oil ratio
    bool writeWgor_() const
    { return summaryKeywords_.count("WGOR") > 0; }

    // write each well's current surface water injection rate
    bool writeWwir_() const
    { return summaryKeywords_.count("WWIR") > 0; }

    // write each well's current surface gas injection rate
    bool writeWgir_() const
    { return summaryKeywords_.count("WGIR") > 0; }

    // write each well's current surface oil injection rate
    bool writeWoir_() const
    { return summaryKeywords_.count("WOIR") > 0; }

    // write each well's current surface water production rate
    bool writeWwpr_() const
    { return summaryKeywords_.count("WWPR") > 0; }

    // write each well's current surface gas production rate
    bool writeWgpr_() const
    { return summaryKeywords_.count("WGPR") > 0; }

    // write each well's current surface oil production rate
    bool writeWopr_() const
    { return summaryKeywords_.count("WOPR") > 0; }

    // write each well's current surface water injection total
    bool writeWwit_() const
    { return summaryKeywords_.count("WWIT") > 0; }

    // write each well's current surface gas injection total
    bool writeWgit_() const
    { return summaryKeywords_.count("WGIT") > 0; }

    // write each well's current surface oil injection total
    bool writeWoit_() const
    { return summaryKeywords_.count("WOIT") > 0; }

    // write each well's current surface water production total
    bool writeWwpt_() const
    { return summaryKeywords_.count("WWPT") > 0; }

    // write each well's current surface gas production total
    bool writeWgpt_() const
    { return summaryKeywords_.count("WGPT") > 0; }

    // write each well's current surface oil production total
    bool writeWopt_() const
    { return summaryKeywords_.count("WOPT") > 0; }

    void addVariables_(const Opm::EclipseState& eclState, const Opm::Schedule& schedule)
    {
        const auto& wellsVector = schedule.getWells();
        for (size_t wellIdx = 0; wellIdx < wellsVector.size(); ++ wellIdx) {
            const auto& eclWell = wellsVector[wellIdx];
            auto& wellInfo = ertWellInfo_[eclWell->name()];

            // the bottom hole and tubing head pressure
            if (writeWbhp_())
                wellInfo.wbhpErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WBHP",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"BARSA",
                                    /*defaultValue=*/0.0);

            if (writeWthp_())
                wellInfo.wthpErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WTHP",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"BARSA",
                                    /*defaultValue=*/0.0);

            // the gas to oil rate
            if (writeWgor_())
                wellInfo.wgorErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WGOR",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"",
                                    /*defaultValue=*/0.0);

            // add injection variables
            if (writeWwir_())
                wellInfo.wwirErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WWIR",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWgir_())
                wellInfo.wgirErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WGIR",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWoir_())
                wellInfo.woirErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WOIR",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            // add production variables
            if (writeWwpr_())
                wellInfo.wwprErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WWPR",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWgpr_())
                wellInfo.wgprErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WGPR",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWopr_())
                wellInfo.woprErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WOPR",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            // add total injected volume variables
            if (writeWwit_())
                wellInfo.wwitErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WWIT",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWgit_())
                wellInfo.wgitErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WGIT",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWoit_())
                wellInfo.woitErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WOIT",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            // add total produced volume variables
            if (writeWwpt_())
                wellInfo.wwptErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WWPT",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWgpt_())
                wellInfo.wgptErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WGPT",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);

            if (writeWopt_())
                wellInfo.woptErtHandle =
                    ecl_sum_add_var(ertSummary_.ertHandle(),
                                    "WOPT",
                                    eclWell->name().c_str(),
                                    /*num=*/0,
                                    /*unit=*/"SM3/DAY",
                                    /*defaultValue=*/0.0);
        }
    }

    // add all quantities which are implied by the ALL summary keyword
    void handleAllKeyword__()
    {
        // these are the keywords implied by ALL which are documented by the Eclipse
        // reference manual
        std::string allKw[] = {
                "FOPR", "GOPR", "WOPR", "FOPT",
                "GOPT", "WOPT", "FOIR", "GOIR",
                "WOIR", "FOIT", "GOIT", "WOIT",
                "FWPR", "GWPR", "WWPR", "FWPT",
                "GWPT", "WWPT", "FWIR", "GWIR",
                "WWIR", "FWIT", "GWIT", "WWIT",
                "FGPR", "GGPR", "WGPR", "FGPT",
                "GGPT", "WGPT", "FGIR", "GGIR",
                "WGIR", "FGIT", "GGIT", "WGIT",
                "FVPR", "GVPR", "WVPR", "FVPT",
                "GVPT", "WVPT", "FVIR", "GVIR",
                "WVIR", "FVIT", "GVIT", "WVIT",
                "FWCT", "GWCT", "WWCT", "FGOR",
                "GGOR", "WGOR", "FWGR", "GWGR",
                "WWGR", "WBHP", "WTHP", "WPI",
                "FOIP", "FOIPL", "FOIPG", "FWIP",
                "FGIP", "FGIPL", "FGIPG", "FPR",
                "FAQR", "FAQRG", "AAQR", "AAQRG",
                "FAQT", "FAQTG", "AAQT", "AAQTG",
        };
        unsigned numSummaryKeywords = sizeof(allKw)/sizeof(allKw[0]);

        for (unsigned kwIdx = 0; kwIdx < numSummaryKeywords; ++kwIdx)
            summaryKeywords_.insert(allKw[kwIdx]);
    }

    // add all quantities which are present in the summary section of the deck
    void addPresentSummaryKeywords_(const Opm::Deck& deck)
    {
        Opm::Section summarySection(deck, "SUMMARY");

        if (summarySection.size() == 0)
            return;

        auto kwIt = summarySection.begin();
        auto kwEndIt = summarySection.end();
        // skip the first keyword as this is "SUMMARY". bug in opm-parser?
        ++kwIt;

        for (; kwIt != kwEndIt; ++kwIt)
            summaryKeywords_.insert((*kwIt).name());
    }

    const Simulator& simulator_;

    std::set<std::string> summaryKeywords_;
    std::map<std::string, ErtWellInfo> ertWellInfo_;

#if HAVE_ERT
    ErtSummary ertSummary_;
#endif
};
} // namespace Ewoms

#endif
