// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Holger Class                                 *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 *
 * \brief Adaption of the BOX scheme to the three-phase three-component flow model.
 *
 * The model is designed for simulating three fluid phases with water, gas, and
 * a liquid contaminant (NAPL - non-aqueous phase liquid)
 */
#ifndef DUMUX_3P3C_MODEL_HH
#define DUMUX_3P3C_MODEL_HH

#include "3p3cproperties.hh"
#include "3p3clocalresidual.hh"
#include "3p3cproblem.hh"

#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCModel
 * \brief Adaption of the BOX scheme to the three-phase three-component flow model.
 *
 * This model implements three-phase three-component flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$ each composed of up to three components
 * \f$\kappa \in \{ water, air, contaminant \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one transport equation for each component is obtained as
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_{\text{mol}, \alpha} x_\alpha^\kappa
 S_\alpha )}{\partial t}
 - \sum\limits_\alpha \text{div} \left\{ \frac{k_{r\alpha}}{\mu_\alpha}
 \varrho_{\text{mol}, \alpha} x_\alpha^\kappa \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\text{mass}, \alpha} \mbox{\bf g}) \right\}
 \nonumber \\
 \nonumber \\
 && - \sum\limits_\alpha \text{div} \left\{ D_{pm}^\kappa \varrho_{\text{mol},
 \alpha } \text{grad}\, x_\alpha^\kappa \right\}
 - q^\kappa = 0 \qquad \forall \kappa , \; \forall \alpha
 \f}
 *
 * Note that these balance equations are molar.
 *
 * The equations are discretized using a fully-coupled vertex
 * centered finite volume (BOX) scheme as spatial scheme and
 * the implicit Euler method as temporal discretization.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations and
 * \f$x^w_\alpha + x^a_\alpha + x^c_\alpha = 1\f$ for the mole fractions.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are dependent on the locally present fluid phases
 * An adaptive primary variable switch is included. The phase state is stored for all nodes
 * of the system. The following cases can be distinguished:
 * <ul>
 *  <li> All three phases are present: Primary variables are two saturations \f$(S_w\f$ and \f$S_n)\f$, and a pressure, in this case \f$p_g\f$. </li>
 *  <li> Only the water phase is present: Primary variables are now the mole fractions of air and contaminant in the water phase \f$(x_w^a\f$ and \f$x_w^c)\f$, as well as the gas pressure, which is, of course, in a case where only the water phase is present, just the same as the water pressure. </li>
 *  <li> Gas and NAPL phases are present: Primary variables \f$(S_n\f$, \f$x_g^w\f$, \f$p_g)\f$. </li>
 *  <li> Water and NAPL phases are present: Primary variables \f$(S_n\f$, \f$x_w^a\f$, \f$p_g)\f$. </li>
 *  <li> Only gas phase is present: Primary variables \f$(x_g^w\f$, \f$x_g^c\f$, \f$p_g)\f$. </li>
 *  <li> Water and gas phases are present: Primary variables \f$(S_w\f$, \f$x_w^g\f$, \f$p_g)\f$. </li>
 * </ul>
 */
template<class TypeTag>
class ThreePThreeCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        pressure0Idx = Indices::pressure0Idx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        cCompIdx = Indices::cCompIdx,
        aCompIdx = Indices::aCompIdx,

        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly

    };

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        verbosity_ = GET_PARAM(TypeTag, int, ThreePThreeCVerbosity);
        numSwitched_ = 0;
    }

    /*!
     * \brief Returns a string with the model's human-readable name
     */
    std::string name() const
    { return "3p3c"; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::ostringstream oss;
        if (pvIdx == Indices::pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (Indices::switch1Idx <= pvIdx && pvIdx < Indices::switch1Idx + numPhases - 1)
            oss << "switch_" << pvIdx - Indices::switch1Idx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    std::string eqName(int eqIdx) const
    { 
        std::ostringstream oss;
        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents) {
            int compIdx = eqIdx - Indices::conti0EqIdx;
            oss << "continuity^" << FluidSystem::componentName(compIdx);
        }
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \param dest Contains the storage of each component for one phase
     * \param phaseIdx The phase index
     */
    void globalPhaseStorage(EqVector &dest, int phaseIdx)
    {
        dest = 0;

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);

            this->localResidual().addPhaseStorage(dest, elemCtx, /*timeIdx=*/0, phaseIdx);
        };

        this->gridView_().comm().sum(dest);
    }

    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();

        auto &sol = this->solution(/*timeIdx=*/0);
        for (int i = 0; i < sol.size(); ++i)
            sol[i].setSwitched(false);
        numSwitched_ = 0;
    };

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global vertex index
     * \param pvIdx The primary variable index
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pressure0Idx == pvIdx)
            return std::min(1.0/this->solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx], 1.0);
        return 1;
    }

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();
        numSwitched_ = 0;
    }

    /*!
     * \brief Return true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    { return numSwitched_ > 0; }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one vertex for the restart file
     * \param vert The vertex
     */
    void serializeEntity(std::ostream &outStream, const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vert);

        int vertIdx = this->dofMapper().map(vert);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vertIdx);

        outStream << this->solution(/*timeIdx=*/0)[vertIdx].phasePresence() << " ";
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param vert The vertex
     */
    void deserializeEntity(std::istream &inStream, const Vertex &vert)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, vert);

        // read phase presence
        int vertIdx = this->dofMapper().map(vert);
        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize vertex " << vertIdx);

        int tmp;
        inStream >> tmp;
        this->solution(/*timeIdx=*/0)[vertIdx].setPhasePresence(tmp);
        this->solution(/*timeIdx=*/1)[vertIdx].setPhasePresence(tmp);
    }

    /*!
     * \brief Do the primary variable switching after a Newton iteration.
     */
    void switchPrimaryVars_()
    {
        numSwitched_ = 0;

        std::vector<bool> visited(this->numDofs(), false);
        ElementContext elemCtx(this->problem_());

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            bool fvElemGeomUpdated = false;
            int numScv = elemIt->template count<dim>();
            for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, scvIdx, dim);

                if (visited[globalIdx])
                    continue;
                visited[globalIdx] = true;

                if (!fvElemGeomUpdated) {
                    fvElemGeomUpdated = true;
                    elemCtx.updateFVElemGeom(*elemIt);
                }

                // compute the volume variables of the current
                // sub-control volume
                auto &priVars = this->solution(/*timeIdx=*/0)[globalIdx];
                elemCtx.updateScvVars(priVars,
                                       scvIdx,
                                       /*timeIdx=*/0);

                // evaluate primary variable switch
                short oldPhasePresence = priVars.phasePresence();

                // set the primary variables and the new phase state
                // from the current fluid state
                primaryVarSwitch_(elemCtx, scvIdx, /*timeIdx=*/0);

                if (oldPhasePresence != priVars.phasePresence())
                {
                    priVars.setSwitched(true);
                    this->jacobianAssembler().markVertexRed(globalIdx);
                    ++ numSwitched_;
                }
                else
                    priVars.setSwitched(false);
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        numSwitched_ = this->gridView_().comm().sum(numSwitched_);

        if (verbosity_ > 0)
            this->problem_().newtonController().endIterMsg() << ", num switched=" << numSwitched_;
    }

protected:
    friend class BoxModel<TypeTag>;

    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    void primaryVarSwitch_(const ElementContext &elemCtx,
                           int spaceIdx,
                           int timeIdx)
    {
        int globalIdx = elemCtx.globalSpaceIndex(spaceIdx, timeIdx);
        const auto &globalPos = elemCtx.pos(spaceIdx, timeIdx);
        const auto &volVars = elemCtx.volVars(spaceIdx, timeIdx);
        const auto &fs = volVars.fluidState();
        auto &priVars = this->solution(timeIdx)[globalIdx];

        // evaluate primary variable switch
        int phasePresence = priVars.phasePresence();
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == threePhases)
        {
            Scalar Smin = 0;
            if (priVars.wasSwitched())
                Smin = -0.01;

            if (fs.saturation(gPhaseIdx) <= Smin)
            {
                // gas phase disappears
                if (verbosity_ > 1)
                    std::cout << "Gas phase disappears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", Sg: "
                              << fs.saturation(gPhaseIdx) << std::endl;
                newPhasePresence = wnPhaseOnly;

                priVars[switch1Idx] = fs.moleFraction(wPhaseIdx, aCompIdx);
            }
            else if (fs.saturation(wPhaseIdx) <= Smin)
            {
                // water phase disappears
                if (verbosity_ > 1)
                    std::cout << "Water phase disappears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", Sw: "
                              << fs.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = gnPhaseOnly;

                priVars[switch1Idx] = fs.moleFraction(gPhaseIdx, wCompIdx);
            }
            else if (fs.saturation(nPhaseIdx) <= Smin)
            {
                // NAPL phase disappears
                if (verbosity_ > 1)
                    std::cout << "NAPL phase disappears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", Sn: "
                              << fs.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wgPhaseOnly;

                priVars[switch2Idx] = fs.moleFraction(gPhaseIdx, cCompIdx);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            bool gasFlag = false;
            bool NAPLFlag = false;
            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = fs.moleFraction(gPhaseIdx, wCompIdx);
            Scalar xag = fs.moleFraction(gPhaseIdx, aCompIdx);
            Scalar xcg = fs.moleFraction(gPhaseIdx, cCompIdx);
            /* take care:
               for xag in case wPhaseOnly we compute xag=henry_air*xaw
               for xwg in case wPhaseOnly we compute xwg=pwsat
               for xcg in case wPhaseOnly we compute xcg=henry_NAPL*xcw
            */

            Scalar xgMax = 1.0;
            if (priVars.wasSwitched())
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xwg + xag + xcg > xgMax)
            {
                // gas phase appears
                if (verbosity_ > 1)
                    std::cout << "gas phase appears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", xwg + xag + xcg: "
                              << xwg + xag + xcg << std::endl;
                gasFlag = true;
            }

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnc = fs.moleFraction(nPhaseIdx, cCompIdx);
            Scalar xna = fs.moleFraction(nPhaseIdx, aCompIdx);
            Scalar xnw = fs.moleFraction(nPhaseIdx, wCompIdx);
            /* take care:
               for xnc in case wPhaseOnly we compute xnc=henry_mesitylene*xcw,
               where a hypothetical gas pressure is assumed for the Henry
               xwn is set to NULL  (all NAPL phase is dirty)
               xan is set to NULL  (all NAPL phase is dirty)
            */

            Scalar xnMax = 1.0;
            if (priVars.wasSwitched())
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fractions would be larger than
            // 100%, NAPL phase appears
            if (xnc + xna + xnw > xnMax)
            {
                // NAPL phase appears
                if (verbosity_ > 1)
                    std::cout << "NAPL phase appears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", xnc: "
                              << xnc << std::endl;
                NAPLFlag = true;
            }

            if (gasFlag && !NAPLFlag)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = 1 - 1e-6;
                priVars[switch2Idx] = 1e-6;
            }
            else if (gasFlag && NAPLFlag)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = 1.0 - 1e-6;
                priVars[switch2Idx] = 1e-6;
            }
            else if (!gasFlag && NAPLFlag)
            {
                newPhasePresence = wnPhaseOnly;
                priVars[switch1Idx] = fs.moleFraction(wPhaseIdx, aCompIdx);
                priVars[switch2Idx] = 1e-6;
            }
        }
        else if (phasePresence == gnPhaseOnly)
        {
            bool NAPLFlag = false;
            bool waterFlag = false;

            Scalar Smin = 0.0;
            if (priVars.wasSwitched())
                Smin = -0.01;

            if (fs.saturation(nPhaseIdx) <= Smin)
            {
                // NAPL phase disappears
                if (verbosity_ > 1)
                    std::cout << "NAPL phase disappears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", Sn: "
                              << fs.saturation(nPhaseIdx) << std::endl;
                NAPLFlag = true;
            }

            // calculate fractions of the hypothetical water phase
            Scalar xww = fs.moleFraction(wPhaseIdx, wCompIdx);
            Scalar xwa = fs.moleFraction(wPhaseIdx, aCompIdx);
            Scalar xwc = fs.moleFraction(wPhaseIdx, cCompIdx);

            Scalar xwMax = 1.0;
            if (priVars.wasSwitched())
                xwMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xww + xwa + xwc > xwMax)
            {
                // water phase appears
                if (verbosity_ > 1)
                    std::cout << "water phase appears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                              << xww << std::endl;
                waterFlag = true;
            }

            if (waterFlag && !NAPLFlag)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = 0.0001;
                priVars[switch2Idx] = fs.saturation(nPhaseIdx);
            }
            else if (waterFlag && NAPLFlag)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = 1e-6;
                priVars[switch2Idx] = fs.moleFraction(gPhaseIdx, cCompIdx);
            }
            else if (!waterFlag && NAPLFlag)
            {
                newPhasePresence = gPhaseOnly;
                priVars[switch1Idx] = fs.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = fs.moleFraction(gPhaseIdx, cCompIdx);
            }
        }
        else if (phasePresence == wnPhaseOnly)
        {
            bool NAPLFlag = false;
            bool gasFlag = false;

            Scalar Smin = 0.0;
            if (priVars.wasSwitched())
                Smin = -0.01;

            if (fs.saturation(nPhaseIdx) <= Smin)
            {
                // NAPL phase disappears
                if (verbosity_ > 1)
                    std::cout << "NAPL phase disappears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", Sn: "
                              << fs.saturation(nPhaseIdx) << std::endl;
                NAPLFlag = true;
            }

            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = fs.moleFraction(gPhaseIdx, wCompIdx);
            Scalar xag = fs.moleFraction(gPhaseIdx, aCompIdx);
            Scalar xcg = fs.moleFraction(gPhaseIdx, cCompIdx);
            /* take care:
               for xag in case wPhaseOnly we compute xag=henry_air*xaw
               for xwg in case wPhaseOnly we compute xwg=pwsat
               for xcg in case wPhaseOnly we compute xcg=henry_NAPL*xcw
            */
            Scalar xgMax = 1.0;
            if (priVars.wasSwitched())
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xwg + xag + xcg > xgMax)
            {
                // gas phase appears
                if (verbosity_ > 1)
                    std::cout << "gas phase appears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", xwg + xag + xcg: "
                              << xwg + xag + xcg << std::endl;
                gasFlag = true;
            }

            if (gasFlag && !NAPLFlag)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = fs.saturation(wPhaseIdx);
                priVars[switch2Idx] = /*Sn=*/1e-6;
            }
            else if (gasFlag && NAPLFlag)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = fs.saturation(wPhaseIdx);
                priVars[switch2Idx] = fs.moleFraction(gPhaseIdx, cCompIdx);
            }
            else if (!gasFlag && NAPLFlag)
            {
                newPhasePresence = wPhaseOnly;
                priVars[switch1Idx] = fs.moleFraction(wPhaseIdx, aCompIdx);
                priVars[switch2Idx] = fs.moleFraction(wPhaseIdx, cCompIdx);
            }
        }
        else if (phasePresence == gPhaseOnly)
        {
            bool NAPLFlag = false;
            bool waterFlag = false;

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnc = fs.moleFraction(nPhaseIdx, cCompIdx);
            Scalar xnw = fs.moleFraction(nPhaseIdx, wCompIdx);
            Scalar xna = fs.moleFraction(nPhaseIdx, aCompIdx);
            /*
              take care:, xnc, if no NAPL phase is there, take xnc=xcg*pg/pcsat
              if this is larger than 1, then NAPL appears
            */

            Scalar xnMax = 1.0;
            if (priVars.wasSwitched())
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fraction would be larger than
            // 100%, NAPL phase appears
            if (xnc + xnw + xna > xnMax)
            {
                // NAPL phase appears
                if (verbosity_ > 1)
                    std::cout << "NAPL phase appears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", xnc: "
                              << xnc << std::endl;
                NAPLFlag = true;
            }
            // calculate fractions of the hypothetical water phase
            Scalar xww = fs.moleFraction(wPhaseIdx, wCompIdx);
            Scalar xwa = fs.moleFraction(wPhaseIdx, aCompIdx);
            Scalar xwc = fs.moleFraction(wPhaseIdx, cCompIdx);
            /*
              take care:, xww, if no water is present, then take xww=xwg*pg/pwsat .
              If this is larger than 1, then water appears
            */
            Scalar xwMax = 1.0;
            if (priVars.wasSwitched())
                xwMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xww + xwa + xwc > xwMax)
            {
                // water phase appears
                if (verbosity_ > 1)
                    std::cout << "water phase appears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                              << xww << std::endl;
                waterFlag = true;
            }
            if (waterFlag && !NAPLFlag)
            {
                newPhasePresence = wgPhaseOnly;
                priVars[switch1Idx] = 1e-6;
                priVars[switch2Idx] = fs.moleFraction(gPhaseIdx, cCompIdx);
            }
            else if (waterFlag && NAPLFlag)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = 1e-6;
                priVars[switch2Idx] = 1e-6;
            }
            else if (!waterFlag && NAPLFlag)
            {
                newPhasePresence = gnPhaseOnly;
                priVars[switch1Idx] = fs.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = 1e-6;
            }
        }
        else if (phasePresence == wgPhaseOnly)
        {
            bool NAPLFlag = false;
            bool gasFlag = false;
            bool waterFlag = false;

            // get the fractions in the hypothetical NAPL phase
            Scalar xnc = fs.moleFraction(nPhaseIdx, cCompIdx);
            Scalar xna = fs.moleFraction(nPhaseIdx, aCompIdx);
            Scalar xnw = fs.moleFraction(nPhaseIdx, wCompIdx);

            // take care: if the NAPL phase is not present, take
            // xnc=xcg*pg/pcsat if this is larger than 1, then NAPL
            // appears
            Scalar xnMax = 1.0;
            if (priVars.wasSwitched())
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fraction would be larger than
            // 100%, NAPL phase appears
            if (xnc + xna + xnw > xnMax)
            {
                if (verbosity_ > 1)
                    // NAPL phase appears
                    std::cout << "NAPL phase appears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", xnc: "
                              << xnc << std::endl;
                NAPLFlag = true;
            }

            Scalar Smin = -1.e-6;
            if (priVars.wasSwitched())
                Smin = -0.01;

            if (fs.saturation(gPhaseIdx) <= Smin)
            {
                // gas phase disappears
                if (verbosity_ > 1)
                    std::cout << "Gas phase disappears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", Sg: "
                              << fs.saturation(gPhaseIdx) << std::endl;
                gasFlag = true;
            }

            Smin = 0.0;
            if (priVars.wasSwitched())
                Smin = -0.01;

            if (fs.saturation(wPhaseIdx) <= Smin)
            {
                // water phase disappears
                if (verbosity_ > 1)
                    std::cout << "Water phase disappears at vertex " << globalIdx
                              << ", coordinates: " << globalPos << ", Sw: "
                              << fs.saturation(wPhaseIdx) << std::endl;
                waterFlag = true;
            }

            if (!gasFlag && NAPLFlag && waterFlag)
            {
                newPhasePresence = gnPhaseOnly;
                priVars[switch1Idx] = fs.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = 1e-6;
            }
            else if (!gasFlag && NAPLFlag && !waterFlag)
            {
                newPhasePresence = threePhases;
                priVars[switch1Idx] = fs.saturation(wPhaseIdx);
                priVars[switch2Idx] = 0.0;
            }
            else if (gasFlag && !NAPLFlag && !waterFlag)
            {
                newPhasePresence = wPhaseOnly;
                priVars[switch1Idx] = fs.moleFraction(wPhaseIdx, aCompIdx);
                priVars[switch2Idx] = fs.moleFraction(wPhaseIdx, cCompIdx);
            }
            else if (!gasFlag && !NAPLFlag && waterFlag)
            {
                newPhasePresence = gPhaseOnly;
                priVars[switch1Idx] = fs.moleFraction(gPhaseIdx, wCompIdx);
                priVars[switch2Idx] = fs.moleFraction(gPhaseIdx, cCompIdx);
            }
        }

        priVars.setPhasePresence(newPhasePresence);
    }

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules meaninful for the model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkPhasePresenceModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
    };

    // number of switches of the phase state in the last Newton
    // iteration
    int numSwitched_;
    
    // verbosity of the model
    int verbosity_;
};

}

#include "3p3cpropertydefaults.hh"

#endif
