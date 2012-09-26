// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2008-2012 by Klaus Mosthaf                                *
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
 * \copydoc Dumux::PvsModel
 */
#ifndef DUMUX_PVS_MODEL_HH
#define DUMUX_PVS_MODEL_HH

#include "pvsproperties.hh"
#include "pvslocalresidual.hh"

#include <dumux/boxmodels/modules/energy/boxmultiphaseenergymodule.hh>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Dumux {

/*!
 * \ingroup PvsModel
 *
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \sum_\alpha \text{div} \left\{{\bf D}_{\alpha, pm}^\kappa \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * PvsIndices::pWsN or PvsIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mass fraction of, e.g., air in the wetting phase \f$X^a_w\f$ is used,
 *      as long as the maximum mass fraction is not exceeded \f$(X^a_w<X^a_{w,max})\f$</li>
 *  <li> Only non-wetting phase is present: The mass fraction of, e.g., water in the non-wetting phase, \f$X^w_n\f$, is used,
 *      as long as the maximum mass fraction is not exceeded \f$(X^w_n<X^w_{n,max})\f$</li>
 * </ul>
 *
 * \todo implement/re-enable molecular diffusion
 */
template<class TypeTag>
class PvsModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef BoxMultiPhaseEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Register all run-time parameters for the immiscible box model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        
        // register runtime parameters of the VTK output modules
        Dumux::BoxVtkPhasePresenceModule<TypeTag>::registerParameters();
        Dumux::BoxVtkMultiPhaseModule<TypeTag>::registerParameters();
        Dumux::BoxVtkCompositionModule<TypeTag>::registerParameters();
        Dumux::BoxVtkTemperatureModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Dumux::BoxVtkEnergyModule<TypeTag>::registerParameters();

        REGISTER_PARAM(TypeTag, int, PvsVerbosity, "The verbosity level of the primary variable switching model");
    }

    /*!
     * \copydoc BoxModel::init
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        verbosity_ = GET_PARAM(TypeTag, int, PvsVerbosity);
        numSwitched_ = 0;
    }

    /*!
     * \copydoc BoxModel::name
     */
    std::string name() const
    { return "pvs"; }

    /*!
     * \copydoc BoxModel::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;
        if (pvIdx == Indices::pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (Indices::switch0Idx <= pvIdx && pvIdx < Indices::switch0Idx + numPhases - 1)
            oss << "switch_" << pvIdx - Indices::switch0Idx;
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc BoxModel::eqName
     */
    std::string eqName(int eqIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::eqName(eqIdx)).empty())
            return s;

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
     * \copydoc ImmiscibleModel::globalPhaseStorage
     */
    void globalPhaseStorage(EqVector &storage, int phaseIdx)
    {
        storage = 0;
        EqVector tmp;

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue; // ignore ghost and overlap elements

            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);

            const auto &fvElemGeom = elemCtx.fvElemGeom(/*timeIdx=*/0);
            
            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
                tmp = 0;
                this->localResidual().addPhaseStorage(tmp, 
                                                      elemCtx,
                                                      scvIdx,
                                                      /*timeIdx=*/0,
                                                      phaseIdx);
                tmp *= 
                    fvElemGeom.subContVol[scvIdx].volume
                    * elemCtx.volVars(scvIdx, /*timeIdx=*/0).extrusionFactor();
                storage += tmp;
            }
        };

        storage = this->gridView_().comm().sum(storage);
    }

    /*!
     * \copydoc BoxModel::updateFailed
     */
    void updateFailed()
    {
        ParentType::updateFailed();
        numSwitched_ = 0;
    }

    /*!
     * \copydoc BoxModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalVertexIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;

        if (Indices::pressure0Idx == pvIdx)
            return std::min(1.0/this->solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx], 1.0);
        return 1; // saturations and mole fractions have a weight of 1
    }

    /*!
     * \copydoc BoxModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalVertexIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;

        int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    /*!
     * \copydoc BoxModel::advanceTimeLevel
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
     * \copydoc BoxModel::serializeEntity
     */
    void serializeEntity(std::ostream &outstream, const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outstream, vert);

        int vertIdx = this->dofMapper().map(vert);
        if (!outstream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vertIdx);

        outstream << this->solution(/*timeIdx=*/0)[vertIdx].phasePresence() << " ";
    }

    /*!
     * \copydoc BoxModel::deserializeEntity
     */
    void deserializeEntity(std::istream &instream, const Vertex &vertex)
    {
        // read primary variables
        ParentType::deserializeEntity(instream, vertex);

        // read phase presence
        int vertIdx = this->dofMapper().map(vertex);
        if (!instream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize vertex " << vertIdx);

        int tmp;
        instream >> tmp;
        this->solution(/*timeIdx=*/0)[vertIdx].setPhasePresence(tmp);
        this->solution(/*timeIdx=*/1)[vertIdx].setPhasePresence(tmp);
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
                const VolumeVariables &volVars = elemCtx.volVars(scvIdx, /*timeIdx=*/0);

                // evaluate primary variable switch
                short oldPhasePresence = priVars.phasePresence();

                // set the primary variables and the new phase state
                // from the current fluid state
                priVars.assignNaive(volVars.fluidState());

                if (oldPhasePresence != priVars.phasePresence())
                {
                    if (verbosity_ > 1)
                        printSwitchedPhases_(elemCtx,
                                             scvIdx,
                                             volVars.fluidState(),
                                             oldPhasePresence,
                                             priVars);
                    this->jacobianAssembler().markVertexRed(globalIdx);
                    ++ numSwitched_;
                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        numSwitched_ = this->gridView_().comm().sum(numSwitched_);

        if (verbosity_ > 0)
            this->problem_().newtonController().endIterMsg() << ", num switched=" << numSwitched_;
    }

private:
    friend class BoxModel<TypeTag>;

    template <class FluidState>
    void printSwitchedPhases_(const ElementContext &elemCtx,
                              int scvIdx,
                              const FluidState &fs,
                              int oldPhasePresence,
                              const PrimaryVariables &newPv) const
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            bool oldPhasePresent = (oldPhasePresence & (1 << phaseIdx)) > 0; 
            bool newPhasePresent = newPv.phaseIsPresent(phaseIdx); 
            if (oldPhasePresent == newPhasePresent)
                continue;
            
            const auto &pos = elemCtx.pos(scvIdx, /*timeIdx=*/0);
            if (oldPhasePresent && !newPhasePresent) {
                std::cout << "'" << FluidSystem::phaseName(phaseIdx)
                          << "' phase disappears at position " << pos 
                          << ". saturation=" << fs.saturation(phaseIdx);
            }
            else {
                Scalar sumx = 0;
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    sumx += fs.moleFraction(phaseIdx, compIdx);

                std::cout << "'" << FluidSystem::phaseName(phaseIdx)
                          << "' phase appears at position " << pos 
                          << " sum x = " << sumx;
            }
        };
        
        std::cout << ", new primary variables: ";
        newPv.print();
        std::cout << "\n";
    }

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules meaninful for the model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkPhasePresenceModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Dumux::BoxVtkEnergyModule<TypeTag>(this->problem_()));
    }

    // number of switches of the phase state in the last Newton
    // iteration
    int numSwitched_;
    
    // verbosity of the model
    int verbosity_;
};

}

#include "pvspropertydefaults.hh"

#endif
