// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::PvsModel
 */
#ifndef EWOMS_PVS_MODEL_HH
#define EWOMS_PVS_MODEL_HH

#include "pvsproperties.hh"

#include <ewoms/models/modules/diffusion/vcfvdiffusionmodule.hh>
#include <ewoms/models/modules/energy/vcfvenergymodule.hh>
#include <ewoms/disc/vcfv/vcfvmodel.hh>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Ewoms {

/*!
 * \ingroup PvsModel
 *
 * \brief A generic compositional multi-phase model using primary-variable switching.
 *
 * This model assumes a flow of \f$M \geq 1\f$ fluid phases
 * \f$\alpha\f$, each of which is assumed to be a mixture \f$N \geq
 * M\f$ chemical species \f$\kappa\f$.
 *
 * By default, the standard multi-phase Darcy approach is used to determine
 * the velocity, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\mathbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
 * although the actual approach which is used can be specified via the
 * \c VelocityModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, VelocityModule, Ewoms::VcfvForchheimerVelocityModule<TypeTag>);
 * \endcode
 *
 * The core of the model is the conservation mass of each component by
 * means of the equation
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \mathrm{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha  \right\}
 * - q^\kappa = 0 \;.
 * \f]
 *
 * To close the system mathematically, \f$M\f$ model equations are
 * also required. This model uses the primary variable switching
 * assumptions, which are given by:
 * \f[
  0 \stackrel{!}{=}
  f_\alpha = \left\{
  \begin{array}{cl}
    S_\alpha & \quad \text{if phase }\alpha\text{ is not present} \\
    1 - \sum_\kappa x_\alpha^\kappa & \quad \text{else}
  \end{array}
  \right.
 * \f]
 *
 * To make this approach applicable, a pseudo primary variable
 * <em>phase presence</em> has to be introduced. Its purpose is to
 * specify for each phase whether it is present or not. It is a
 * <em>pseudo</em> primary variable because it is not directly considered when
 * linearizing the system in the Newton method, but after each Newton
 * iteration, it gets updated like the "real" primary variables.  The
 * following rules are used for this update procedure:
 *
 * <ul>

 * <li>If phase \f$\alpha\f$ is present according to the pseudo
 *     primary variable, but \f$S_\alpha < 0\f$ after the Newton
 *     update, consider the phase \f$\alpha\f$ disappeared for the
 *     next iteration and use the set of primary variables which
 *     correspond to the new phase presence.</li>

 * <li>If phase \f$\alpha\f$ is not present according to the pseudo
 *     primary variable, but the sum of the component mole fractions
 *     in the phase is larger than 1, i.e. \f$\sum_\kappa
 *     x_\alpha^\kappa > 1\f$, consider the phase \f$\alpha\f$ present
 *     in the the next iteration and update the set of primary
 *     variables to make it consistent with the new phase
 *     presence.</li>
 *
 * <li>In all other cases don't modify the phase presence for phase
 *     \f$\alpha\f$.</li>
 *
 * </ul>
 *
 * The model always requires \f$N\f$ primary variables, but their
 * interpretation is dependent on the phase presence:
 *
 * <ul>
 *
 * <li>The first primary variable is always interpreted as the
 *      pressure of the phase with the lowest index \f$PV_0 =
 *      p_0\f$.</li>
 *
 * <li>Then, \f$M - 1\f$ "switching primary variables" follow, which
 *     are interpreted depending in the presence of the first
 *     \f$M-1\f$ phases: If phase \f$\alpha\f$ is present, its
 *     saturation \f$S_\alpha = PV_i\f$ is used as primary variable;
 *     if it is not present, the mole fraction \f$PV_i =
 *     x_{\alpha^\star}^\alpha\f$ of the component with index
 *     \f$\alpha\f$ in the phase with the lowest index that is present
 *     \f$\alpha^\star\f$ is used instead.</li>
 *
 * <li>Finally, the mole fractions of the \f$N-M\f$ components with
 *     the largest index in the phase with the lowest index that is
 *     present \f$x_{\alpha^\star}^\kappa\f$ are used as primary
 *     variables.</li>
 *
 * </ul>
 *
 * This model is then discretized using a fully-coupled vertex
 * centered finite volume scheme as spatial and the implicit Euler
 * method as temporal discretization.
 */
template<class TypeTag>
class PvsModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef VcfvModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };


    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef VcfvEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Register all run-time parameters for the immiscible VCVF discretization.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VcfvVtkPhasePresenceModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkMultiPhaseModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkCompositionModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkTemperatureModule<TypeTag>::registerParameters();

        if (enableDiffusion)
            Ewoms::VcfvVtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Ewoms::VcfvVtkEnergyModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, int, PvsVerbosity, "The verbosity level of the primary variable switching model");
    }

    /*!
     * \copydoc VcfvModel::init
     */
    void init(Problem &problem)
    {
        verbosity_ = EWOMS_GET_PARAM(TypeTag, int, PvsVerbosity);
        numSwitched_ = 0;

        ParentType::init(problem);

        intrinsicPermeability_.resize(this->numDofs());
    }

    /*!
     * \copydoc VcfvModel::name
     */
    std::string name() const
    { return "pvs"; }

    /*!
     * \copydoc VcfvModel::primaryVarName
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
        else if (Indices::switch0Idx + numPhases - 1 <= pvIdx &&
                 pvIdx < Indices::switch0Idx + numComponents - 1)
            oss << "auxMoleFrac^" << FluidSystem::componentName(pvIdx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc VcfvModel::eqName
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
     * \copydoc VcfvModel::updateFailed
     */
    void updateFailed()
    {
        ParentType::updateFailed();
        numSwitched_ = 0;
    }

    /*!
     * \copydoc VcfvModel::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        referencePressure_ = this->solution(/*timeIdx=*/0)[/*vertexIdx=*/0][/*pvIdx=*/Indices::pressure0Idx];
    }

    /*!
     * \copydoc VcfvModel::updatePVWeights
     */
    void updatePVWeights(const ElementContext &elemCtx) const
    {
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
            const auto &K = elemCtx.volVars(scvIdx, /*timeIdx=*/0).intrinsicPermeability();

            int globalIdx = elemCtx.globalSpaceIndex(scvIdx, /*timeIdx=*/0);
            intrinsicPermeability_[globalIdx] = K[0][0];
        }
    }

    /*!
     * \copydoc VcfvModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalVertexIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;

        if (Indices::pressure0Idx == pvIdx) {
            // use a pressure gradient of 1e2 Pa/m for liquid water as
            // a reference
            Scalar KRef = intrinsicPermeability_[globalVertexIdx];
            static const Scalar muRef = 1e-3;
            static const Scalar pGradRef = 1e-2; // [Pa / m]
            Scalar r = std::pow(this->boxVolume(globalVertexIdx), 1.0/dimWorld);

            return std::max(10/referencePressure_, pGradRef * KRef/muRef / r);
        }

        if (Indices::switch0Idx <= pvIdx && pvIdx < Indices::switch0Idx + numPhases - 1) {
            int phaseIdx = pvIdx -Indices::switch0Idx;

            if (!this->solution(/*timeIdx=*/0)[globalVertexIdx].phaseIsPresent(phaseIdx))
                // for saturations, the weight is always 1
                return 1;

            // for saturations, the PvsMoleSaturationsBaseWeight
            // property determines the weight
            return GET_PROP_VALUE(TypeTag, PvsSaturationsBaseWeight);
        }

        // for mole fractions, the PvsMoleFractionsBaseWeight
        // property determines the weight
        return GET_PROP_VALUE(TypeTag, PvsMoleFractionsBaseWeight);
    }

    /*!
     * \copydoc VcfvModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalVertexIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;

        int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numComponents);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    /*!
     * \copydoc VcfvModel::advanceTimeLevel
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
     * \copydoc VcfvModel::serializeEntity
     */
    void serializeEntity(std::ostream &outstream, const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outstream, vert);

        int vertIdx = this->dofMapper().map(vert);
        if (!outstream.good())
            OPM_THROW(std::runtime_error, "Could not serialize vertex " << vertIdx);

        outstream << this->solution(/*timeIdx=*/0)[vertIdx].phasePresence() << " ";
    }

    /*!
     * \copydoc VcfvModel::deserializeEntity
     */
    void deserializeEntity(std::istream &instream, const Vertex &vertex)
    {
        // read primary variables
        ParentType::deserializeEntity(instream, vertex);

        // read phase presence
        int vertIdx = this->dofMapper().map(vertex);
        if (!instream.good())
            OPM_THROW(std::runtime_error,
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
            this->problem_().newtonMethod().endIterMsg() << ", num switched=" << numSwitched_;
    }

private:
    friend class VcfvModel<TypeTag>;

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
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkPhasePresenceModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkTemperatureModule<TypeTag>(this->problem_()));
        if (enableDiffusion)
            this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkDiffusionModule<TypeTag>(this->problem_()));
        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkEnergyModule<TypeTag>(this->problem_()));
    }

    mutable Scalar referencePressure_;
    mutable std::vector<Scalar> intrinsicPermeability_;

    // number of switches of the phase state in the last Newton
    // iteration
    int numSwitched_;

    // verbosity of the model
    int verbosity_;
};

}

#include "pvspropertydefaults.hh"

#endif
