// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#include "immisciblelocalresidual.hh"

#include <ewoms/models/modules/energymodule.hh>

#include <dune/common/unused.hh>

#include <sstream>
#include <string>

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
class ImmiscibleModel : public Ewoms::MultiPhaseBaseModel<TypeTag>
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
     * \copydoc FvBaseDiscretization::init
     */
    void init()
    {
        ParentType::init();

        intrinsicPermeability_.resize(this->numDof());
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    const char *name() const
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
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::phaseIdxParam
     */
    void globalPhaseStorage(EqVector &storage, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        storage = 0;
        EqVector tmp;

        ElementContext elemCtx(this->problem_);
        ElementIterator elemIt = this->gridView_.template begin<0>();
        const ElementIterator elemEndIt = this->gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue; // ignore ghost and overlap elements

            elemCtx.updateStencil(*elemIt);
            elemCtx.updateVolVars(/*timeIdx=*/0);

            const auto &stencil = elemCtx.stencil(/*timeIdx=*/0);

            for (int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
                const auto &scv = stencil.subControlVolume(dofIdx);
                const auto &volVars = elemCtx.volVars(dofIdx, /*timeIdx=*/0);

                tmp = 0;
                this->localResidual().addPhaseStorage(tmp,
                                                      elemCtx,
                                                      dofIdx,
                                                      /*timeIdx=*/0,
                                                      phaseIdx);
                tmp *= scv.volume()*volVars.extrusionFactor();
                storage += tmp;
            }
        };

        storage = this->gridView_.comm().sum(storage);
    }

    /*!
     * \copydoc FvBaseDiscretization::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        referencePressure_ = this->solution(
            /*timeIdx=*/0)[/*vertexIdx=*/0][/*pvIdx=*/Indices::pressure0Idx];
    }

    /*!
     * \copydoc FvBaseDiscretization::updatePVWeights
     */
    void updatePVWeights(const ElementContext &elemCtx) const
    {
        for (int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
            int globalIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            const auto &K = elemCtx.volVars(dofIdx, /*timeIdx=*/0).intrinsicPermeability();
            intrinsicPermeability_[globalIdx] = K[0][0];
        }
    }

    /*!
     * \copydetails FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalDofIdx, int pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(asImp_(), globalDofIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;
        if (Indices::pressure0Idx == pvIdx) {
            // use a pressure gradient of 1e2 Pa/m for liquid water as
            // a reference
            Scalar KRef = intrinsicPermeability_[globalDofIdx];
            static const Scalar muRef = 1e-3;
            static const Scalar pGradRef = 1e-2; // [Pa / m]
            Scalar r = std::pow(this->boxVolume(globalDofIdx), 1.0/dimWorld);

            return std::max(10 / referencePressure_, pGradRef * KRef / muRef / r);
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

        DUNE_UNUSED int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);

        // make all kg equal
        return 1.0;
    }

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Ewoms::VtkEnergyModule<TypeTag>(this->problem_));
    }

private:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    mutable Scalar referencePressure_;
    mutable std::vector<Scalar> intrinsicPermeability_;
};
} // namespace Ewoms

#include "immisciblepropertydefaults.hh"

#endif
