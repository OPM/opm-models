// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::ImmiscibleModel
 */
#ifndef EWOMS_IMMISCIBLE_MODEL_HH
#define EWOMS_IMMISCIBLE_MODEL_HH

#include "immiscibleproperties.hh"
#include "immisciblelocalresidual.hh"

#include <ewoms/models/modules/energy/vcfvenergymodule.hh>

#include <sstream>
#include <string>

namespace Ewoms {
/*!
 * \ingroup ImmiscibleVcfvModel
 * \brief A fully-implicit multi-phase flow model which assumes
 *        immiscibility of the phases.
 *
 * This model implements multi-phase flow of \f$M > 0\f$ immiscible
 * fluids \f$\alpha\f$. By default, the standard multi-phase Darcy
 * approach is used to determine the velocity, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
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
 * \frac{\partial\;\phi S_\alpha \rho_\alpha }{\partial t}
 * - \text{div} \left\{ \rho_\alpha \mathbf{v}_\alpha  \right\}
 * - q_\alpha = 0 \;.
 * \f]
 *
 * These equations are discretized by a fully-coupled vertex centered
 * finite volume scheme as spatial and the implicit Euler method as
 * time discretization.
 *
 * The model uses the following primary variables:
 * - The pressure \f$p_0\f$ in Pascal of the phase with the lowest index
 * - The saturations \f$S_\alpha\f$ of the \f$M - 1\f$ phases that
 *   exhibit the lowest indices
 * - The absolute temperature \f$T\f$ in Kelvin if energy is conserved
 *   via the energy equation
 */
template<class TypeTag >
class ImmiscibleModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef VcfvModel<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { numComponents = FluidSystem::numComponents };

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef VcfvEnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Register all run-time parameters for the immiscible VCVF discretization.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VcfvVtkMultiPhaseModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkTemperatureModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Ewoms::VcfvVtkEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc VcfvModel::name
     */
    const char *name() const
    { return "immiscible"; }

    /*!
     * \copydoc VcfvModel::primaryVarName
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
        else if (Indices::saturation0Idx <= pvIdx && pvIdx < Indices::saturation0Idx + numPhases - 1) {
            int phaseIdx = pvIdx - Indices::saturation0Idx;
            oss << "saturation_" << FluidSystem::phaseName(phaseIdx);
        }
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

        if (Indices::conti0EqIdx <= eqIdx && eqIdx < Indices::conti0EqIdx + numComponents)
            oss << "conti_" << FluidSystem::phaseName(eqIdx - Indices::conti0EqIdx);
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
                const auto &scv = fvElemGeom.subContVol[scvIdx];
                const auto &volVars = elemCtx.volVars(scvIdx, /*timeIdx=*/0);

                tmp = 0;
                this->localResidual().addPhaseStorage(tmp,
                                                      elemCtx,
                                                      scvIdx,
                                                      /*timeIdx=*/0,
                                                      phaseIdx);
                tmp *= scv.volume*volVars.extrusionFactor();
                storage += tmp;
            }
        };

        storage = this->gridView_().comm().sum(storage);
    }

    /*!
     * \copydetails VcfvModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(asImp_(), globalVertexIdx, pvIdx);
        if (tmp > 0)
            // energy related quantity
            return tmp;
        if (Indices::pressure0Idx == pvIdx) {
            // use a pressure gradient of 1e3 Pa/m an intrinisc
            // permeability of 1e-12 as reference (basically, a highly
            // permeable sand stone filled with liquid water.)
            static constexpr Scalar KRef = 1e-12; // [m^2]
            static constexpr Scalar pGradRef = 1e3; // [Pa / m]
            Scalar V = this->boxVolume(globalVertexIdx);

            return std::max(1e-5, pGradRef * KRef / V);
        }
        return 1.0;
    }

    /*!
     * \copydetails VcfvModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(asImp_(), globalVertexIdx, eqIdx);
        if (tmp > 0)
            // energy related equation
            return tmp;

        DUNE_UNUSED int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);

        // make all kg equal
        return 1.0;
    }

protected:
    friend class VcfvModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules available on all model
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkTemperatureModule<TypeTag>(this->problem_()));

        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkEnergyModule<TypeTag>(this->problem_()));
    }

private:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};
}

#include "immisciblepropertydefaults.hh"

#endif
