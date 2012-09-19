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
 * \copydoc Dumux::ImmiscibleModel
 */
#ifndef DUMUX_IMMISCIBLE_MODEL_HH
#define DUMUX_IMMISCIBLE_MODEL_HH

#include "immiscibleproperties.hh"
#include "immisciblelocalresidual.hh"

#include <sstream>
#include <string>

namespace Dumux {
/*!
 * \ingroup ImmiscibleBoxModel
 * \brief A isothermal multi-phase flow model using the box scheme.
 *
 * This model implements multi-phase flow of immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * These equations are discretized by a fully-coupled vertex centered finite volume
 * (box) scheme as spatial and the implicit Euler method as time
 * discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to N. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either

 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */
template<class TypeTag >
class ImmiscibleModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef BoxModel<TypeTag> ParentType;

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
    
public:
    /*!
     * \copydoc BoxModel::name
     */
    const char *name() const
    { return "immiscible"; }

    /*!
     * \copydoc BoxModel::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    { 
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
     * \copydoc BoxModel::eqName
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
     * \copydetails BoxModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (pvIdx == Indices::pressure0Idx) {
            Scalar absPv = std::abs(this->solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx]);
            return std::min(1.0/absPv, 1.0);
        }
        return 1.0;
    }

    /*!
     * \copydetails BoxModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        int DUNE_UNUSED compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);

        // make all kg equal
        return 1.0;
    }

protected:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules available on all model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));

        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Dumux::BoxVtkEnergyModule<TypeTag>(this->problem_()));
    }
};
}

#include "immisciblepropertydefaults.hh"

#endif
