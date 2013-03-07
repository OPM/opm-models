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
 * \copydoc Ewoms::FlashModel
 */
#ifndef EWOMS_FLASH_MODEL_HH
#define EWOMS_FLASH_MODEL_HH

#include "flashproperties.hh"

#include <ewoms/models/modules/diffusion/vcfvdiffusionmodule.hh>
#include <ewoms/models/modules/energy/vcfvenergymodule.hh>
#include <ewoms/disc/vcfv/vcfvmodel.hh>

#include <sstream>
#include <string>

namespace Ewoms {

/*!
 * \ingroup FlashModel
 *
 * \brief A compositional multi-phase model based on flash-calculations
 *
 * \brief A generic compositional multi-phase model using primary-variable switching.
 *
 * This model assumes a flow of \f$M \geq 1\f$ fluid phases
 * \f$\alpha\f$, each of which is assumed to be a mixture \f$N \geq
 * M\f$ chemical species (denoted by the upper index \f$\kappa\f$).
 *
 * By default, the standard multi-phase Darcy approach is used to determine
 * the velocity, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
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
 * - \sum_\alpha \text{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha  \right\}
 * - q^\kappa = 0 \;.
 * \f]
 *
 * To determine the quanties that occur in the equations above, this
 * model uses <i>flash calculations</i>. A flash solver starts with
 * the total mass or molar mass per volume for each component and,
 * calculates the compositions, saturations and pressures of all
 * phases at a given temperature. For this the flash solver has to use
 * some model assumptions internally. (Often these are the same
 * primary variable switching or NCP assumptions as used by the other
 * fully implicit compositional multi-phase models provided by eWoms.)
 *
 * Using flash calculations for the flow model has some disadvantages:
 * - The accuracy of the flash solver needs to be sufficient to
 *   calculate the parital derivatives using numerical differentiation
 *   which are required for the Newton scheme.
 * - Flash calculations tend to be quite computationally expensive and
 *   are often numerically unstable.
 *
 * It is thus adviced to increase the target tolerance of the Newton
 * scheme or a to use type for scalar values which exhibits higher
 * precision than the standard \c double (e.g. \c quad) if this model
 * ought to be used.
 *
 * The model uses the following primary variables:
 * - The total molar concentration of each component:
 *   \f$c^\kappa = \sum_\alpha S_\alpha x_\alpha^\kappa \rho_{mol, \alpha}\f$
 * - The absolute temperature $T$ in Kelvins if the energy equation enabled.
 */
template<class TypeTag>
class FlashModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef VcfvModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

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
        Ewoms::VcfvVtkMultiPhaseModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkCompositionModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkTemperatureModule<TypeTag>::registerParameters();

        if (enableDiffusion)
            Ewoms::VcfvVtkDiffusionModule<TypeTag>::registerParameters();

        if (enableEnergy)
            Ewoms::VcfvVtkEnergyModule<TypeTag>::registerParameters();

        REGISTER_PARAM(TypeTag, Scalar, FlashTolerance, "The maximum tolerance for the flash solver to consider the solution converged");
    }

    /*!
     * \copydoc VcfvModel::name
     */
    std::string name() const
    { return "flash"; }

    /*!
     * \copydoc VcfvModel::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        const std::string &tmp = EnergyModule::primaryVarName(pvIdx);
        if (tmp != "")
            return tmp;

        std::ostringstream oss;
        if (Indices::cTot0Idx <= pvIdx && pvIdx < Indices::cTot0Idx + numComponents)
            oss << "c_tot," << FluidSystem::componentName(/*compIdx=*/pvIdx - Indices::cTot0Idx);
        else
            assert(false);

        return oss.str();
    }

    /*!
     * \copydoc VcfvModel::eqName
     */
    std::string eqName(int eqIdx) const
    {
        const std::string &tmp = EnergyModule::eqName(eqIdx);
        if (tmp != "")
            return tmp;

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
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);

            const auto &fvElemGeom = elemCtx.fvElemGeom(/*timeIdx=*/0);

            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
                if (elemIt->partitionType() != Dune::InteriorEntity)
                    continue; // ignore ghost and overlap elements

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
     * \copydoc VcfvModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalVertexIdx, pvIdx);
        if (tmp > 0)
            return tmp;

        int compIdx = pvIdx - Indices::cTot0Idx;

        // make all kg equal. also, divide the weight of all total
        // compositions by 100 to make the relative errors more
        // comparable to the ones of the other models (at 10% porosity
        // the medium is fully saturated with water at atmospheric
        // conditions if 100 kg/m^3 are present!)
        return FluidSystem::molarMass(compIdx) / 100.0;
    }

    /*!
     * \copydoc VcfvModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        Scalar tmp = EnergyModule::eqWeight(*this, globalVertexIdx, eqIdx);
        if (tmp > 0)
            return tmp;

        int compIdx = eqIdx - Indices::conti0EqIdx;

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

private:
    friend class VcfvModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules meaninful for the model
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkTemperatureModule<TypeTag>(this->problem_()));
        if (enableDiffusion)
            this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkDiffusionModule<TypeTag>(this->problem_()));
        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkEnergyModule<TypeTag>(this->problem_()));
    }
};

}

#include "flashpropertydefaults.hh"

#endif
