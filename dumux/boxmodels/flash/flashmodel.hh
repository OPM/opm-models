// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010-2012 by Bernd Flemisch                               *
 *   Copyright (C) 2012 by Klaus Mosthaf                                     *
 *   Copyright (C) 2010-2011 by Melanie Darcis                               *
 *   Copyright (C) 2012 by Philipp Nuske                                     *
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
 * \brief The flash calculation based compositional flow model.
 */
#ifndef DUMUX_FLASH_MODEL_HH
#define DUMUX_FLASH_MODEL_HH

#include "flashproperties.hh"
#include "flashlocalresidual.hh"
#include <dumux/boxmodels/modules/energy/multiphaseenergymodule.hh>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Dumux
{
/*!
 * \ingroup FlashModel
 * \brief A compositional model based on flash-calculations
 */
template<class TypeTag>
class FlashModel: public GET_PROP_TYPE(TypeTag, BaseModel)
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
     * \brief Initialize the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);
    }

    /*!
     * \brief Returns a string with the model's human-readable name
     */
    std::string name() const
    { return "flash"; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
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
     * \brief Given an equation index, return a human readable name.
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
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     */
    void globalPhaseStorage(EqVector &dest, int phaseIdx)
    {
        dest = 0;
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
                dest += tmp;
            }
        };

        dest = this->gridView_().comm().sum(dest);
    }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global vertex index
     * \param pvIdx The primary variable index
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        Scalar tmp = EnergyModule::primaryVarWeight(*this, globalVertexIdx, pvIdx);
        if (tmp > 0)
            return tmp;

        int compIdx = pvIdx - Indices::cTot0Idx;

        // make all kg equal. also, divide the weight of all primary
        // variables by 100 to make the relative errors more
        // comparable to the ones of the other models
        return FluidSystem::molarMass(compIdx) / 1000.0;
    }

    /*!
     * \brief Returns the relative weight of an equation
     *
     * \param globalVertexIdx The global index of the vertex
     * \param eqIdx The index of the primary variable
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

protected:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules meaninful for the model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));

        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Dumux::BoxVtkEnergyModule<TypeTag>(this->problem_()));
    }
};

}

#include "flashpropertydefaults.hh"

#endif
