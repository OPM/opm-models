// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Dumux::BlackOilModel
 */
#ifndef DUMUX_BLACK_OIL_MODEL_HH
#define DUMUX_BLACK_OIL_MODEL_HH

#include "blackoilproperties.hh"
#include "blackoillocalresidual.hh"

#include <sstream>
#include <string>

namespace Dumux {

/*!
 * \ingroup BlackOilBoxModel
 * \brief A fully-implicit black-oil flow model using the box scheme.
 */
template<class TypeTag >
class BlackOilModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = FluidSystem::numComponents };
    
public:
    /*!
     * \copydoc BoxModel::name
     */
    const char *name() const
    { return "blackoil"; }

    /*!
     * \copydoc BoxModel::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::ostringstream oss;
        
        if (pvIdx == Indices::pressure0Idx) {
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        }
        else if (Indices::saturation0Idx <= pvIdx && pvIdx <= Indices::saturation0Idx + numPhases - 1) {
            int phaseIdx = pvIdx - Indices::saturation0Idx;
            oss << "saturation_" << FluidSystem::phaseName(phaseIdx);;
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
     * \copydoc BoxModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (pvIdx == Indices::pressure0Idx) {
            Scalar absPv = std::abs(this->solution(/*timeIdx=*/1)[globalVertexIdx][pvIdx]);
            return std::min(1.0/absPv, 1.0);
        }
        return 1;
    }

    /*!
     * \copydoc BoxModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

protected:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules available on all model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
    }
};
}

#include "blackoilpropertydefaults.hh"

#endif
