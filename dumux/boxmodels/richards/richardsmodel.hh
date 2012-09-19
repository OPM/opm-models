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
 * \copydoc Dumux::RichardsModel
 */
#ifndef DUMUX_RICHARDS_MODEL_HH
#define DUMUX_RICHARDS_MODEL_HH

#include "richardslocalresidual.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

#include <sstream>
#include <string>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 *
 * \brief This model which implements a variant of the Richards
 *        equation for quasi-twophase flow.
 *
 * In the unsaturated zone, Richards' equation is frequently used to
 * approximate the water distribution above the groundwater level. It
 * can be derived from the two-phase equations, i.e.
 \f[
 \frac{\partial\;\phi S_\alpha \rho_\alpha}{\partial t}
 -
 \text{div} \left\{
 \rho_\alpha \frac{k_{r\alpha}}{\mu_\alpha}\; \mathbf{K}
 \textbf{grad}\left[
 p_\alpha - g\rho_\alpha
 \right]
 \right\}
 =
 q_\alpha,
 \f]
 * where \f$\alpha \in \{w, n\}\f$ is the fluid phase,
 * \f$\rho_\alpha\f$ is the fluid density, \f$S_\alpha\f$ is the fluid
 * saturation, \f$\phi\f$ is the porosity of the soil,
 * \f$k_{r\alpha}\f$ is the relative permeability for the fluid,
 * \f$\mu_\alpha\f$ is the fluid's dynamic viscosity, \f$\mathbf{K}\f$ is the
 * intrinsic permeability, \f$p_\alpha\f$ is the fluid pressure and
 * \f$g\f$ is the potential of the gravity field.
 *
 * In contrast to the full two-phase model, the Richards model assumes
 * gas as the non-wetting fluid and that it exhibits a much lower
 * viscosity than the (liquid) wetting phase. (For example at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term
 * typically is much larger for the gas phase than for the wetting
 * phase. For this reason, the Richards model assumes that
 * \f$\frac{k_{rn}}{\mu_n}\f$ is infinitly large. This implies that
 * the pressure of the gas phase is equivalent to the static pressure
 * distribution and that therefore, mass conservation only needs to be
 * considered for the wetting phase.
 *
 * The model thus choses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 \f[
 S_w = p_c^{-1}(p_n - p_w)
 \f]
 * holds, where \f$p_n\f$ is a given reference pressure. Nota bene,
 * that the last step assumes that the capillary
 * pressure-saturation curve can be uniquely inverted, so it is not
 * possible to set the capillary pressure to zero when using the
 * Richards model!
 */
template<class TypeTag >
class RichardsModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };

public:
    /*!
     * \copydoc BoxModel::name
     */
    std::string name() const
    { return "richards"; }

    /*!
     * \copydoc BoxModel::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    { 
        std::ostringstream oss;
        if (pvIdx == Indices::pwIdx)
            oss << "pressure_" << FluidSystem::phaseName(wPhaseIdx);
        else
            assert(0);
       
        return oss.str();
    }

    /*!
     * \copydoc BoxModel::eqName
     */
    std::string eqName(int eqIdx) const
    { 
        std::ostringstream oss;
        if (eqIdx == Indices::contiEqIdx)
            oss << "continuity_" << FluidSystem::phaseName(wPhaseIdx);
        else
            assert(0);

        return oss.str();
    }

    /*!
     * \copydoc BoxModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pwIdx == pvIdx)
            return 1e-6;
        return 1;
    }

    /*!
     * \copydoc BoxModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        int DUNE_UNUSED compIdx = eqIdx - Indices::contiEqIdx;
        assert(0 <= compIdx && compIdx <= FluidSystem::numPhases);

        // make all kg equal
        return 1.0;
    }


private:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules available on all model
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
    }
};
}

#include "richardspropertydefaults.hh"

#endif
