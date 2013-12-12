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
 * \copydoc Ewoms::RichardsModel
 */
#ifndef EWOMS_RICHARDS_MODEL_HH
#define EWOMS_RICHARDS_MODEL_HH

#include "richardslocalresidual.hh"

#include <dune/common/unused.hh>

#include <sstream>
#include <string>

namespace Ewoms {

/*!
 * \ingroup RichardsModel
 *
 * \brief This model implements a variant of the Richards equation for
 *        quasi-twophase flow.
 *
 * In the unsaturated zone, Richards' equation is frequently used to
 * approximate the water distribution above the groundwater level. It
 * can be derived from the two-phase equations, i.e.
 * \f[
 * \frac{\partial\;\phi S_\alpha \rho_\alpha}{\partial t}
 * -
 * \mathrm{div} \left\{
 * \rho_\alpha \frac{k_{r\alpha}}{\mu_\alpha}\; \mathbf{K}\;
 * \mathbf{grad}\left[
 * p_\alpha - g\rho_\alpha
 * \right]
 * \right\}
 * =
 * q_\alpha,
 * \f]
 * where \f$\alpha \in \{w, n\}\f$ is the index of the fluid phase,
 * \f$\rho_\alpha\f$ is the fluid density, \f$S_\alpha\f$ is the fluid
 * saturation, \f$\phi\f$ is the porosity of the soil,
 * \f$k_{r\alpha}\f$ is the relative permeability for the fluid,
 * \f$\mu_\alpha\f$ is the fluid's dynamic viscosity, \f$\mathbf{K}\f$
 * is the intrinsic permeability tensor, \f$p_\alpha\f$ is the fluid
 * phase pressure and \f$g\f$ is the potential of the gravity field.
 *
 * In contrast to the "full" two-phase model, the Richards model
 * assumes that the non-wetting fluid is gas and that it thus exhibits
 * a much lower viscosity than the (liquid) wetting phase. (This
 * assumption is quite realistic in many applications: For example, at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term
 * typically is much larger for the gas phase than for the wetting
 * phase. Using this reasoning, the Richards model assumes that
 * \f$\frac{k_{rn}}{\mu_n}\f$ is infinitely large compared to the same
 * term of the liquid phase. This implies that the pressure of the gas
 * phase is equivalent to the static pressure distribution and that
 * therefore, mass conservation only needs to be considered for the
 * liquid phase.
 *
 * The model thus choses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 * \f[
 * S_w = p_c^{-1}(p_n - p_w)
 * \f]
 * holds, where \f$p_n\f$ is a reference pressure given by the
 * problem's \c referencePressure() method. Nota bene, that the last
 * step assumes that the capillary pressure-saturation curve can be
 * uniquely inverted, i.e. it is not possible to set the capillary
 * pressure to zero if the Richards model ought to be used!
 */
template <class TypeTag>
class RichardsModel : public GET_PROP_TYPE(TypeTag, Discretization)
{
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { dimWorld = GridView::dimensionworld };

public:
    RichardsModel(Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \copydoc FvBaseDiscretization::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
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
    std::string name() const
    { return "richards"; }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarName
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::ostringstream oss;
        if (pvIdx == Indices::pressureWIdx)
            oss << "pressure_" << FluidSystem::phaseName(wPhaseIdx);
        else
            assert(0);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::eqName
     */
    std::string eqName(int eqIdx) const
    {
        std::ostringstream oss;
        if (eqIdx == Indices::contiWEqIdx)
            oss << "continuity_" << FluidSystem::phaseName(wPhaseIdx);
        else
            assert(0);

        return oss.str();
    }

    /*!
     * \copydoc FvBaseDiscretization::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalDofIdx, int pvIdx) const
    {
        if (Indices::pressureWIdx == pvIdx) {
            // use a pressure gradient of 1e2 Pa/m for liquid water as
            // a reference
            Scalar KRef = intrinsicPermeability_[globalDofIdx];
            static const Scalar muRef = 1e-3;
            static const Scalar pGradRef = 1e-2; // [Pa / m]
            Scalar r = std::pow(this->boxVolume(globalDofIdx), 1.0/dimWorld);

            return std::max(10 / referencePressure_, pGradRef * KRef / muRef / r);
        }

        return 1;
    }

    /*!
     * \copydoc FvBaseDiscretization::eqWeight
     */
    Scalar eqWeight(int globalDofIdx, int eqIdx) const
    {
        int DUNE_UNUSED compIdx = eqIdx - Indices::contiWEqIdx;
        assert(0 <= compIdx && compIdx <= FluidSystem::numPhases);

        // make all kg equal
        return 1.0;
    }

    /*!
     * \copydoc FvBaseDiscretization::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        referencePressure_ = this->solution(/*timeIdx=*/0)[/*dofIdx=*/0][/*pvIdx=*/Indices::pressureWIdx];
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
     * \copydoc FvBaseDiscretization::phaseIsConsidered
     */
    bool phaseIsConsidered(int phaseIdx) const
    { return phaseIdx == wPhaseIdx; }

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();
    }

    mutable Scalar referencePressure_;
    mutable std::vector<Scalar> intrinsicPermeability_;
};
} // namespace Ewoms

#include "richardspropertydefaults.hh"

#endif
