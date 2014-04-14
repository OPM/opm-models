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

#include "richardsproperties.hh"
#include "richardsindices.hh"
#include "richardslocalresidual.hh"
#include "richardsfluxvariables.hh"
#include "richardsratevector.hh"
#include "richardsboundaryratevector.hh"
#include "richardsprimaryvariables.hh"
#include "richardsvolumevariables.hh"
#include "richardsnewtonmethod.hh"

#include <ewoms/models/common/multiphasebasemodel.hh>

#include <opm/material/components/NullComponent.hpp>
#include <opm/material/fluidsystems/LiquidPhase.hpp>
#include <opm/material/fluidsystems/GasPhase.hpp>
#include <opm/material/fluidsystems/2pImmiscibleFluidSystem.hpp>

#include <dune/common/unused.hh>

#include <sstream>
#include <string>

namespace Ewoms {
template <class TypeTag>
class RichardsModel;
}

namespace Opm {
namespace Properties {
//! The type tag for problems discretized using the Richards model
NEW_TYPE_TAG(Richards, INHERITS_FROM(MultiPhaseBaseModel));

//! By default, assume that the first phase is the liquid one
SET_INT_PROP(Richards, LiquidPhaseIndex, 0);

//! The local residual operator
SET_TYPE_PROP(Richards,
              LocalResidual,
              Ewoms::RichardsLocalResidual<TypeTag>);

//! The global model used
SET_TYPE_PROP(Richards, Model, Ewoms::RichardsModel<TypeTag>);

//! the RateVector property
SET_TYPE_PROP(Richards, RateVector, Ewoms::RichardsRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(Richards, BoundaryRateVector, Ewoms::RichardsBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(Richards, PrimaryVariables, Ewoms::RichardsPrimaryVariables<TypeTag>);

//! The class for the volume averaged quantities
SET_TYPE_PROP(Richards, VolumeVariables, Ewoms::RichardsVolumeVariables<TypeTag>);

//! The class for the quantities required for the flux calculation
SET_TYPE_PROP(Richards, FluxVariables, Ewoms::RichardsFluxVariables<TypeTag>);

//! The class of the Newton method
SET_TYPE_PROP(Richards, NewtonMethod, Ewoms::RichardsNewtonMethod<TypeTag>);

//! The class with all index definitions for the model
SET_TYPE_PROP(Richards, Indices, Ewoms::RichardsIndices);

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. Please be aware that you
 * should be careful to use the Richards model in conjunction with
 * liquid non-wetting phases. This is only meaningful if the viscosity
 * of the liquid phase is _much_ lower than the viscosity of the
 * wetting phase.
 */
SET_PROP(Richards, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::LiquidPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

/*!
 * \brief The wetting phase used.
 *
 * By default we use the null-phase, i.e. this has to be defined by
 * the problem for the program to work. This doed not need to be
 * specified by the problem for the Richards model to work because the
 * Richards model does not conserve the non-wetting phase.
 */
SET_PROP(Richards, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

public:
    typedef Opm::GasPhase<Scalar, Opm::NullComponent<Scalar> > type;
};

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the immiscible twophase fluid system. The
 * actual fluids used are specified using in the problem definition by
 * the WettingPhase and NonwettingPhase properties. Be aware that
 * using different fluid systems in conjunction with the Richards
 * model only makes very limited sense.
 */
SET_PROP(Richards, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Opm::FluidSystems::TwoPImmiscible<Scalar, WettingPhase,
                                              NonwettingPhase> type;
};

} // namespace Properties
} // namespace Opm

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
class RichardsModel
    : public MultiPhaseBaseModel<TypeTag>
{
    typedef MultiPhaseBaseModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { dimWorld = GridView::dimensionworld };

public:
    RichardsModel(Simulator &simulator)
        : ParentType(simulator)
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
    }

    /*!
     * \copydoc FvBaseDiscretization::name
     */
    static std::string name()
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
            return 10 / referencePressure_;
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

        // find the a reference pressure. The first degree of freedom
        // might correspond to non-interior entities which would lead
        // to an undefined value, so we have to iterate...
        for (size_t dofIdx = 0; dofIdx < this->numDof(); ++ dofIdx) {
            if (this->dofTotalVolume(dofIdx) > 0) {
                referencePressure_ =
                    this->solution(/*timeIdx=*/0)[dofIdx][/*pvIdx=*/Indices::pressureWIdx];
                break;
            }
        }
    }

    /*!
     * \copydoc FvBaseDiscretization::phaseIsConsidered
     */
    bool phaseIsConsidered(int phaseIdx) const
    { return phaseIdx == wPhaseIdx; }

    void registerOutputModules_()
    {
        ParentType::registerOutputModules_();
    }

    mutable Scalar referencePressure_;
};
} // namespace Ewoms

#endif
