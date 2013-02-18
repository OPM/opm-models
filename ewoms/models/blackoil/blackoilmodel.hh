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
 * \copydoc Ewoms::BlackOilModel
 */
#ifndef EWOMS_BLACK_OIL_MODEL_HH
#define EWOMS_BLACK_OIL_MODEL_HH

#include "blackoilproperties.hh"
#include "blackoillocalresidual.hh"

#include <sstream>
#include <string>

namespace Ewoms {

/*!
 * \ingroup BlackOilVcfvModel
 * \brief A fully-implicit black-oil flow model using the VCVF discretization.
 *
 * The black-oil model is a three-phase, three-component model widely
 * used for oil reservoir simulation.  The phases are denoted by lower
 * index \f$\alpha \in \{ w, g, o \}\f$ ("water", "gas" and "oil") and
 * the components by upper index \f$\kappa \in \{ W, G, O \}\f$
 * ("Water", "Gas" and "Oil"). The model assumes partial miscibility:
 *
 * - Water and the gas phases are immisicible and are assumed to be
 *   only composed of the water and gas components respectively-
 * - The oil phase is assumed to be a mixture of the gas and the oil
 *  components.
 *
 * The densities of the phases are determined by so-called
 * <i>formation volume factors</i>:
 *
 * \f[
 * B_\alpha := \frac{\varrho_\alpha(1\,\text{bar})}{\varrho_\alpha(p_\alpha)}
 * \f]
 *
 * Since the gas and water phases are assumed to be immiscible, this
 * is sufficint to calculate their density. For the formation volume
 * factor of the the oil phase \f$B_o\f$ determines the density of
 * *saturated* oil, i.e. the density of the oil phase if some gas
 * phase is present.
 *
 * The composition of the oil phase is given by the <i>gas formation
 * factor</i> \f$R_s\f$, which defined as the volume of gas at
 * atmospheric pressure that is dissolved in saturated oil at a given
 * pressure:
 *
 * \f[
 * R_s := \frac{x_o^G(p)\,\varrho_{mol,o}(p)}{\varrho_g(1\,\text{bar})}\;.
 * \f]
 *
 * This allows to calculate all quantities required for the
 * mass-conservation equations for each component, i.e.
 *
 * \f[
 * \sum_\alpha \frac{\partial\;\phi c_\alpha^\kappa S_\alpha }{\partial t}
 * - \sum_\alpha \text{div} \left\{ c_\alpha^\kappa \mathbf{v}_\alpha  \right\}
 * - q^\kappa = 0 \;,
 * \f]
 * where \f$\mathrm{v}_\alpha\f$ is the filter velocity of the phase
 * \f$\alpha\f$.
 *
 * By default \f$\mathrm{v}_\alpha\f$ is determined by using the
 * standard multi-phase Darcy approach, i.e.
 * \f[ \mathbf{v}_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right) \;, \f]
 * although the actual approach which is used can be specified via the
 * \c VelocityModule property. For example, the velocity model can by
 * changed to the Forchheimer approach by
 * \code
 * SET_TYPE_PROP(MyProblemTypeTag, VelocityModule, Ewoms::VcfvForchheimerVelocityModule<TypeTag>);
 * \endcode
 *
 * The primary variables used by this model are:
 * - The pressure of the phase with the lowest index
 * - The two saturations of the phases with the lowest indices
 */
template<class TypeTag >
class BlackOilModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef VcfvModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = FluidSystem::numComponents };

public:
    /*!
     * \brief Register all run-time parameters for the immiscible VCVF discretization.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VcfvVtkBlackOilModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkMultiPhaseModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkCompositionModule<TypeTag>::registerParameters();
        Ewoms::VcfvVtkTemperatureModule<TypeTag>::registerParameters();
    }

    /*!
     * \copydoc VcfvModel::init
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);
    }

    /*!
     * \copydoc VcfvModel::name
     */
    const char *name() const
    { return "blackoil"; }

    /*!
     * \copydoc VcfvModel::primaryVarName
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
     * \copydoc VcfvModel::eqName
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
     * \copydoc VcfvModel::primaryVarWeight
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pressure0Idx == pvIdx) {
            // use a pressure gradient of 1e3 Pa/m an intrinisc
            // permeability of 1e-12 as reference (basically, a highly
            // permeable sand stone filled with liquid water.)
            static constexpr Scalar KRef = 1e-12; // [m^2]
            static constexpr Scalar pGradRef = 1e3; // [Pa / m]
            Scalar r = std::pow(this->boxVolume(globalVertexIdx), 1.0/dimWorld);

            return std::max(1/referencePressure_, pGradRef * KRef / r);
        }
        return 1;
    }

    /*!
     * \copydoc VcfvModel::eqWeight
     */
    Scalar eqWeight(int globalVertexIdx, int eqIdx) const
    {
        int compIdx = eqIdx - Indices::conti0EqIdx;
        assert(0 <= compIdx && compIdx <= numPhases);

        // make all kg equal
        return FluidSystem::molarMass(compIdx);
    }

    /*!
     * \copydoc VcfvModel::updateBegin
     */
    void updateBegin()
    {
        ParentType::updateBegin();

        referencePressure_ = this->solution(/*timeIdx=*/0)[/*vertexIdx=*/0][/*pvIdx=*/Indices::pressure0Idx];
    }

protected:
    friend class VcfvModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules available on all model
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkBlackOilModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Ewoms::VcfvVtkTemperatureModule<TypeTag>(this->problem_()));
    }

    mutable Scalar referencePressure_;
};
}

#include "blackoilpropertydefaults.hh"

#endif
