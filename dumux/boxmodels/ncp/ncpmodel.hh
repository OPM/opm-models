// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
#ifndef DUMUX_NCP_MODEL_HH
#define DUMUX_NCP_MODEL_HH

#include "ncpproperties.hh"
#include "energy/ncpvolumevariablesenergy.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dune/common/fvector.hh>

#include <sstream>
#include <string>
#include <vector>
#include <array>

namespace Dumux
{
/*!
 * \ingroup NcpModel
 * \brief A fully implicit model for M-phase, N-component flow using
 *        vertex centered finite volumes.
 *
 * This model implements a \f$M\f$-phase flow of a fluid mixture
 * composed of \f$N\f$ chemical species. The phases are denoted by
 * lower index \f$\alpha \in \{ 1, \dots, M \}\f$. All fluid phases
 * are mixtures of \f$N \geq M - 1\f$ chemical species which are
 * denoted by the upper index \f$\kappa \in \{ 1, \dots, N \} \f$.
 *
 * The standard multi-phase Darcy law is used as the equation for
 * the conservation of momentum:
 * \f[
     v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \boldsymbol{K}
     \left(
       \text{grad}\left(p_\alpha - \varrho_{\alpha} g\right)
     \right)
     \f]
 *
 * By inserting this into the equations for the conservation of the
 * mass of each component, one gets one mass-continuity equation for
 * each component \f$\kappa\f$
 * \f[
 \sum_{\kappa} \left(
    \phi \frac{\partial \varrho_\alpha x_\alpha^\kappa S_\alpha}{\partial t}
    -
    \mathrm{div}\;
    \left\{
       \frac{\varrho_\alpha}{\overline M_\alpha} x_\alpha^\kappa
       \frac{k_{r\alpha}}{\mu_\alpha} \boldsymbol{K}
       \text{grad}\left( p_\alpha - \varrho_{\alpha} g\right)
    \right\}
    \right)
    = q^\kappa
    \f]
 * with \f$\overline M_\alpha\f$ being the average molar mass of the
 * phase \f$\alpha\f$: \f[ \overline M_\alpha = \sum_\kappa M^\kappa
 * \; x_\alpha^\kappa \f]
 *
 * For the missing \f$M\f$ model assumptions, the model assumes that
 * if a fluid phase is not present, the sum of the mole fractions of
 * this fluid phase is smaller than \f$1\f$, i.e.  
 * \f[
 * \forall \alpha: S_\alpha = 0 \implies \sum_\kappa x_\alpha^\kappa \leq 1
 * \f]
 *
 * Also, if a fluid phase may be present at a given spatial location
 * its saturation must be positive:
 * \f[ \forall \alpha: \sum_\kappa x_\alpha^\kappa = 1 \implies S_\alpha \geq 0 \f]
 *
 * Since at any given spatial location, a phase is always either
 * present or not present, one of the strict equalities on the
 * right hand side is always true, i.e.
 * \f[ \forall \alpha: S_\alpha \left( \sum_\kappa x_\alpha^\kappa - 1 \right) = 0 \f]
 * always holds.
 *
 * These three equations constitute a non-linear complementarity
 * problem, which can be solved using so-called non-linear
 * complementarity functions \f$\Phi(a, b)\f$ which have the property
 * \f[\Phi(a,b) = 0 \iff a \geq0 \land b \geq0  \land a \cdot b = 0 \f]
 *
 * Several non-linear complementarity functions have been suggested,
 * e.g. the Fischer-Burmeister function
 * \f[ \Phi(a,b) = a + b - \sqrt{a^2 + b^2} \;. \f]
 * This model uses
 * \f[ \Phi(a,b) = \min \{a,  b \}\;, \f]
 * because of its piecewise linearity.
 *
 * These equations are then discretized using a fully-implicit vertex
 * centered finite volume scheme (often known as 'box'-scheme) for
 * spatial discretization and the implicit Euler method as temporal
 * discretization.
 *
 * The model assumes local thermodynamic equilibrium and uses the
 * following primary variables:
 * - The component fugacities \f$f^1, \dots, f^{N}\f$
 * - The pressure of the first phase \f$p_1\f$
 * - The saturations of the first \f$M-1\f$ phases \f$S_1, \dots, S_{M-1}\f$
 * - Temperature \f$T\f$ if the energy equation is enabled
 */
template<class TypeTag>
class NcpModel : public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    enum { fugacity00Idx = Indices::fugacity00Idx };
    enum { pressure0Idx = Indices::pressure0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { ncp0EqIdx = Indices::ncp0EqIdx };

    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    typedef NcpVolumeVariablesEnergy<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Apply the initial conditions to the model.
     *
     * \param prob The object representing the problem which needs to
     *             be simulated.
     */
    void init(Problem &prob)
    {
        ParentType::init(prob);
        minActivityCoeff_.resize(this->numDofs());
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     */
    void globalPhaseStorage(EqVector &dest, int phaseIdx)
    {
        dest = 0;

        ElementContext elemCtx(this->problem_());
        ElementIterator elemIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            elemCtx.updateFVElemGeom(*elemIt);
            elemCtx.updateScvVars(/*timeIdx=*/0);

            this->localResidual().addPhaseStorage(dest, 
                                                  elemCtx,
                                                  /*timeIdx=*/0,
                                                  phaseIdx);
        };

        dest = this->gridView_().comm().sum(dest);
    }

    /*!
     * \brief Returns a string with the model's human-readable name
     */
    const char *name() const
    { return "ncp"; }

    /*!
     * \brief Given an primary variable index, return a human readable name.
     */
    std::string primaryVarName(int pvIdx) const
    {
        std::string s;
        if (!(s = EnergyModule::primaryVarName(pvIdx)).empty())
            return s;

        std::ostringstream oss;
        if (pvIdx == pressure0Idx)
            oss << "pressure_" << FluidSystem::phaseName(/*phaseIdx=*/0);
        else if (saturation0Idx <= pvIdx && pvIdx < saturation0Idx + (numPhases - 1))
            oss << "saturation_" << FluidSystem::phaseName(/*phaseIdx=*/pvIdx - saturation0Idx);
        else if (fugacity00Idx <= pvIdx && pvIdx < fugacity00Idx + numComponents)
            oss << "fugacity_" << FluidSystem::phaseName(/*phaseIdx=*/0) << "^" << FluidSystem::componentName(pvIdx - fugacity00Idx);
        else
            assert(false);
        
        return oss.str();
    }

    /*!
     * \brief Given an equation index, return a human readable name.
     */
    std::string eqName(int eqIdx) const
    { 
        std::string s;
        if (!(s = EnergyModule::eqName(eqIdx)).empty())
            return s;

        std::ostringstream oss;
        if (conti0EqIdx <= eqIdx && eqIdx < conti0EqIdx + numComponents)
            oss << "continuity^" << FluidSystem::componentName(eqIdx - conti0EqIdx);
        else if (ncp0EqIdx <= eqIdx && eqIdx < ncp0EqIdx + numPhases)
            oss << "ncp_" << FluidSystem::phaseName(/*phaseIdx=*/eqIdx - ncp0EqIdx);
        else
            assert(false);
        
        return oss.str();
    }
    
    /*!
     * \brief Update the weights of all primary variables within an
     *        element given the complete set of volume variables
     *
     * \param element The DUNE codim 0 entity
     * \param volVars All volume variables for the element
     */
    void updatePVWeights(const ElementContext &elemCtx) const
    {
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
            int globalIdx = elemCtx.globalSpaceIndex(scvIdx, /*timeIdx=*/0);

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                minActivityCoeff_[globalIdx][compIdx] = 1e100;
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    const auto &fs = elemCtx.volVars(scvIdx, /*timeIdx=*/0).fluidState();

                    minActivityCoeff_[globalIdx][compIdx] =
                        std::min(minActivityCoeff_[globalIdx][compIdx],
                                 fs.fugacityCoefficient(phaseIdx, compIdx)
                                 * fs.pressure(phaseIdx) );
                };
            };
        };
    }

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param vertIdx The global index of the control volume
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int vertIdx, int pvIdx) const
    {
        if (fugacity00Idx <= pvIdx && pvIdx < fugacity00Idx + numComponents) {
            int compIdx = pvIdx - fugacity00Idx;

            Valgrind::CheckDefined(minActivityCoeff_[vertIdx][compIdx]);
            return 1.0 / minActivityCoeff_[vertIdx][compIdx];
        }

        return ParentType::primaryVarWeight(vertIdx, pvIdx);
    }

    /*!
     * \brief Returns the smallest activity coefficient of a component for the most
     *        current solution at a vertex.
     */
    Scalar minActivityCoeff(int vertIdx, int compIdx) const
    {
        return minActivityCoeff_[vertIdx][compIdx];
    }

protected:
    friend class BoxModel<TypeTag>;

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        this->vtkOutputModules_.push_back(new Dumux::BoxVtkMultiPhaseModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkCompositionModule<TypeTag>(this->problem_()));
        this->vtkOutputModules_.push_back(new Dumux::BoxVtkTemperatureModule<TypeTag>(this->problem_()));
        if (enableEnergy)
            this->vtkOutputModules_.push_back(new Dumux::BoxVtkEnergyModule<TypeTag>(this->problem_()));
    }

    mutable std::vector<ComponentVector> minActivityCoeff_;
};

}

#include "ncppropertydefaults.hh"

#endif
