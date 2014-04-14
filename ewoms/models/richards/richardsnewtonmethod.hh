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
 * \copydoc Ewoms::RichardsNewtonMethod
 */
#ifndef EWOMS_RICHARDS_NEWTON_METHOD_HH
#define EWOMS_RICHARDS_NEWTON_METHOD_HH

#include "richardsproperties.hh"

#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup Newton
 *
 * \brief A Richards model specific Newton method.
 */
template <class TypeTag>
class RichardsNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { pressureWIdx = Indices::pressureWIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { nonWettingPhaseIdx = 1 - wPhaseIdx };

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    RichardsNewtonMethod(Simulator &simulator) : ParentType(simulator)
    {}

    // HACK necessary for GCC 4.4
/*
protected:
    friend class NewtonMethod<TypeTag>;
    friend ParentType;
*/

    /*!
     * \copydoc FvBaseNewtonMethod::update_
     */
    void update_(SolutionVector &uCurrentIter,
                 const SolutionVector &uLastIter,
                 const GlobalEqVector &deltaU,
                 const GlobalEqVector &previousResidual)
    {
        const auto &assembler = this->simulator_.model().jacobianAssembler();
        const auto &simulator = this->simulator_;
        const auto &problem = simulator.problem();

        ParentType::update_(uCurrentIter, uLastIter, deltaU, previousResidual);

        // do not clamp anything after 5 iterations
        if (this->numIterations_ > 4)
            return;

        // clamp saturation change to at most 20% per iteration
        ElementContext elemCtx(simulator);

        ElementIterator elemIt = simulator.gridView().template begin<0>();
        const ElementIterator &elemEndIt = simulator.gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (assembler.elementColor(*elemIt) == JacobianAssembler::Green)
                // don't look at green elements, since they
                // probably have not changed much anyways
                continue;

            elemCtx.updateStencil(*elemIt);

            for (int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
                int globI = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                if (assembler.dofColor(globI) == JacobianAssembler::Green)
                    // don't limit at green DOFs, since they
                    // probably have not changed much anyways
                    continue;

                // calculate the old wetting phase saturation
                const MaterialLawParams &matParams =
                    problem.materialLawParams(elemCtx, dofIdx, /*timeIdx=*/0);

                Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;

                // set the temperatures
                Scalar T = problem.temperature(elemCtx, dofIdx, /*timeIdx=*/0);
                fs.setTemperature(T);

                /////////
                // calculate the phase pressures of the previous iteration
                /////////

                // first, we have to find the minimum capillary pressure
                // (i.e. Sw = 0)
                fs.setSaturation(wPhaseIdx, 1.0);
                fs.setSaturation(nonWettingPhaseIdx, 0.0);
                PhaseVector pC;
                MaterialLaw::capillaryPressures(pC, matParams, fs);

                // non-wetting pressure can be larger than the
                // reference pressure if the medium is fully
                // saturated by the wetting phase
                Scalar pWOld = uLastIter[globI][pressureWIdx];
                Scalar pNOld =
                    std::max(problem.referencePressure(elemCtx, dofIdx, /*timeIdx=*/0),
                             pWOld + (pC[nonWettingPhaseIdx] - pC[wPhaseIdx]));

                /////////
                // find the saturations of the previous iteration
                /////////
                fs.setPressure(wPhaseIdx, pWOld);
                fs.setPressure(nonWettingPhaseIdx, pNOld);

                PhaseVector satOld;
                MaterialLaw::saturations(satOld, matParams, fs);
                satOld[wPhaseIdx] = std::max<Scalar>(0.0, satOld[wPhaseIdx]);

                /////////
                // find the wetting phase pressures which
                // corrospond to a 20% increase and a 20% decrease
                // of the wetting saturation
                /////////
                fs.setSaturation(wPhaseIdx, satOld[wPhaseIdx] - 0.2);
                fs.setSaturation(nonWettingPhaseIdx, 1.0 - (satOld[wPhaseIdx] - 0.2));
                MaterialLaw::capillaryPressures(pC, matParams, fs);
                Scalar pwMin = pNOld - (pC[nonWettingPhaseIdx] - pC[wPhaseIdx]);

                fs.setSaturation(wPhaseIdx, satOld[wPhaseIdx] + 0.2);
                fs.setSaturation(nonWettingPhaseIdx, 1.0 - (satOld[wPhaseIdx] + 0.2));
                MaterialLaw::capillaryPressures(pC, matParams, fs);
                Scalar pwMax = pNOld - (pC[nonWettingPhaseIdx] - pC[wPhaseIdx]);

                /////////
                // clamp the result to the minimum and the maximum
                // pressures we just calculated
                /////////
                Scalar pW = uCurrentIter[globI][pressureWIdx];
                pW = std::max(pwMin, std::min(pW, pwMax));
                uCurrentIter[globI][pressureWIdx] = pW;

                this->model_().invalidateVolumeVariablesCacheEntry(globI, /*timeIdx=*/0);
            }
        }
    }
};
} // namespace Ewoms

#endif
