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
 * \copydoc Ewoms::RichardsNewtonMethod
 */
#ifndef EWOMS_RICHARDS_NEWTON_METHOD_HH
#define EWOMS_RICHARDS_NEWTON_METHOD_HH

#include "richardsproperties.hh"

#include <ewoms/disc/vcfv/vcfvnewtonmethod.hh>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>

#include <dune/common/fvector.hh>

namespace Ewoms {

/*!
 * \ingroup Newton
 *
 * \brief A Richards model specific Newton method.
 */
template <class TypeTag>
class RichardsNewtonMethod : public VcfvNewtonMethod<TypeTag>
{
    typedef VcfvNewtonMethod<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { dim = GridView::dimension };
    enum { pressureWIdx = Indices::pressureWIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { wPhaseIdx = GET_PROP_VALUE(TypeTag, LiquidPhaseIndex) };
    enum { nPhaseIdx = 1 - wPhaseIdx };

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

public:
    RichardsNewtonMethod(Problem &problem) : ParentType(problem)
    {}

protected:
    friend class NewtonMethod<TypeTag>;
    friend class VcfvNewtonMethod<TypeTag>;

    /*!
     * \copydoc VcfvNewtonMethod::update_
     */
    void update_(SolutionVector &uCurrentIter, const SolutionVector &uLastIter,
                 const GlobalEqVector &deltaU)
    {
        const auto &assembler = this->problem().model().jacobianAssembler();
        const auto &problem = this->problem();

        ParentType::update_(uCurrentIter, uLastIter, deltaU);

        if (!this->enableLineSearch_()) {
            // do not clamp anything after 5 iterations
            if (this->numIterations_ > 4)
                return;

            // clamp saturation change to at most 20% per iteration
            ElementContext elemCtx(problem);

            ElementIterator elemIt = problem.gridView().template begin<0>();
            const ElementIterator &elemEndIt
                = problem.gridView().template end<0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                if (assembler.elementColor(*elemIt) == JacobianAssembler::Green)
                    // don't look at green elements, since they
                    // probably have not changed much anyways
                    continue;

                elemCtx.updateFVElemGeom(*elemIt);

                for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
                    int globI = problem.vertexMapper().map(*elemIt, scvIdx, dim);
                    if (assembler.vertexColor(globI) == JacobianAssembler::Green)
                        // don't limit at green vertices, since they
                        // probably have not changed much anyways
                        continue;

                    // calculate the old wetting phase saturation
                    const MaterialLawParams &matParams
                        = problem.materialLawParams(elemCtx, scvIdx,
                                                    /*timeIdx=*/0);

                    Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;

                    // set the temperatures
                    Scalar T = elemCtx.problem().temperature(elemCtx, scvIdx,
                                                             /*timeIdx=*/0);
                    fs.setTemperature(T);

                    /////////
                    // calculate the phase pressures of the previous iteration
                    /////////

                    // first, we have to find the minimum capillary pressure
                    // (i.e. Sw = 0)
                    fs.setSaturation(wPhaseIdx, 1.0);
                    fs.setSaturation(nPhaseIdx, 0.0);
                    PhaseVector pC;
                    MaterialLaw::capillaryPressures(pC, matParams, fs);

                    // non-wetting pressure can be larger than the
                    // reference pressure if the medium is fully
                    // saturated by the wetting phase
                    Scalar pWOld = uLastIter[globI][pressureWIdx];
                    Scalar pNOld
                        = std::max(problem.referencePressure(elemCtx, scvIdx,
                                                             /*timeIdx=*/0),
                                   pWOld + (pC[nPhaseIdx] - pC[wPhaseIdx]));

                    /////////
                    // find the saturations of the previous iteration
                    /////////
                    fs.setPressure(wPhaseIdx, pWOld);
                    fs.setPressure(nPhaseIdx, pNOld);

                    PhaseVector satOld;
                    MaterialLaw::saturations(satOld, matParams, fs);
                    satOld[wPhaseIdx] = std::max<Scalar>(0.0, satOld[wPhaseIdx]);

                    /////////
                    // find the wetting phase pressures which
                    // corrospond to a 20% increase and a 20% decrease
                    // of the wetting saturation
                    /////////
                    fs.setSaturation(wPhaseIdx, satOld[wPhaseIdx] - 0.2);
                    fs.setSaturation(nPhaseIdx, 1.0 - (satOld[wPhaseIdx] - 0.2));
                    MaterialLaw::capillaryPressures(pC, matParams, fs);
                    Scalar pwMin = pNOld - (pC[nPhaseIdx] - pC[wPhaseIdx]);

                    fs.setSaturation(wPhaseIdx, satOld[wPhaseIdx] + 0.2);
                    fs.setSaturation(nPhaseIdx, 1.0 - (satOld[wPhaseIdx] + 0.2));
                    MaterialLaw::capillaryPressures(pC, matParams, fs);
                    Scalar pwMax = pNOld - (pC[nPhaseIdx] - pC[wPhaseIdx]);

                    /////////
                    // clamp the result to the minimum and the maximum
                    // pressures we just calculated
                    /////////
                    Scalar pW = uCurrentIter[globI][pressureWIdx];
                    pW = std::max(pwMin, std::min(pW, pwMax));
                    uCurrentIter[globI][pressureWIdx] = pW;
                }
            }
        }
    }
};
} // namespace Ewoms

#endif
