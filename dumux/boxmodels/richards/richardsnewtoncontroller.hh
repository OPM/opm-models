// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
/*!
 * \file
 * \brief A newton solver specific to the Richards problem.
 */
#ifndef DUMUX_RICHARDS_NEWTON_CONTROLLER_HH
#define DUMUX_RICHARDS_NEWTON_CONTROLLER_HH

#include "richardsproperties.hh"

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup Newton
 * \brief A Richards model specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * and can thus do update smarter than the plain Newton controller.
 */
template <class TypeTag>
class RichardsNewtonController : public NewtonController<TypeTag>
{
    typedef NewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, RichardsIndices) Indices;
    enum {
        dim = GridView::dimension,
        pwIdx = Indices::pwIdx,
    };

public:
    /*!
     * \brief Constructor
     */
    RichardsNewtonController(const Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Update the current solution of the newton method
     *
     * This is basically the step
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param uCurrentIter The solution after the current Newton iteration \f$ u^{k+1} \f$
     * \param uLastIter The solution after the last Newton iteration \f$ u^k \f$
     * \param deltaU The vector of differences between the last
     *               iterative solution and the next one \f$ \Delta u^k \f$
     */
    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const SolutionVector &deltaU)
    {
        ParentType::newtonUpdate(uCurrentIter, uLastIter, deltaU);

        const auto &assembler = this->model_().jacobianAssembler();
        if (!GET_PARAM(TypeTag, bool, NewtonUseLineSearch))
        {
            // do not clamp anything after 5 iterations
            if (this->numSteps_ > 4)
                return;

            // clamp saturation change to at most 20% per iteration
            ElementContext elemCtx(this->problem_());
            
            ElementIterator elemIt = this->gridView_().template begin<0>();
            const ElementIterator &elemEndIt = this->gridView_().template end<0>();
            for (; elemIt != elemEndIt; ++elemIt)
            {
                if (assembler.elementColor(*elemIt) == JacobianAssembler::Green)
                    // don't look at green elements, since they
                    // probably have not changed much anyways
                    continue;

                elemCtx.updateFVElemGeom(*elemIt);
                
                for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
                    int globI = this->problem_().vertexMapper().map(*elemIt, scvIdx, dim);
                    if (assembler.vertexColor(globI) == JacobianAssembler::Green)
                        // don't limit at green vertices, since they
                        // probably have not changed much anyways
                        continue;

                    // calculate the old wetting phase saturation
                    const SpatialParameters &spatialParams = this->problem_().spatialParameters();
                    const MaterialLawParams &matParams = spatialParams.materialLawParams(elemCtx, scvIdx);

                    Scalar pcMin = MaterialLaw::pC(matParams, /*Sw=*/1.0);
                    Scalar pW = uLastIter[globI][pwIdx];
                    Scalar pN = std::max(this->problem_().referencePressure(elemCtx, scvIdx),
                                         pW + pcMin);
                    Scalar pcOld = pN - pW;
                    Scalar SwOld = std::max<Scalar>(0.0, MaterialLaw::Sw(matParams, pcOld));

                    // convert into minimum and maximum wetting phase
                    // pressures
                    Scalar pwMin = pN - MaterialLaw::pC(matParams, SwOld - 0.2);
                    Scalar pwMax = pN - MaterialLaw::pC(matParams, SwOld + 0.2);

                    // clamp the result
                    pW = uCurrentIter[globI][pwIdx];
                    pW = std::max(pwMin, std::min(pW, pwMax));
                    uCurrentIter[globI][pwIdx] = pW;
                }
            }
        }
    }
};
}

#endif
