// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase box model.
 */
#ifndef DUMUX_1P_LOCAL_RESIDUAL_HH
#define DUMUX_1P_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxlocalresidual.hh>

#include "1pvolumevariables.hh"
#include "1pfluxvariables.hh"
#include "1pproperties.hh"

namespace Dumux
{
/*!
 * \ingroup OnePBoxModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase box model.
 */
template<class TypeTag>
class OnePLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;


    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { contiEqIdx = Indices::contiEqIdx };

public:
    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the OneP
     *        model.
     *
     * This function should not include the source and sink terms.
     *  \param storage The phase mass within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);

        // partial time derivative of the wetting phase mass
        storage[contiEqIdx] =
            volVars.fluidState().density(/*phaseIdx=*/0)
            * volVars.porosity();
    }


    /*!
     * \brief Evaluate the mass flux over a face of a sub-control
     *        volume.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param faceIdx The index of the SCV face
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx, 
                     int timeIdx) const
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const auto &up = elemCtx.volVars(fluxVars.upstreamIdx(/*phaseIdx=*/0), timeIdx);
        const auto &dn = elemCtx.volVars(fluxVars.downstreamIdx(/*phaseIdx=*/0), timeIdx);


        flux[contiEqIdx] =
            fluxVars.filterVelocityNormal(/*phaseIdx=*/0)
            * ( fluxVars.upstreamWeight(/*phaseIdx=*/0)
                * up.fluidState().density(/*phaseIdx=*/0)
                +
                fluxVars.downstreamWeight(/*phaseIdx=*/0)
                * dn.fluidState().density(/*phaseIdx=*/0));
    }

    /*!
     * \brief Calculate the source term of the equation.
     *
     * \param q The source/sink in the SCV
     * \param localVertexIdx The index of the SCV
     *
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        elemCtx.problem().source(source, elemCtx, scvIdx, timeIdx);
    }

private:
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
};

}

#endif
