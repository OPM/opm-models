// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Element-wise calculation of the residual for the Richards box model.
 */
#ifndef DUMUX_RICHARDS_LOCAL_RESIDUAL_HH
#define DUMUX_RICHARDS_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxlocalresidual.hh>

#include "richardsvolumevariables.hh"

#include "richardsfluxvariables.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the residual for the Richards box model.
 */
template<class TypeTag>
class RichardsLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVariables) ElementVariables;

    typedef typename GET_PROP_TYPE(TypeTag, RichardsIndices) Indices;
    enum {
        contiEqIdx = Indices::contiEqIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dimWorld = GridView::dimensionworld};

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    /*!
     * \brief Constructor. Sets the upwind weight.
     */
    RichardsLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
    };

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the Richards
     *        model.
     *
     * This function should not include the source and sink terms.
     *
     * \param result Stores the average mass per unit volume for each phase \f$\mathrm{[kg/m^3]}\f$
     * \param scvIdx The sub control volume index of the current element
     * \param usePrevSol Calculate the storage term of the previous solution
     *                   instead of the model's current solution.
     */
    void computeStorage(PrimaryVariables &result, 
                        const ElementVariables &elemVars,
                        int scvIdx, 
                        int historyIdx) const
    {
        const VolumeVariables &volVars = elemVars.volVars(scvIdx, historyIdx);

        // partial time derivative of the wetting phase mass
        result[contiEqIdx] =
            volVars.density(wPhaseIdx)
            * volVars.saturation(wPhaseIdx)
            * volVars.porosity();
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     *
     *
     * \param flux Stores the total mass fluxes over a sub-control volume face
     *             of the current element \f$\mathrm{[kg/s]}\f$
     * \param scvfIdx The sub control volume face index inside the current
     *                element
     */
    void computeFlux(PrimaryVariables &flux,
                     const ElementVariables &elemVars,
                     int scvfIdx) const
    {
        const auto &fluxVarsEval = elemVars.evalPointFluxVars(scvfIdx);
        const auto &fluxVars = elemVars.fluxVars(scvfIdx);

        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up = elemVars.volVars(fluxVarsEval.upstreamIdx());
        const VolumeVariables &dn = elemVars.volVars(fluxVarsEval.upstreamIdx());

        flux[contiEqIdx] =
            fluxVars.normalFlux()
            *
            ((    massUpwindWeight_)*up.density(wPhaseIdx)*up.mobility()
             +
             (1 - massUpwindWeight_)*dn.density(wPhaseIdx)*dn.mobility());
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param q Stores the average source term of all phases inside a
     *          sub-control volume of the current element \f$\mathrm{[kg/(m^3 * s)]}\f$
     * \param scvIdx The sub control volume index inside the current
     *               element
     */
    void computeSource(PrimaryVariables &q,
                       const ElementVariables &elemVars,
                       int scvIdx) const
    {
        elemVars.problem().source(q,
                                  elemVars,
                                  scvIdx);
    }

private:
    Scalar massUpwindWeight_;
};

}

#endif
