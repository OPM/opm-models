// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
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
 * \brief This file contains the data which is required to calculate
 *        all fluxes (mass of components and energy) over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities at the integration point, etc.
 */
#ifndef DUMUX_2P2CNI_FLUX_VARIABLES_HH
#define DUMUX_2P2CNI_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPTwoCNIModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes (mass of components and energy) over a face of a finite
 *        volume for the non-isothermal two-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */
template <class TypeTag>
class TwoPTwoCNIFluxVariables : public TwoPTwoCFluxVariables<TypeTag>
{
    typedef TwoPTwoCFluxVariables<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> Vector;

public:
    void update(const ElementContext &elemCtx, int scvfIdx)
    {
        ParentType::update(elemCtx, scvfIdx);

        const auto &scvf = elemCtx.fvElemGeom().subContVolFace[scvfIdx];
        // calculate temperature gradient using finite element
        // gradients
        Vector temperatureGrad;
        Vector tmp;
        temperatureGrad = Scalar(0.0);
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            const auto &feGrad = scvf.grad[scvIdx];
            const auto &volVars = elemCtx.volVars(scvIdx, /*historyIdx=*/0);

            tmp = feGrad;
            tmp *= volVars.fluidState().temperature(/*phaseIdx=*/0);
            temperatureGrad += tmp;
        }

        // scalar product of temperature gradient and scvf normal
        temperatureGradNormal_ = 0.0;
        for (int i = 0; i < dimWorld; ++ i)
            temperatureGradNormal_ += scvf.normal[i]*temperatureGrad[i];

        const auto &volVarsInside = elemCtx.volVars(this->insideIdx());
        const auto &volVarsOutside = elemCtx.volVars(this->outsideIdx());

        // arithmetic mean
        heatConductivity_ =
            0.5 * (volVarsInside.heatConductivity()
                   +
                   volVarsOutside.heatConductivity());
        Valgrind::CheckDefined(heatConductivity_);
    }

    /*!
     * \brief The temperature gradient times the face normal [K m^2 / m]
     */
    Scalar temperatureGradNormal() const
    { return temperatureGradNormal_; }

    /*!
     * \brief The total heat conductivity at the face \f$\mathrm{[W/m^2 / (K/m)]}\f$
     */
    Scalar heatConductivity() const
    { return heatConductivity_; }

private:
    Scalar temperatureGradNormal_;
    Scalar heatConductivity_;
};

} // end namepace

#endif
