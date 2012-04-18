// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 *
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_1P2C_FLUX_VARIABLES_HH
#define DUMUX_1P2C_FLUX_VARIABLES_HH

#include "1p2cproperties.hh"

#include <dumux/common/math.hh>
#include <dumux/boxmodels/common/boxmultiphasefluxvariables.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the fluxes of the fluid phases over a face of a
 *        finite volume for the one-phase, two-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class OnePTwoCFluxVariables : public BoxMultiPhaseFluxVariables<TypeTag>
{
    typedef BoxMultiPhaseFluxVariables<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    /*!
     * \brief Caclulates the quantities required on a sub-control
     *        volume face for the 2p box model.
     */
    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        ParentType::update(elemCtx, scvfIdx, timeIdx);

        calculateDiffusiveGradients_(elemCtx, scvfIdx, timeIdx);
        calculateDiffCoeffPM_(elemCtx, scvfIdx, timeIdx);
    };

    /*!
    * \brief The binary diffusion coefficient for each fluid phase in the porous medium.
    */
    Scalar porousDiffCoeff() const
    {
        // TODO: tensorial diffusion coefficients
        return diffCoeffPM_;
    };

    /*!
     * \brief Return molar density \f$\mathrm{[mol/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar molarDensity(int phaseIdx = 0) const
    {
        assert(phaseIdx == 0);
        return molarDensity_;
    }

    /*!
     * \brief The molar concentration gradient of a component in a phase.
     *
     * \param compIdx The index of the considered component
     */
    const Vector &moleFracGrad(int phaseIdx = 0, int compIdx = 1) const
    {
        // The 1p2c model is supposed to query only the mole fraction
        // gradient of the second component of the first phase!
        assert(phaseIdx == 0 && compIdx == 1);

        return moleFracGrad_;
    };

private:
    void calculateDiffusiveGradients_(const ElementContext &elemCtx,
                                      int scvfIdx,
                                      int timeIdx)
    {
        // reset all gradients to 0
        moleFracGrad_ = Scalar(0);

        // reset all scalar values to 0
        molarDensity_ = 0;

        typedef typename FVElementGeometry::SubControlVolumeFace Scvf;
        const Scvf &scvf = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx];
        const auto &problem = elemCtx.problem();

        if (!problem.useTwoPointGradient(elemCtx, scvfIdx, timeIdx))
        {
            // calculate gradients
            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++ scvIdx)
            {
                // FE gradient at vertex idx
                const Vector &feGrad = scvf.grad[scvIdx];
                const auto &fs = elemCtx.volVars(scvIdx, timeIdx).fluidState();

                Vector tmp;

                // the mole fraction gradient of the 2nd component [1/m]
                tmp = feGrad;
                tmp *= fs.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);
                moleFracGrad_ += tmp;

                // the molar and the mass density of the fluid
                Scalar feValue = scvf.shapeValue[scvIdx];
                molarDensity_ += feValue * fs.molarDensity(/*phaseIdx=*/0);
            }
        }
        else {
            // use two-point gradients
            const GlobalPosition &posI = elemCtx.element().geometry().corner(this->insideIdx());
            const GlobalPosition &posJ = elemCtx.element().geometry().corner(this->outsideIdx());

            // tmp = pos_i - pos_j
            Vector tmp;
            for (int i=0; i < Vector::dimension; ++i)
                tmp[i] = posI[i] - posJ[i];

            Scalar dist = tmp.two_norm();

            const auto &fsI = elemCtx.volVars(this->insideIdx(), timeIdx).fluidState();
            const auto &fsJ = elemCtx.volVars(this->outsideIdx(), timeIdx).fluidState();

            tmp = scvf.normal;
            tmp /= scvf.normal.two_norm()*dist;

            moleFracGrad_ = tmp;
            moleFracGrad_ *= fsJ.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1) - fsI.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);

            molarDensity_ = (fsJ.molarDensity(/*phaseIdx=*/0) + fsI.molarDensity(/*phaseIdx=*/0))/2;
        }
    }

    /*!
     * \brief Calculation of the effective diffusion coefficient
     *
     * \param problem The considered problem file
     * \param element The considered element of the grid
     * \param elemDat The parameters stored in the considered element
     */
    void calculateDiffCoeffPM_(const ElementContext &elemCtx,
                               int scvfIdx,
                               int timeIdx)
    {
        const auto &volVarsI = elemCtx.volVars(this->insideIdx(), timeIdx);
        const auto &volVarsJ = elemCtx.volVars(this->outsideIdx(), timeIdx);

        // Diffusion coefficient in the porous medium
        diffCoeffPM_
            = harmonicMean(volVarsI.porosity()
                           * volVarsI.tortuosity()
                           * volVarsI.diffCoeff(),
                           volVarsJ.porosity()
                           * volVarsJ.tortuosity()
                           * volVarsJ.diffCoeff());
    }

    //! mole fraction gradient
    Vector moleFracGrad_;

    //! the effective diffusion coefficent in the porous medium
    Scalar diffCoeffPM_;

    //! molar density of the fluid at the integration point
    Scalar molarDensity_;
};

} // end namepace

#endif
