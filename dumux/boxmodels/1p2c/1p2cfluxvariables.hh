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
#include <dumux/common/valgrind.hh>

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
class OnePTwoCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;

    typedef typename GET_PROP_TYPE(TypeTag, OnePTwoCIndices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
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
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief Caclulates the quantities required on a sub-control
     *        volume face for the 2p box model.
     */
    void update(const ElementContext &elemCtx, int scvfIdx)
    {
        insideScvIdx_ = elemCtx.fvElemGeom().subContVolFace[scvfIdx].i;
        outsideScvIdx_ = elemCtx.fvElemGeom().subContVolFace[scvfIdx].j;

        extrusionFactor_ =
            (elemCtx.volVars(insideScvIdx_).extrusionFactor() 
             + elemCtx.volVars(outsideScvIdx_).extrusionFactor()) / 2;

        calculateGradients_(elemCtx, scvfIdx);
        calculateVelocities_(elemCtx, scvfIdx);
        calculateDiffCoeffPM_(elemCtx, scvfIdx);
    };

    /*!
     * \brief Return the extrusion factor of the SCVF.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

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
     * \brief Return density \f$\mathrm{[kg/m^3]}\f$ of a phase at the integration
     *        point.
     */
    Scalar density(int phaseIdx = 0) const
    {
        assert(phaseIdx == 0);
        return density(phaseIdx);
    }

    /*!
     * \brief Return the concentration gradient.
     *
     * \param compIdx The index of the considered component
     */
    const Vector &molarityGrad(int phaseIdx = 0, int compIdx = 1) const
    {
        // The 1p2c model is supposed to query only the
        // concentration gradient of the second component of the
        // first phase!
        assert(phaseIdx == 0 && compIdx == 1);

        return molarityGrad_;
    };
    
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
    
    /*!
     * \brief The mass fraction gradient of a component in a phase.
     *
     * \param compIdx The index of the considered component
     */
    const Vector &massFracGrad(int phaseIdx = 0, int compIdx = 1) const
    {
        // The 1p2c model is supposed to query only the mass fraction
        // gradient of the second component of the first phase!
        assert(phaseIdx == 0 && compIdx == 1);
        
        return massFracGrad_;
    };

    /*!
     * \brief Return a phase's pressure potential gradient.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &potentialGrad(int phaseIdx = 0) const
    { 
        assert(phaseIdx == 0);
        return potentialGrad_;
    }
    
    /*!
     * \brief Return the fluid's filter velocity.
     *
     * For low Reynolds numbers that's the Darcy velocity, for higher
     * ones advanced velocity models like the one suggested by
     * Forchheimer kick in.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Vector &filterVelocity(int phaseIdx) const
    {
        assert(phaseIdx == 0);
        return filterVelocity_;
    }

    /*!
     * \brief Return a phase's pressure potential gradient times
     *        intrinsic permeability times times the mobility timesthe
     *        normal of the sub control volume face times the area of
     *        the SCVF.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar filterVelocityNormal(int phaseIdx) const
    { 
        assert(phaseIdx == 0); // this is a single phase model!
        return filterVelocityNormal_;
    }

    /*!
     * \brief Return the local index of the control volume which is on
     *        the "inside" of the sub-control volume face.
     */
    short insideIdx() const
    { return insideScvIdx_; }

    /*!
     * \brief Return the local index of the control volume which is on
     *        the "outside" of the sub-control volume face.
     */
    short outsideIdx() const
    { return outsideScvIdx_; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the downstream
     *                 direction is requested.
     */
    short downstreamIdx(int phaseIdx) const
    { 
        assert(phaseIdx == 0); // this is a single phase model!
        return (filterVelocityNormal_ > 0)?outsideScvIdx_:insideScvIdx_;
    }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase for which the upstream
     *                 direction is requested.
     */
    short upstreamIdx(int phaseIdx) const
    {
        assert(phaseIdx == 0); // this is a single phase model!
        return (filterVelocityNormal_ > 0)?insideScvIdx_:outsideScvIdx_;
    }

    /*!
     * \brief Return the weight of the upstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar upstreamWeight(int phaseIdx) const
    { return 1.0; }

    /*!
     * \brief Return the weight of the downstream control volume
     *        for a given phase as a function of the normal flux.
     *
     * \param phaseIdx The index of the fluid phase
     */
    Scalar downstreamWeight(int phaseIdx) const
    { return 0.0; }

private:
    void calculateGradients_(const ElementContext &elemCtx,
                             int scvfIdx)
    {
        // reset all gradients to 0
        potentialGrad_ = Scalar(0);
        molarityGrad_ = Scalar(0);
        moleFracGrad_ = Scalar(0);
        massFracGrad_ = Scalar(0);
        
        // reset all scalar values to 0
        density_ = 0;
        molarDensity_ = 0;
        
        typedef typename FVElementGeometry::SubControlVolumeFace Scvf;
        const Scvf &scvf = elemCtx.fvElemGeom().subContVolFace[scvfIdx];
        const auto &spatialParams = elemCtx.problem().spatialParameters();

        if (!spatialParams.useTwoPointGradient(elemCtx, scvfIdx))
        {
            // calculate gradients
            for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++ scvIdx)
            {
                // FE gradient at vertex idx
                const Vector &feGrad = scvf.grad[scvIdx];
                const auto &fs = elemCtx.volVars(scvIdx, /*historyIdx=*/0).fluidState();
                
                Vector tmp;

                // the pressure gradient [Pa/m]
                tmp = feGrad;
                tmp *= fs.pressure(/*phaseIdx=*/0);
                potentialGrad_ += tmp;
                
                // the molarity gradient of the 2nd component [mol/m^3/m]
                tmp = feGrad;
                tmp *= fs.molarity(/*phaseIdx=*/0, /*compIdx=*/1);
                molarityGrad_ += tmp;
                
                // the mole fraction gradient of the 2nd component [1/m]
                tmp = feGrad;
                tmp *= fs.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);
                moleFracGrad_ += tmp;
                
                // the mass fraction gradient of the 2nd component [1/m]
                tmp = feGrad;
                tmp *= fs.massFraction(/*phaseIdx=*/0, /*compIdx=*/1);
                massFracGrad_ += tmp;
                
                // the molar and the mass density of the fluid
                Scalar feValue = scvf.shapeValue[scvIdx];
                density_ += feValue * fs.density(/*phaseIdx=*/0);
                molarDensity_ += feValue * fs.molarDensity(/*phaseIdx=*/0);
            }
        }
        else {
            // use two-point gradients
            const GlobalPosition &posI = elemCtx.element().geometry().corner(insideScvIdx_);
            const GlobalPosition &posJ = elemCtx.element().geometry().corner(outsideScvIdx_);

            // tmp = pos_i - pos_j
            Vector tmp;
            for (int i=0; i < Vector::dimension; ++i)
                tmp[i] = posI[i] - posJ[i];

            Scalar dist = tmp.two_norm();

            const auto &fsI = elemCtx.volVars(insideScvIdx_, /*historyIdx=*/0).fluidState();
            const auto &fsJ = elemCtx.volVars(outsideScvIdx_, /*historyIdx=*/0).fluidState();

            tmp = scvf.normal;
            tmp /= scvf.normal.two_norm()*dist;

            potentialGrad_ = tmp;
            potentialGrad_ *= fsJ.pressure(/*phaseIdx=*/0) - fsI.pressure(/*phaseIdx=*/0);

            molarityGrad_ = tmp;
            molarityGrad_ *= fsJ.molarity(/*phaseIdx=*/0, /*compIdx=*/1) - fsI.molarity(/*phaseIdx=*/0, /*compIdx=*/1);

            moleFracGrad_ = tmp;
            moleFracGrad_ *= fsJ.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1) - fsI.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);

            massFracGrad_ = tmp;
            massFracGrad_ *= fsJ.massFraction(/*phaseIdx=*/0, /*compIdx*/1) - fsI.massFraction(/*phaseIdx=*/0, /*compIdx=*/1);
            
            density_ = (fsJ.density(/*phaseIdx=*/0) + fsI.density(/*phaseIdx=*/0))/2;
            molarDensity_ = (fsJ.molarDensity(/*phaseIdx=*/0) + fsI.molarDensity(/*phaseIdx=*/0))/2;
        }
        
        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity))
        {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector g(elemCtx.problem().gravity(elemCtx, insideScvIdx_));
            g += elemCtx.problem().gravity(elemCtx, outsideScvIdx_);
            g /= 2;
            
            // make gravity acceleration a force
            Vector f(g);
            f *= density_;

            // calculate the final potential gradient
            potentialGrad_ -= f;
        }
    }

    void calculateVelocities_(const ElementContext &elemCtx, 
                              int scvfIdx)
    {
        const SpatialParameters &spatialParams = elemCtx.problem().spatialParameters();

        // calculate the intrinsic permeability
        Tensor K;
        spatialParams.meanK(K,
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                insideScvIdx_),
                            spatialParams.intrinsicPermeability(elemCtx,
                                                                outsideScvIdx_));

        const Vector &normal = elemCtx.fvElemGeom().subContVolFace[scvfIdx].normal;
        
        // calculate the flux in the normal direction of the
        // current sub control volume face:
        //
        // v = - (K grad p) * n
        //
        // (the minus comes from the Darcy law which states that
        // the flux is from high to low pressure potentials.)                           
        K.mv(potentialGrad_, filterVelocity_);
        // velocities go along negative pressure gradients
        filterVelocity_ *= -1;

        // scalar product with the face normal
        filterVelocityNormal_ = 0.0;
        for (int i = 0; i < Vector::size; ++i) 
            filterVelocityNormal_ += filterVelocity_[i]*normal[i];
        
        // multiply both with the upstream mobility
        const auto &up = elemCtx.volVars(upstreamIdx(/*phaseIdx=*/0), /*historyIdx=*/0);
        filterVelocityNormal_ *= up.mobility(/*phaseIdx=*/0);
        filterVelocity_ *= up.mobility(/*phaseIdx=*/0);
    }

    /*!
     * \brief Calculation of the effective diffusion coefficient
     *
     * \param problem The considered problem file
     * \param element The considered element of the grid
     * \param elemDat The parameters stored in the considered element
     */
    void calculateDiffCoeffPM_(const ElementContext &elemCtx, 
                               int scvfIdx)
    {
        const auto &volVarsI = elemCtx.volVars(insideScvIdx_, /*historyIdx=*/0);
        const auto &volVarsJ = elemCtx.volVars(outsideScvIdx_, /*historyIdx=*/0);

        // Diffusion coefficient in the porous medium
        diffCoeffPM_
            = harmonicMean(volVarsI.porosity()
                           * volVarsI.tortuosity()
                           * volVarsI.diffCoeff(),
                           volVarsJ.porosity()
                           * volVarsJ.tortuosity()
                           * volVarsJ.diffCoeff());
    }

    // local indices of the inside and the outside sub-control volumes
    int insideScvIdx_;
    int outsideScvIdx_;

    // extrusion factor
    Scalar extrusionFactor_;

    //! pressure potential gradient
    Vector potentialGrad_;

    //! concentratrion gradient
    Vector molarityGrad_;

    //! mole fraction gradient
    Vector moleFracGrad_;

    //! mass fraction gradient
    Vector massFracGrad_;

    //! the effective diffusion coefficent in the porous medium
    Scalar diffCoeffPM_;

    //! mass density of the fluid at the integration point
    Scalar density_;

    //! molar density of the fluid at the integration point
    Scalar molarDensity_;

    // velocities
    Vector filterVelocity_;

    // normal velocity
    Scalar filterVelocityNormal_;
};

} // end namepace

#endif
