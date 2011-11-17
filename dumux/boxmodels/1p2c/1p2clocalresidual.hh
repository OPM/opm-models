// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Katherina Baber                                   *
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
 * \brief Element-wise calculation the local Jacobian for the single-phase,
 *        two-component model in the BOX scheme.
 */
#ifndef DUMUX_ONEP_TWOC_LOCAL_RESIDUAL_HH
#define DUMUX_ONEP_TWOC_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dumux/boxmodels/1p2c/1p2cproperties.hh>
#include <dumux/boxmodels/1p2c/1p2cvolumevariables.hh>
#include <dumux/boxmodels/1p2c/1p2cfluxvariables.hh>

#warning "TODO: implement outflow condition properly"
//#include <dumux/boxmodels/1p2c/1p2cboundaryvariables.hh>

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dumux
{
/*!
 *
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxLocalResidual
 * \brief Calculate the local Jacobian for the single-phase,
 *        two-component model in the BOX scheme.
 *
 *  This class is used to fill the gaps in BoxLocalResidual for the 1P-2C flow.
 */
template<class TypeTag>
class OnePTwoCLocalResidual : public BoxLocalResidual<TypeTag>
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, OnePTwoCIndices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,
        // indices of the equations
        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx
    };

public:
    /*!
     * \brief Evaluate the amount of all conservation quantities
     *        (e.g. phase mass) within a finite volume.
     *
     *        \param result The mass of the component within the sub-control volume
     *        \param scvIdx The index of the considered face of the sub control volume
     *        \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(EqVector &result,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        // retrieve the volume variables for the SCV at the specified
        // point in time
        const VolumeVariables &volVars = elemCtx.volVars(scvIdx, timeIdx);
        const auto &fs = volVars.fluidState();

        result = 0;

        // storage term of continuity equation- molefractions
        //careful: molarDensity changes with moleFrac!
        result[contiEqIdx] +=
            fs.molarDensity(/*phaseIdx=*/0)
            * volVars.porosity();
        // storage term of the transport equation - molefractions
        result[transEqIdx] +=
            fs.molarDensity(/*phaseIdx=*/0)
            * fs.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1)
            * volVars.porosity();
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     *
     *        \param flux The flux over the SCV (sub-control-volume) face for each component
     *        \param faceId The index of the considered face of the sub control volume
     *        \param onBoundary If the considered face exists at the boundary
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        flux = 0;
        asImp_().computeAdvectiveFlux(flux, elemCtx, scvfIdx, timeIdx);
        asImp_().computeDiffusiveFlux(flux, elemCtx, scvfIdx, timeIdx);
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume or a boundary segment.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeAdvectiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    {
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);
        const FluxVariables &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx, timeIdx);

        ////////
        // advective fluxes of all components in all phases
        ////////

        // data attached to upstream and the downstream vertices
        // of the current phase
        const VolumeVariables &up =
            elemCtx.volVars(evalPointFluxVars.upstreamIdx(/*phaseIdx=*/0), timeIdx);
        const VolumeVariables &dn =
            elemCtx.volVars(evalPointFluxVars.downstreamIdx(/*phaseIdx=*/0), timeIdx);
        const auto &fsUp = up.fluidState();
        const auto &fsDn = dn.fluidState();

        Scalar rhoUp = fsUp.molarDensity(/*phaseIdx=*/0);
        Scalar rhoDn = fsDn.molarDensity(/*phaseIdx=*/0);
        Scalar xUp = fsUp.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);
        Scalar xDn = fsDn.moleFraction(/*phaseIdx=*/0, /*compIdx=*/1);

        // total mass/molar flux
        flux[contiEqIdx] +=
            fluxVars.filterVelocityNormal(/*phaseIdx=*/0) *
            (fluxVars.upstreamWeight(/*phaseIdx=*/0)*rhoUp
             +
             fluxVars.downstreamWeight(/*phaseIdx=*/0)*rhoDn);

        // advective flux of the second component
        flux[transEqIdx] +=
            fluxVars.filterVelocityNormal(/*phaseIdx=*/0) *
            (fluxVars.upstreamWeight(/*phaseIdx=*/0)*rhoUp*xUp
             +
             fluxVars.downstreamWeight(/*phaseIdx=*/0)*rhoDn*xDn);
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * \param flux The diffusive flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV
     */
    void computeDiffusiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx) const
    {
        const auto &normal = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].normal;
        const FluxVariables &fluxVars = elemCtx.fluxVars(scvfIdx, timeIdx);

        // molar formulation
        
        // diffusive flux of the second component
        Scalar tmp = 0;
        for (int i = 0; i < normal.dimension; ++ i)
            tmp += fluxVars.moleFracGrad(/*phaseIdx=*/0, /*compIdx=*/1)[i]*normal[i];
        tmp *= -1;
        tmp *= fluxVars.porousDiffCoeff() * fluxVars.molarDensity();
        
#if 0
        // dispersive flux of second component - molefraction
        Vector normalDisp;
        fluxVars.dispersionTensor().mv(normal, normalDisp);
        tmp -= fluxVars.molarDensityAtIP()*
            (normalDisp * fluxVars.moleFracGrad());
#endif
        
        flux[transEqIdx] += tmp;
    }

    /*!
     * \brief Calculate the source term of the equation
     *        \param q The source/sink in the SCV for each component
     *        \param localVertexIdx The index of the vertex of the sub control volume
     *
     */
    void computeSource(RateVector &q,
                       const ElementContext &elemCtx,
                       int scvIdx, 
                       int timeIdx) const
    {
        elemCtx.problem().source(q,
                                  elemCtx,
                                  scvIdx,
                                  timeIdx);
    }

protected:
#warning "TODO: implement outflow boundary conditions"
#if 0
    friend class BoxModel<TypeTag>;

    void evalBoundary_()
    {
        ParentType::evalBoundary_();

        if (this->bcTypes_().hasOutflow())
        {
            evalOutflow_();
        }

        if (this->bcTypes_().hasDirichlet())
            this->evalDirichlet_();
    }

    void evalOutflux(PrimaryVariables &outFlux,
                     const BoundaryFluxVariables &fluxVars,
                     const ElementContext &elemCtx,
                     int scvIdx,
                     int boundaryIdx) const
    {
        DUNE_THROW(Dune::NotImplemented, "evalOutflux()");
    };

    /*!
     * \brief Add all Outflow boundary conditions to the local
     *        residual.
     */
    void evalOutflow_()
    {
        Dune::GeometryType geoType = this->elem_().geometry().type();

        typedef typename Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
        typedef typename Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        IntersectionIterator isIt = this->gridView_().ibegin(this->elem_());
        const IntersectionIterator &endIt = this->gridView_().iend(this->elem_());
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all vertices of the current
            // face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElem.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 ++faceVertIdx)
            {
                int elemVertIdx = refElem.subEntity(faceIdx,
                                                    1,
                                                    faceVertIdx,
                                                    dim);

                int boundaryFaceIdx =
                    this->fvElemGeom_().boundaryFaceIndex(faceIdx, faceVertIdx);

                // add the residual of all vertices of the boundary
                // segment
                evalOutflowSegment_(isIt,
                                    elemVertIdx,
                                    boundaryFaceIdx);
            }
        }
    }

    /*!
     * \brief Add Outflow boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalOutflowSegment_(const IntersectionIterator &isIt,
                             int scvIdx,
                             int boundaryFaceIdx)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables flux(0.0);
        const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);

        // deal with neumann boundaries
        if (bcTypes.hasOutflow())
        {
            const BoundaryVariables boundaryVars(this->problem_(),
                                                 this->elem_(),
                                                 this->fvElemGeom_(),
                                                 boundaryFaceIdx,
                                                 this->curVolVars_(),
                                                 scvIdx);

            const VolumeVariables& vertVars = this->curVolVars_()[scvIdx];

            // mass balance
            if (bcTypes.isOutflow(contiEqIdx))
            {
                if(!useMolarFormulation) //use massfractions
                {
                    flux[contiEqIdx] += boundaryVars.KmvpNormal()*vertVars.density()/vertVars.viscosity();
                }
                else //use molefractions
                {
                    flux[contiEqIdx] += boundaryVars.KmvpNormal()*vertVars.molarDensity()/vertVars.viscosity();
                }
            }

            // component transport
            if (bcTypes.isOutflow(transEqIdx))
            {
                if(!useMolarFormulation)//use massfractions
                {
                    // advective flux
                    flux[transEqIdx]+= boundaryVars.KmvpNormal()*vertVars.density()/vertVars.viscosity()
                        *vertVars.fluidState().massFraction(phaseIdx, /*compIdx=*/1);

                    // diffusive flux of comp1 component in phase0
                    Scalar tmp = 0;
                    for (int i = 0; i < Vector::size; ++ i)
                        tmp += boundaryVars.massFracGrad(/*compIdx=*/1)[i]*boundaryVars.boundaryFace().normal[i];
                    tmp *= -1;

                    tmp *= boundaryVars.porousDiffCoeff()*boundaryVars.densityAtIP();
                    flux[transEqIdx] += tmp;//* FluidSystem::molarMass(/*compIdx=*/1);
                }
                else //use molefractions
                {
                    // advective flux
                    flux[transEqIdx]+= boundaryVars.KmvpNormal()*vertVars.molarDensity()/vertVars.viscosity()
                        *vertVars.fluidState().moleFraction(phaseIdx, /*compIdx=*/1);

                    // diffusive flux of comp1 component in phase0
                    Scalar tmp = 0;
                    for (int i = 0; i < Vector::size; ++ i)
                        tmp += boundaryVars.moleFracGrad(/*compIdx=*/1)[i]*boundaryVars.boundaryFace().normal[i];
                    tmp *= -1;

                    tmp *= boundaryVars.porousDiffCoeff()*boundaryVars.molarDensityAtIP();
                    flux[transEqIdx] += tmp;
                }
            }
            this->residual_[scvIdx] += flux;
        }
    }
#endif

    Implementation &asImp_()
    { return *static_cast<Implementation *> (this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *> (this); }
};

}

#endif
