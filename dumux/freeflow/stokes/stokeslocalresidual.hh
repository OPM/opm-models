// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Katherina Baber, Klaus Mosthaf               *
 *   Copyright (C) 2008-2009 by Bernd Flemisch, Andreas Lauser               *
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
 * \brief The local residual function for problems using the
 *        Stokes box model.
 */

#ifndef DUMUX_STOKES_LOCAL_RESIDUAL_HH
#define DUMUX_STOKES_LOCAL_RESIDUAL_HH

#include "stokesvolumevariables.hh"
#include "stokesfluxvariables.hh"
#include "stokesproperties.hh"

#include <dumux/boxmodels/common/boxmodel.hh>
#include <dumux/boxmodels/common/boxneumanncontext.hh>

#include <dune/grid/common/grid.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dumux
{
/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxLocalResidual
 * \brief The local residual function for problems using the
 *        Stokes box model.
 */
template<class TypeTag>
class StokesLocalResidual : public GET_PROP_TYPE(TypeTag, BaseLocalResidual)
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, BaseLocalResidual) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum {
        dim = GridView::dimension,
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        phaseIdx = GET_PROP_VALUE(TypeTag, StokesPhaseIndex)
    };
    enum {
        massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance
        momentum0Idx = Indices::momentum0Idx //!< Index of the x-component of the momentum balance
    };
    enum { pressureIdx = Indices::pressureIdx }; //!< Index of the pressure in a solution vector

    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
    typedef Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;

    typedef Dune::FieldVector<Scalar, dim> ScalarGradient;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dumux::BoxNeumannContext<TypeTag> NeumannContext;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

 public:
    /*!
     * \brief Constructor. Sets the upwind weight and the stabilization parameters.
     */
    StokesLocalResidual()
    {
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        massUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
        stabilizationAlpha_ = GET_PARAM(TypeTag, Scalar, StabilizationAlpha);
        stabilizationBeta_ = GET_PARAM(TypeTag, Scalar, StabilizationBeta);
    };

    /*!
     * \brief Evaluates the amount of all conservation quantities
     *        (mass and momentum) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param result The mass of the component within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(EqVector &result,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);
        const auto &fs = volVars.fluidState();

        // mass balance
        result = 0.0;
        result[massBalanceIdx] = fs.density(phaseIdx);

        // momentum balance
        for (int dimIdx = 0; dimIdx < dim; ++ dimIdx)
            result[momentum0Idx + dimIdx] = 
                fs.density(phaseIdx)
                * volVars.velocity()[dimIdx];
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a sub-control volume. The face may be within
     *        an element (SCV face) or on the boundary. The advective and
     *        the diffusive fluxes are computed.
     *
     * \param flux The flux over the SCV (sub-control-volume) face
     * \param faceIdx The index of the SCV face (may also be a boundary face)
     * \param onBoundary Indicates, if the flux is evaluated on a boundary face. If it is true,
     *        the created fluxVars object contains boundary variables evaluated at the IP of the
     *        boundary face
     */
    void computeFlux(RateVector &flux, 
                     const ElementContext &elemCtx,
                     int faceIdx,
                     int timeIdx,
                     bool onBoundary=false) const
    {
        flux = 0.0;
        asImp_()->computeAdvectiveFlux(flux, elemCtx, faceIdx, timeIdx, onBoundary);
        Valgrind::CheckDefined(flux);
        asImp_()->computeDiffusiveFlux(flux, elemCtx, faceIdx, timeIdx, onBoundary);
        Valgrind::CheckDefined(flux);
    }

    /*!
     * \brief Evaluates the advective fluxes over
     *        a face of a sub-control volume.
     *
     * \param flux The advective flux over the sub-control-volume face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeAdvectiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int faceIdx,
                              int timeIdx,
                              bool onBoundary) const
    {
        const FluxVariables &fluxVars = 
            onBoundary
            ? elemCtx.boundaryFluxVars(faceIdx, timeIdx)
            : elemCtx.fluxVars(faceIdx, timeIdx);

        // if the momentum balance has a dirichlet b.c., the mass balance
        // is replaced, thus we do not need to calculate outflow fluxes here
        if (fluxVars.onBoundary())
            if (momentumBalanceDirichlet_(elemCtx.boundaryTypes(fluxVars.upstreamIdx(), timeIdx)))
                return;

        // data attached to upstream and the downstream vertices
        const VolumeVariables &up = elemCtx.volVars(fluxVars.upstreamIdx(), timeIdx);
        const VolumeVariables &dn = elemCtx.volVars(fluxVars.downstreamIdx(), timeIdx);

        const auto &fsUp = up.fluidState();
        const auto &fsDn = dn.fluidState();
        
        // mass balance with upwinded density
        FieldVector massBalanceResidual = fluxVars.velocityAtIP();
        massBalanceResidual *= 
            massUpwindWeight_ * fsUp.density(phaseIdx)
            +
            (1.-massUpwindWeight_) * fsDn.density(phaseIdx);

        if (!onBoundary)
        {
            // stabilization of the mass balance
            // with 0.5*alpha*(V_i + V_j)*grad P
            FieldVector stabilizationTerm = fluxVars.pressureGradAtIP();
            const auto &fvGeom = elemCtx.fvElemGeom(timeIdx);
            Scalar tmp = fvGeom.subContVol[fluxVars.upstreamIdx()].volume;
            tmp += fvGeom.subContVol[fluxVars.upstreamIdx()].volume;
            tmp /= 2;

            stabilizationTerm *= stabilizationAlpha_*tmp;
            massBalanceResidual += stabilizationTerm;
        }

        flux[massBalanceIdx] +=
            massBalanceResidual
            * fluxVars.normal();

        // momentum balance - pressure is evaluated as volume term
        // at the center of the SCV in computeSource
        // viscosity is upwinded

        // compute symmetrized gradient for the momentum flux:
        // mu (grad v + (grad v)^t)
        Dune::FieldMatrix<Scalar, dim, dim> symmVelGrad = fluxVars.velocityGradAtIP();
        for (int i=0; i<dim; ++i)
            for (int j=0; j<dim; ++j)
                symmVelGrad[i][j] += fluxVars.velocityGradAtIP()[j][i];

        FieldVector velGradComp(0.);
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
        {
            velGradComp = symmVelGrad[dimIdx];

            // TODO: dilatation term has to be accounted for in outflow, coupling, neumann
            //            velGradComp[dimIdx] += 2./3*fluxVars.velocityDivAtIP;
            velGradComp *= fluxVars.viscosityAtIP();

            flux[momentum0Idx + dimIdx] -=
                velGradComp*fluxVars.normal();

            // gravity is accounted for in computeSource; alternatively:
            //            Scalar gravityTerm = fluxVars.densityAtIP *
            //                    this->problem_().gravity()[dim-1] *
            //                    fluxVars.face().ipGlobal[dim-1]*
            //                    fluxVars.face().normal[dimIdx];
            //            flux[momentumXIdx + dimIdx] -=
            //                    gravityTerm;

        }
    }


    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        a SCV face or a boundary face.
     *
     * It doesn't do anything in the Stokes model but is used by the
     * transport and non-isothermal models to calculate diffusive and
     * conductive fluxes.
     *
     * \param flux The diffusive flux over the SCV face or boundary face for each component
     * \param fluxVars The flux variables at the current SCV/boundary face
     */
    void computeDiffusiveFlux(RateVector &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx,
                              int timeIdx,
                              bool onBoundary) const
    { }

    /*!
     * \brief Calculate the source term of all equations.
     *        The pressure gradient at the center of a SCV is computed
     *        and the gravity term evaluated.
     *
     * \param q The source/sink in the sub control volume for each component
     * \param localVertexIdx The local index of the sub-control volume
     */
    void computeSource(RateVector &q,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        const auto &volVars = elemCtx.volVars(scvIdx, timeIdx);
        const auto &fs = volVars.fluidState();

        // retrieve the source term intrinsic to the problem
        elemCtx.problem().source(q, elemCtx, scvIdx, timeIdx);

        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        const Scalar alphaH2 = 
            stabilizationAlpha_
            * elemCtx.fvElemGeom(timeIdx).subContVol[scvIdx].volume
            * elemCtx.volVars(scvIdx, timeIdx).extrusionFactor();
        q[massBalanceIdx] *= alphaH2; // stabilization of the source term

        // pressure gradient at the center of the SCV,
        // the pressure is discretized as volume term,
        // while -mu grad v is calculated in computeFlux
        ScalarGradient pressureGradAtSCVCenter(0.0);
        ScalarGradient grad(0.0);

        for (int vertexIdx = 0; vertexIdx < elemCtx.numScv(); vertexIdx++)
        {
            grad = elemCtx.fvElemGeom(timeIdx).subContVol[scvIdx].gradCenter[vertexIdx];
            Valgrind::CheckDefined(grad);
            grad *= elemCtx.volVars(vertexIdx, timeIdx).fluidState().pressure(phaseIdx);

            pressureGradAtSCVCenter += grad;
        }

        // add the component of the pressure gradient to the respective part
        // of the momentum equation and take the gravity term into account
        // signs are inverted, since q is subtracted
        const auto &g = elemCtx.problem().gravity(elemCtx, scvIdx, timeIdx);
        for (int dimIdx=0; dimIdx<dim; ++dimIdx)
        {
            q[momentum0Idx + dimIdx] -= pressureGradAtSCVCenter[dimIdx];
            q[momentum0Idx + dimIdx] += fs.density(phaseIdx)*g[dimIdx];
        }
    }

    /*!
     * \brief The Stokes model needs a modified treatment of the boundary conditions as
     *        the common box models
     */
    void evalBoundary_(LocalBlockVector &residual,
                       LocalBlockVector &storageTerm,
                       const ElementContext &elemCtx,
                       int timeIdx) const
    {
        assert(residual.size() == elemCtx.numScv());
        
        /*
        if (!elemCtx.onBoundary())
            return;
        */

        NeumannContext neumannCtx(elemCtx);
        const auto &elem = elemCtx.element();

        // loop over vertices of the element
        const ReferenceElement &refElem = ReferenceElements::general(elem.geometry().type());
        for (int vertexIdx = 0; vertexIdx < elemCtx.numScv(); vertexIdx++)
        {
            // consider only SCVs on the boundary
            if (elemCtx.fvElemGeom(timeIdx).subContVol[vertexIdx].inner)
                continue;

            // important at corners of the grid
            FieldVector momentumResidual(0.0);
            FieldVector averagedNormal(0.0);
            int numberOfOuterFaces = 0;
            // evaluate boundary conditions for the intersections of
            // the current element
            const auto &bcTypes = elemCtx.boundaryTypes(vertexIdx, timeIdx);

            const auto &gridView = elemCtx.gridView();
            auto &isIt = neumannCtx.intersectionIt();

            isIt = gridView.ibegin(elem);
            const auto &endIt = gridView.iend(elem);
            for (; isIt != endIt; ++isIt)
            {
                // handle only intersections on the boundary
                if (!isIt->boundary())
                    continue;

                // assemble the boundary for all vertices of the current face
                const int faceIdx = isIt->indexInInside();
                const int numFaceVertices = refElem.size(faceIdx, 1, dim);

                // loop over the single vertices on the current face
                for (int faceVertIdx = 0; faceVertIdx < numFaceVertices; ++faceVertIdx)
                {
                    // only evaluate, if we consider the same face vertex as in the outer
                    // loop over the element vertices
                    if (refElem.subEntity(faceIdx, 1, faceVertIdx, dim)
                        != vertexIdx)
                        continue;

                    const int boundaryFaceIdx = elemCtx.fvElemGeom(timeIdx).boundaryFaceIndex(faceIdx, faceVertIdx);
                    const auto &boundaryVars = neumannCtx.fluxVars(boundaryFaceIdx, timeIdx);

                    // the computed residual of the momentum equations is stored
                    // into momentumResidual for the replacement of the mass balance
                    // in case of Dirichlet conditions for the momentum balance;
                    // the fluxes at the boundary are added in the second step
                    if (momentumBalanceDirichlet_(bcTypes))
                    {
                        FieldVector muGradVelNormal(0.);
                        const FieldVector &boundaryFaceNormal =
                            boundaryVars.normal();

                        boundaryVars.velocityGradAtIP().umv(boundaryFaceNormal, muGradVelNormal);
                        muGradVelNormal *= boundaryVars.viscosityAtIP();

                        for (int i=0; i < elemCtx.numScv(); i++)
                            Valgrind::CheckDefined(residual[i]);

                        for (int dimIdx=0; dimIdx < dim; ++dimIdx)
                            momentumResidual[dimIdx] = residual[vertexIdx][momentum0Idx+dimIdx];

                        //Sign is right!!!: boundary flux: -mu grad v n
                        //but to compensate outernormal -> residual - (-mu grad v n)
                        momentumResidual += muGradVelNormal;
                        averagedNormal += boundaryFaceNormal;
                    }

                    // evaluate fluxes at a single boundary segment
                    asImp_()->evalNeumannSegment_(residual, neumannCtx, vertexIdx, boundaryFaceIdx, timeIdx);
                    asImp_()->evalOutflowSegment_(residual, neumannCtx, vertexIdx, boundaryFaceIdx, timeIdx);

                    // count the number of outer faces to determine, if we are on
                    // a corner point and if an interpolation should be done
                    numberOfOuterFaces++;
                } // end loop over face vertices
            } // end loop over intersections

            if(!bcTypes.isDirichlet(massBalanceIdx))
            {
                // replace the mass balance by the sum of the residua of the momentum balance
                if (momentumBalanceDirichlet_(bcTypes)) {
                    replaceMassbalanceResidual_(residual, momentumResidual, averagedNormal, vertexIdx);
                }
                else {
                    // de-stabilize (remove alpha*grad p - alpha div f
                    // from computeFlux on the boundary)
                    removeStabilizationAtBoundary_(residual, elemCtx, vertexIdx);
                }
            }
            if (numberOfOuterFaces == 2)
                interpolateCornerPoints_(residual, elemCtx, vertexIdx);
        } // end loop over element vertices

        // evaluate the dirichlet conditions of the element
        if (elemCtx.hasDirichlet())
            asImp_()->evalDirichlet_(residual, 
                                     storageTerm, 
                                     elemCtx,
                                     timeIdx);
    }

protected:
    /*!
     * \brief Evaluate and add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalNeumannSegment_(LocalBlockVector &residual,
                             const NeumannContext &neumannCtx,
                             int scvIdx,
                             int boundaryFaceIdx,
                             int timeIdx) const
    {
        const auto &bcTypes = neumannCtx.boundaryTypes(boundaryFaceIdx, timeIdx);
        const auto &is = neumannCtx.intersection(boundaryFaceIdx);

        if (bcTypes.hasNeumann())
        {
            // call evalNeumannSegment_() of the base class first
            ParentType::evalNeumannSegment_(residual, neumannCtx, scvIdx, boundaryFaceIdx, timeIdx);

            // temporary vector to store the Neumann boundary fluxes
            PrimaryVariables values(0.0);
            if (momentumBalanceHasNeumann_(bcTypes))
            {
                // Neumann BC of momentum equation needs special treatment
                // mathematically Neumann BC: p n - mu grad v n = q
                // boundary terms: -mu grad v n
                // implement q * A (from evalBoundarySegment) - p n(unity) A
                FieldVector pressureCorrection(neumannCtx.fluxVars(boundaryFaceIdx, timeIdx).normal());
                const auto &volVars = neumannCtx.volVars(boundaryFaceIdx, timeIdx);
                const auto &boundaryVars = neumannCtx.fluxVars(boundaryFaceIdx, timeIdx);
                pressureCorrection *= volVars.fluidState().pressure(phaseIdx);
                for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                    if(bcTypes.isNeumann(momentum0Idx + dimIdx))
                        residual[scvIdx][momentum0Idx  + dimIdx] += pressureCorrection[momentum0Idx + dimIdx];

                // beta-stabilization at the boundary
                // in case of neumann conditions for the momentum equation;
                // calculate  mu grad v t t
                // center in the face of the reference element
                FieldVector tangent;
                if (stabilizationBeta_ != 0)
                {
                    const FieldVector& elementUnitNormal = is.centerUnitOuterNormal();
                    tangent[0] = elementUnitNormal[1];  //TODO: 3D
                    tangent[1] = -elementUnitNormal[0];
                    FieldVector tangentialVelGrad;
                    boundaryVars.velocityGradAtIP().mv(tangent, tangentialVelGrad);
                    tangentialVelGrad *= boundaryVars.viscosityAtIP();

                    residual[scvIdx][massBalanceIdx] -=
                        stabilizationBeta_*0.5*
                        volVars.fluidState().pressure(phaseIdx);
                    residual[scvIdx][massBalanceIdx] -= 
                        stabilizationBeta_*0.5*
                        (tangentialVelGrad*tangent);

                    for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                        residual[scvIdx][massBalanceIdx] -= stabilizationBeta_*0.5
                            * values[momentum0Idx + dimIdx]*elementUnitNormal[dimIdx];
                }
                Valgrind::CheckDefined(residual);
            }
        }
    }

    /*!
     * \brief Evaluate outflow boundary conditions for a single SCV face on the boundary.
     */
    void evalOutflowSegment_(LocalBlockVector &residual,
                             const NeumannContext &neumannCtx,
                             int scvIdx,
                             int boundaryFaceIdx,
                             int timeIdx) const
    {
        const auto &bcTypes = neumannCtx.boundaryTypes(boundaryFaceIdx, timeIdx);
        const auto &is = neumannCtx.intersection(boundaryFaceIdx);

        if (bcTypes.hasOutflow())
        {
            RateVector values(0.0);
            const auto &boundaryVars = neumannCtx.fluxVars(boundaryFaceIdx, timeIdx);

            asImp_()->computeFlux(values, neumannCtx.elemContext(), boundaryFaceIdx, timeIdx, /*onBoundary=*/true);
            Valgrind::CheckDefined(values);

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                if (!bcTypes.isOutflow(eqIdx) )
                    continue;
                // do not calculate outflow for the mass balance
                // if the momentum balance is dirichlet -
                // it is replaced in that case
                if (eqIdx==massBalanceIdx && momentumBalanceDirichlet_(bcTypes))
                    continue;
                // deduce outflow
                residual[scvIdx][eqIdx] += values[eqIdx];
            }

            // beta-stabilization at the boundary in case of outflow condition
            // for the momentum balance
            if(momentumBalanceOutflow_(bcTypes) && stabilizationBeta_ != 0)
            {
                // calculate  mu grad v t t for beta-stabilization
                // center in the face of the reference element
                FieldVector tangent;
                const FieldVector& elementUnitNormal = is.centerUnitOuterNormal();
                tangent[0] = elementUnitNormal[1];
                tangent[1] = -elementUnitNormal[0];
                FieldVector tangentialVelGrad;
                boundaryVars.velocityGradAtIP().mv(tangent, tangentialVelGrad);

                residual[scvIdx][massBalanceIdx] -= 0.5*stabilizationBeta_
                    * boundaryVars.viscosityAtIP()
                    * (tangentialVelGrad*tangent);
            }
        }
    }

    /*!
     * \brief Remove the alpha stabilization at boundaries.
     */
    void removeStabilizationAtBoundary_(LocalBlockVector &residual,
                                        const ElementContext &elemCtx, 
                                        int vertexIdx) const
    {
        if (stabilizationAlpha_ != 0)
        {
            FluxVariables fluxVars;
            const auto &fvElemGeom = elemCtx.fvElemGeom(/*timeIdx=*/0);

            // loop over the edges of the element
            for (int faceIdx = 0; faceIdx < fvElemGeom.numEdges; faceIdx++)
            {
                const auto &scvf = fvElemGeom.subContVolFace[faceIdx];
                
                const int i = scvf.i;
                const int j = scvf.j;

                if (i != vertexIdx && j != vertexIdx)
                    continue;
                fluxVars.update(elemCtx, faceIdx, /*timeIdx=*/0);

                Scalar tmp = 0.5*(fvElemGeom.subContVol[i].volume + fvElemGeom.subContVol[j].volume);
                const Scalar alphaH2 = stabilizationAlpha_*tmp;
                Scalar stabilizationTerm = 
                    fluxVars.pressureGradAtIP()
                    * fluxVars.normal()
                    * alphaH2;
                
                if (vertexIdx == i)
                    residual[i][massBalanceIdx] += stabilizationTerm;
                if (vertexIdx == j)
                    residual[j][massBalanceIdx] -= stabilizationTerm;
            }

            // destabilize source term
            RateVector q(0.0);
            elemCtx.problem().source(q, elemCtx, vertexIdx, /*timeIdx=*/0);

            const Scalar alphaH2 = 
                stabilizationAlpha_ 
                * elemCtx.fvElemGeom(/*timeIdx=*/0).subContVol[vertexIdx].volume;
            
            residual[vertexIdx][massBalanceIdx] += 
                alphaH2
                * q[massBalanceIdx]
                * elemCtx.fvElemGeom(/*timeIdx=*/0).subContVol[vertexIdx].volume
                * elemCtx.volVars(vertexIdx, /*timeIdx=*/0).extrusionFactor();
        }
    }

    /*!
     * \brief Interpolate the pressure at corner points of the grid, thus taking the degree of freedom there. 
     * 		  This is required due to stability reasons.
     */
    void interpolateCornerPoints_(LocalBlockVector &residual, 
                                  const ElementContext &elemCtx,
                                  int vertexIdx) const
    {
        const BoundaryTypes &bcTypes = elemCtx.boundaryTypes(vertexIdx, /*timeIdx=*/0);

        // TODO: 3D
        if (bcTypes.isCouplingInflow(massBalanceIdx) || bcTypes.isCouplingOutflow(massBalanceIdx))
        {
            if (vertexIdx == 0 || vertexIdx == 3)
                residual[vertexIdx][massBalanceIdx] =
                    elemCtx.primaryVars(/*scvIdx=*/0, /*timeIdx=*/0)[pressureIdx]
                    - elemCtx.primaryVars(/*scvIdx=*/3, /*timeIdx=*/0)[pressureIdx];
            if (vertexIdx == 1 || vertexIdx == 2)
                residual[vertexIdx][massBalanceIdx] =
                    elemCtx.primaryVars(/*scvIdx=*/1, /*timeIdx=*/0)[pressureIdx]
                    - elemCtx.primaryVars(/*scvIdx=*/2, /*timeIdx=*/0)[pressureIdx];
        }
        else
        {
            if (!bcTypes.isDirichlet(massBalanceIdx)) // do nothing in case of dirichlet
                residual[vertexIdx][massBalanceIdx] =
                    elemCtx.primaryVars(/*scvIdx=*/0, /*timeIdx=*/0)[pressureIdx]
                    + elemCtx.primaryVars(/*scvIdx=*/3, /*timeIdx=*/0)[pressureIdx]
                    - elemCtx.primaryVars(/*scvIdx=*/1, /*timeIdx=*/0)[pressureIdx]
                    - elemCtx.primaryVars(/*scvIdx=*/2, /*timeIdx=*/0)[pressureIdx];
        }
    }

    /*!
     * \brief Replace the local residual of the mass balance equation by
     *        the sum of the residuals of the momentum balance equation.
     */
    void replaceMassbalanceResidual_(LocalBlockVector &residual,
                                     const FieldVector& momentumResidual,
                                     FieldVector& averagedNormal,
                                     const int vertexIdx) const
    {
        assert(averagedNormal.two_norm() != 0.0);

        // divide averagedNormal by its length
        averagedNormal /= averagedNormal.two_norm();
        // replace the mass balance by the sum of the residuals of the momentum balances
        residual[vertexIdx][massBalanceIdx] = momentumResidual*averagedNormal;
    }

    /*!
     * \brief Returns true, if all boundary conditions for the momentum balance
     *        at the considered vertex are Dirichlet.
     */
    bool momentumBalanceDirichlet_(const BoundaryTypes& bcTypes) const
    {
        for (int dimIdx=0; dimIdx < dim; ++dimIdx)
            if (!bcTypes.isDirichlet(momentum0Idx + dimIdx))
                return false;
        return true;
    }

    /*!
     * \brief Returns true, if at least one boundary condition of the momentum balance is Neumann.
     */
    bool momentumBalanceHasNeumann_(const BoundaryTypes& bcTypes) const
    {
        for (int dimIdx=0; dimIdx < dim; ++dimIdx)
            if (bcTypes.isNeumann(momentum0Idx + dimIdx))
                return true;
        return false;
    }

    /*!
     * \brief Returns true, if all boundary conditions for the momentum balance are outflow.
     */
    bool momentumBalanceOutflow_(const BoundaryTypes& bcTypes) const
    {
        for (int dimIdx=0; dimIdx < dim; ++dimIdx)
            if (!bcTypes.isOutflow(momentum0Idx + dimIdx))
                return false;
        return true;
    }

    Scalar massUpwindWeight_;
    Scalar stabilizationAlpha_;
    Scalar stabilizationBeta_;

    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }
};

}

#endif
