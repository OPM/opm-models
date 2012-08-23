// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
 * \ingroup BoxModel
 *
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#include <dumux/common/valgrind.hh>

#include <dune/istl/bvector.hh>
#include <dune/grid/common/geometry.hh>

#include <dune/common/fvector.hh>

#include "boxproperties.hh"
#include "boxboundarycontext.hh"
#include "boxconstraintscontext.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \ingroup BoxLocalResidual
 *
 * \brief Element-wise caculation of the residual matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::ctype CoordScalar;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dumux::BoxConstraintsContext<TypeTag> ConstraintsContext;
    typedef Dumux::BoxBoundaryContext<TypeTag> BoundaryContext;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

    // copying the local residual class is not a good idea
    BoxLocalResidual(const BoxLocalResidual &)
    {}

public:
    BoxLocalResidual()
    { }

    ~BoxLocalResidual()
    { }

    /*!
     * \brief Return the result of the eval() call using internal
     *        storage.
     */
    const LocalBlockVector &residual() const
    { return internalResidual_; }

    /*!
     * \brief Return the result of the eval() call using internal
     *        storage.
     *
     * \param scvIdx The index of the sub-control volume index for
     *               which the local residual is of interest.
     */
    const VectorBlock &residual(int scvIdx) const
    { return internalResidual_[scvIdx]; }

    /*!
     * \brief Return the storage term calculated using the last call
     *        to eval() using internal storage.
     */
    const LocalBlockVector &storageTerm() const
    { return internalStorageTerm_; }
  
    /*!
     * \brief Return the storage term calculated using the last call
     *        to eval() using internal storage.
     *
     * \param scvIdx The index of the sub-control volume index for
     *               which the local storage term is of interest.
     */
    const VectorBlock &storageTerm(int scvIdx) const
    { return internalStorageTerm_[scvIdx]; }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        conservation equations from zero and store the results
     *        internally.
     *
     * The results can be requested afterwards using the residual()
     * and storageTerm() methods.
     *
     * \param problem The problem which is to be solved.
     * \param elem The grid element for which the local
     *             residual should be calculated.
     */
    void eval(const Problem &problem, const Element &elem)
    {
        ElementContext elemCtx(problem);
        elemCtx.updateAll(elem);
        eval(elemCtx);
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        conservation equations from zero and store the results
     *        internally.
     *
     * The results can be requested afterwards using the residual()
     * and storageTerm() methods.
     *
     * \param elemCtx The element execution context for which the
     *                local residual should be calculated.
     */
    void eval(const ElementContext &elemCtx)
    {
        int numScv = elemCtx.numScv();
        internalResidual_.resize(numScv);
        internalStorageTerm_.resize(numScv);
        asImp_().eval(internalResidual_, internalStorageTerm_, elemCtx);
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        conservation equations from zero.
     *
     * \param residual Stores the residual vector
     * \param storageTerm Stores the derivative of the storage term to time
     * \param elemCtx All the element's secondary variables required to calculate the local residual
     */
    void eval(LocalBlockVector &residual,
              LocalBlockVector &storageTerm,
              const ElementContext &elemCtx) const
    {
        residual = 0.0;
        storageTerm = 0.0;

        // evaluate the flux terms
        asImp_().evalFluxes(residual, elemCtx, /*timeIdx=*/0);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemCtx.fvElemGeom(/*timeIdx=*/0).numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // HAVE_VALGRIND

        // evaluate the storage and the source terms
        asImp_().evalVolumeTerms_(residual, storageTerm, elemCtx);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemCtx.fvElemGeom(/*timeIdx=*/0).numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // !defined NDEBUG && HAVE_VALGRIND

        // evaluate the boundary conditions
        asImp_().evalBoundary_(residual, storageTerm, elemCtx, /*timeIdx=*/0);

        // evaluate the constraint DOFs
        asImp_().evalConstraints_(residual, storageTerm, elemCtx, /*timeIdx=*/0);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemCtx.fvElemGeom(/*timeIdx=*/0).numVertices; i++) {
            Valgrind::CheckDefined(residual[i]);
        }
#endif // HAVE_VALGRIND
    }

    /*!
     * \brief Calculate the amount of all conservation quantities
     *        stored in all element's sub-control volumes for a given
     *        history index.
     *
     * This is used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \param elemCtx The element execution context for which the
     *                local residual should be calculated.
     * \param timeIdx The index for time discretition for which the
     *                local storage term ought to be calculated.
     */
    void evalStorage(const ElementContext &elemCtx,
                     int timeIdx)
    {
        int numScv = elemCtx.numScv();
        internalStorageTerm_.resize(numScv);
        evalStorage(internalStorageTerm_,
                    elemCtx,
                    timeIdx);
    }

    /*!
     * \brief Calculate the amount of all conservation quantities
     *        stored in all element's sub-control volumes for a given
     *        history index.
     *
     * This is used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \param storage A Dune::BlockVector<EqVector> which stores the local
     *                storage term.
     * \param elemCtx The element execution context for which the
     *                local storage term should be calculated.
     * \param timeIdx The index for time discretition for which the
     *                local storage term ought to be calculated.
     */
    void evalStorage(LocalBlockVector &storage,
                     const ElementContext &elemCtx,
                     int timeIdx) const
    {
        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (int scvIdx=0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            storage[scvIdx] = 0.0;
            asImp_().computeStorage(storage[scvIdx],
                                    elemCtx,
                                    scvIdx,
                                    timeIdx);
            storage[scvIdx] *=
                elemCtx.fvElemGeom(timeIdx).subContVol[scvIdx].volume
                * elemCtx.volVars(scvIdx, timeIdx).extrusionFactor();
        }
    }

    /*!
     * \brief Add the flux term to a local residual.
     *
     * \param residual A Dune::BlockVector<EqVector> which stores the local
     *                residual.
     * \param elemCtx The element execution context for which the
     *                local residual should be calculated.
     * \param timeIdx The index for time discretition for which the
     *                local residual ought to be calculated.
     */
    void evalFluxes(LocalBlockVector &residual,
                    const ElementContext &elemCtx,
                    int timeIdx) const
    {
        RateVector flux;

        // calculate the mass flux over the sub-control volume faces
        for (int scvfIdx = 0;
             scvfIdx < elemCtx.fvElemGeom(timeIdx).numEdges;
             scvfIdx++)
        {
            int i = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].i;
            int j = elemCtx.fvElemGeom(timeIdx).subContVolFace[scvfIdx].j;

            Valgrind::SetUndefined(flux);
            asImp_().computeFlux(flux, /*context=*/elemCtx, scvfIdx, timeIdx);
            flux *= elemCtx.fluxVars(scvfIdx, timeIdx).extrusionFactor();
            Valgrind::CheckDefined(flux);

            // The balance equation for a finite volume is given by
            //
            // dStorage/dt = Flux + Source
            //
            // where the 'Flux' and the 'Source' terms represent the
            // mass per second which _ENTER_ the finite
            // volume. Re-arranging this, we get
            //
            // dStorage/dt - Source - Flux = 0
            //
            // Since the mass flux as calculated by computeFlux() goes
            // _OUT_ of sub-control volume i and _INTO_ sub-control
            // volume j, we need to add the flux to finite volume i
            // and subtract it from finite volume j
            residual[i] += flux;
            residual[j] -= flux;
        }
    }

protected:
    /*!
     * \brief Evaluate the boundary conditions of an element.
     */
    void evalBoundary_(LocalBlockVector &residual,
                       LocalBlockVector &storageTerm,
                       const ElementContext &elemCtx,
                       int timeIdx) const
    {
        if (!elemCtx.onBoundary())
            return;

        const Element &elem = elemCtx.element();
        Dune::GeometryType geoType = elem.geometry().type();
        const ReferenceElement &refElement = ReferenceElements::general(geoType);

        BoundaryContext boundaryCtx(elemCtx);
        const GridView &gridView = elemCtx.gridView();
        IntersectionIterator &isIt = boundaryCtx.intersectionIt();
        const IntersectionIterator &endIt = gridView.iend(elem);
        for (; isIt != endIt; ++isIt)
        {
            // handle only faces on the boundary
            if (!isIt->boundary())
                continue;

            // Assemble the boundary for all vertices of the current
            // face
            int faceIdx = isIt->indexInInside();
            int numFaceVerts = refElement.size(faceIdx, 1, dim);
            for (int faceVertIdx = 0;
                 faceVertIdx < numFaceVerts;
                 ++faceVertIdx)
            {
                int boundaryFaceIdx =
                    elemCtx.fvElemGeom(timeIdx).boundaryFaceIndex(faceIdx, faceVertIdx);

                // add the residual of all vertices of the boundary
                // segment
                evalBoundarySegment_(residual,
                                     boundaryCtx,
                                     boundaryFaceIdx,
                                     timeIdx);
            }
        }
    }

    /*!
     * \brief Evaluate all boundary conditions for a single
     *        sub-control volume face to the local residual.
     */
    void evalBoundarySegment_(LocalBlockVector &residual,
                              const BoundaryContext &boundaryCtx,
                              int boundaryFaceIdx,
                              int timeIdx) const
    {
        BoundaryRateVector values;

        Valgrind::SetUndefined(values);
        boundaryCtx.problem().boundary(values,
                                       boundaryCtx,
                                       boundaryFaceIdx,
                                       timeIdx);
        Valgrind::CheckDefined(values);

        int scvIdx = boundaryCtx.insideScvIndex(boundaryFaceIdx, timeIdx);
        values *=
            boundaryCtx.fvElemGeom(timeIdx).boundaryFace[boundaryFaceIdx].area
            * boundaryCtx.elemContext().volVars(scvIdx, timeIdx).extrusionFactor();
        
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            residual[scvIdx][eqIdx] += values[eqIdx];
        }
    }

    /*!
     * \brief Set the values of the constraint volumes of the current element.
     */
    void evalConstraints_(LocalBlockVector &residual,
                          LocalBlockVector &storageTerm,
                          const ElementContext &elemCtx,
                          int timeIdx) const
    {
        if (!GET_PROP_VALUE(TypeTag, EnableConstraints))
            return;

        const auto &problem = elemCtx.problem();
        Constraints constraints;
        ConstraintsContext constraintsCtx(elemCtx);
        for (int scvIdx = 0; scvIdx < constraintsCtx.numScv(); ++scvIdx) {
            // ask the problem for the constraint values
            constraints.reset();
            problem.constraints(constraints, elemCtx, scvIdx, timeIdx);

            if (!constraints.isConstraint())
                continue;

            // enforce the constraints
            const PrimaryVariables &priVars = elemCtx.primaryVars(scvIdx, timeIdx);
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!constraints.isConstraint(eqIdx))
                    continue;

                int pvIdx = constraints.eqToPvIndex(eqIdx);

                assert(0 <= pvIdx && pvIdx < numEq);
                Valgrind::CheckDefined(constraints[pvIdx]);

                residual[scvIdx][eqIdx] = priVars[pvIdx] - constraints[pvIdx];
                storageTerm[scvIdx][eqIdx] = 0.0;
            };
        };
    }

    /*!
     * \brief Add the change in the storage terms and the source term
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    void evalVolumeTerms_(LocalBlockVector &residual,
                          LocalBlockVector &storageTerm,
                          const ElementContext &elemCtx) const
    {
        EqVector tmp, tmp2;
        RateVector sourceRate;

        // evaluate the volume terms (storage + source terms)
        for (int scvIdx=0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            Scalar extrusionFactor =
                elemCtx.volVars(scvIdx, /*timeIdx=*/0).extrusionFactor();
            Scalar scvVolume = 
                elemCtx.fvElemGeom(/*timeIdx=*/0).subContVol[scvIdx].volume * extrusionFactor;

            // mass balance within the element. this is the
            // \f$\frac{m}{\partial t}\f$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            asImp_().computeStorage(tmp,
                                    elemCtx,
                                    scvIdx,
                                    /*timeIdx=*/0);
            asImp_().computeStorage(tmp2,
                                    elemCtx,
                                    scvIdx,
                                    /*timeIdx=*/1);

            tmp -= tmp2;
            tmp *= scvVolume / elemCtx.problem().timeManager().timeStepSize();

            storageTerm[scvIdx] += tmp;
            residual[scvIdx] += tmp;

            // subtract the source term from the residual
            asImp_().computeSource(sourceRate, elemCtx, scvIdx, /*timeIdx=*/0);
            sourceRate *= scvVolume;
            residual[scvIdx] -= sourceRate;

            // make sure that only defined quantities were used
            // to calculate the residual.
            Valgrind::CheckDefined(storageTerm[scvIdx]);
            Valgrind::CheckDefined(residual[scvIdx]);
        }
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param source The source/sink term [kg/m^3] in the sub control
     *               volume for each component
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, scvIdx, timeIdx);
        Valgrind::CheckDefined(source);
    }

private:
    Implementation &asImp_()
    {
      assert(static_cast<Implementation*>(this) != 0);
      return *static_cast<Implementation*>(this);
    }

    const Implementation &asImp_() const
    {
      assert(static_cast<const Implementation*>(this) != 0);
      return *static_cast<const Implementation*>(this);
    }

    LocalBlockVector internalResidual_;
    LocalBlockVector internalStorageTerm_;
};

}

#endif
