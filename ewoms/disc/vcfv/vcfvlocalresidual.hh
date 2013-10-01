// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::VcfvLocalResidual
 */
#ifndef EWOMS_VCFV_LOCAL_RESIDUAL_HH
#define EWOMS_VCFV_LOCAL_RESIDUAL_HH

#include <opm/material/Valgrind.hpp>

#include <dune/istl/bvector.hh>
#include <dune/grid/common/geometry.hh>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include "vcfvproperties.hh"
#include "vcfvboundarycontext.hh"
#include "vcfvconstraintscontext.hh"

namespace Ewoms {

/*!
 * \ingroup VcfvModel
 *
 * \brief Element-wise caculation of the residual matrix for models
 *        based on the VCVF discretization .
 *
 * \copydetails Doxygen::typeTagTParam
 */
template<class TypeTag>
class VcfvLocalResidual
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::ctype CoordScalar;

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    typedef typename Dune::ReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<CoordScalar, dim> ReferenceElement;
#else
    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;
#endif

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

    typedef Ewoms::VcfvConstraintsContext<TypeTag> ConstraintsContext;
    typedef Ewoms::VcfvBoundaryContext<TypeTag> BoundaryContext;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

    // copying the local residual class is not a good idea
    VcfvLocalResidual(const VcfvLocalResidual &)
    {}

public:
    VcfvLocalResidual()
    { }

    ~VcfvLocalResidual()
    { }

    /*!
     * \brief Register all run-time parameters for the local residual.
     */
    static void registerParameters()
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
     * \copydetails Doxygen::vcfvScvIdxParam
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
     * \copydetails Doxygen::vcfvScvIdxParam
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
     * \copydetails Doxygen::problemParam
     * \copydetails Doxygen::elementParam
     */
    void eval(const Problem &problem, const Element &element)
    {
        ElementContext elemCtx(problem);
        elemCtx.updateAll(element);
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
     * \copydetails Doxygen::vcfvElemCtxParam
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
     * \copydetails Doxygen::residualParam
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::vcfvElemCtxParam
     */
    void eval(LocalBlockVector &residual,
              LocalBlockVector &storage,
              const ElementContext &elemCtx) const
    {
        residual = 0.0;
        storage = 0.0;

        // evaluate the flux terms
        asImp_().evalFluxes(residual, elemCtx, /*timeIdx=*/0);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemCtx.fvElemGeom(/*timeIdx=*/0).numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // HAVE_VALGRIND

        // evaluate the storage and the source terms
        asImp_().evalVolumeTerms_(residual, storage, elemCtx);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemCtx.fvElemGeom(/*timeIdx=*/0).numVertices; i++)
            Valgrind::CheckDefined(residual[i]);
#endif // !defined NDEBUG && HAVE_VALGRIND

        // evaluate the boundary conditions
        asImp_().evalBoundary_(residual, storage, elemCtx, /*timeIdx=*/0);

        // evaluate the constraint DOFs
        asImp_().evalConstraints_(residual, storage, elemCtx, /*timeIdx=*/0);

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
     * \copydetails Doxygen::vcfvElemCtxParam
     * \copydetails Doxygen::timeIdxParam
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
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::vcfvElemCtxParam
     * \copydetails Doxygen::timeIdxParam
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
     * \copydetails Doxygen::residualParam
     * \copydetails Doxygen::vcfvElemCtxParam
     * \copydetails Doxygen::timeIdxParam
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

    /////////////////////////////
    // The following methods _must_ be overloaded by the actual flow
    // models!
    /////////////////////////////

    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a finite sub-control volume.
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::vcfvScvCtxParams
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                   "Not implemented: The local residual " << Dune::className<Implementation>()
                   << " does not implement the required method 'computeStorage()'");
    };

    /*!
     * \brief Evaluates the total mass flux of all conservation
     *        quantities over a face of a sub-control volume.
     *
     * \copydetails Doxygen::areaFluxParam
     * \copydetails Doxygen::vcfvScvfCtxParams
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                  "Not implemented: The local residual " << Dune::className<Implementation>()
                  << " does not implement the required method 'computeFlux()'");
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \copydoc Doxygen::sourceParam
     * \copydoc Doxygen::vcfvScvCtxParams
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                  "Not implemented: The local residual " << Dune::className<Implementation>()
                  << " does not implement the required method 'computeSource()'");
    }

protected:
    /*!
     * \brief Evaluate the boundary conditions of an element.
     */
    void evalBoundary_(LocalBlockVector &residual,
                       LocalBlockVector &storage,
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
                          LocalBlockVector &storage,
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
                storage[scvIdx][eqIdx] = 0.0;
            };
        };
    }

    /*!
     * \brief Add the change in the storage terms and the source term
     *        to the local residual of all sub-control volumes of the
     *        current element.
     */
    void evalVolumeTerms_(LocalBlockVector &residual,
                          LocalBlockVector &storage,
                          const ElementContext &elemCtx) const
    {
        EqVector tmp(0), tmp2(0);
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

            storage[scvIdx] += tmp;
            residual[scvIdx] += tmp;

            // subtract the source term from the residual
            asImp_().computeSource(sourceRate, elemCtx, scvIdx, /*timeIdx=*/0);
            sourceRate *= scvVolume;
            residual[scvIdx] -= sourceRate;

            // make sure that only defined quantities were used
            // to calculate the residual.
            Valgrind::CheckDefined(storage[scvIdx]);
            Valgrind::CheckDefined(residual[scvIdx]);
        }
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
