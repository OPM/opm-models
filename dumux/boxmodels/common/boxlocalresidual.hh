// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2011 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#include <dune/istl/bvector.hh>
#include <dune/grid/common/geometry.hh>

#include <dumux/common/valgrind.hh>

#include "boxproperties.hh"
#include "boxboundarycontext.hh"
#include "boxneumanncontext.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise caculation of the residual matrix for models
 *        based on the box scheme .
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual
{
private:
    typedef BoxLocalResidual<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::ctype CoordScalar;

    typedef typename Dune::GenericReferenceElements<CoordScalar, dim> ReferenceElements;
    typedef typename Dune::GenericReferenceElement<CoordScalar, dim> ReferenceElement;

    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Element::Geometry Geometry;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef Dumux::BoxNeumannContext<TypeTag> NeumannContext;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

    // copying the local residual class is not a good idea
    BoxLocalResidual(const BoxLocalResidual &)
    {};

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
     */
    void eval(const Problem &problem, const Element &elem)
    {
        ElementContext elemCtx(problem);
        elemCtx.updateAll(elem);
        eval(elemCtx);
    };

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        conservation equations from zero and store the results
     *        internally.
     *
     * The results can be requested afterwards using the residual()
     * and storageTerm() methods.
     */
    void eval(const ElementContext &elemCtx)
    {
        int numScv = elemCtx.numScv();
        internalResidual_.resize(numScv);
        internalStorageTerm_.resize(numScv);
        asImp_().eval(internalResidual_, internalStorageTerm_, elemCtx);
    };

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
        if (!elemCtx.onBoundary()) {
            return;
        }

        if (elemCtx.hasNeumann())
            asImp_().evalNeumann_(residual, elemCtx, timeIdx);

        if (elemCtx.hasDirichlet())
            asImp_().evalDirichlet_(residual, storageTerm, elemCtx, timeIdx);
    }

    /*!
     * \brief Set the values of the Dirichlet boundary control volumes
     *        of the current element.
     */
    void evalDirichlet_(LocalBlockVector &residual,
                        LocalBlockVector &storageTerm,
                        const ElementContext &elemCtx,
                        int timeIdx) const
    {
        PrimaryVariables tmp(0);
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); ++scvIdx) {
            const BoundaryTypes &bcTypes = elemCtx.boundaryTypes(scvIdx, timeIdx);
            if (!bcTypes.hasDirichlet())
                continue;

            // ask the problem for the dirichlet values
            Valgrind::SetUndefined(tmp);
            asImp_().computeDirichlet_(tmp, elemCtx, scvIdx, timeIdx);
            const PrimaryVariables &priVars =
                elemCtx.primaryVars(scvIdx, timeIdx);

            // set the dirichlet conditions
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!bcTypes.isDirichlet(eqIdx))
                    continue;
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);

                assert(0 <= pvIdx && pvIdx < numEq);
                Valgrind::CheckDefined(tmp[pvIdx]);

                residual[scvIdx][eqIdx] = priVars[pvIdx] - tmp[pvIdx];
                storageTerm[scvIdx][eqIdx] = 0.0;
            };
        };
    }

    /*!
     * \brief Add all Neumann boundary conditions to the local
     *        residual.
     */
    void evalNeumann_(LocalBlockVector &residual,
                      const ElementContext &elemCtx,
                      int timeIdx) const
    {
        const Element &elem = elemCtx.element();
        Dune::GeometryType geoType = elem.geometry().type();
        const ReferenceElement &refElem = ReferenceElements::general(geoType);

        NeumannContext neumannVars(elemCtx);
        const GridView &gridView = elemCtx.gridView();
        IntersectionIterator &isIt = neumannVars.intersectionIt();
        const IntersectionIterator &endIt = gridView.iend(elem);
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
                int scvIdx = refElem.subEntity(/*entityIdx=*/faceIdx,
                                               /*entityCodim=*/1,
                                               /*subEntityIdx=*/faceVertIdx,
                                               /*subEntityCodim=*/dim);

                int boundaryFaceIdx =
                    elemCtx.fvElemGeom(timeIdx).boundaryFaceIndex(faceIdx, faceVertIdx);

                // add the residual of all vertices of the boundary
                // segment
                evalNeumannSegment_(residual,
                                    neumannVars,
                                    scvIdx,
                                    boundaryFaceIdx,
                                    timeIdx);
            }
        }
    }

    void computeDirichlet_(PrimaryVariables &values,
                           const ElementContext &elemCtx,
                           int scvIdx,
                           int timeIdx) const
    { elemCtx.problem().dirichlet(values, elemCtx, scvIdx, timeIdx); }


    /*!
     * \brief Add Neumann boundary conditions for a single sub-control
     *        volume face to the local residual.
     */
    void evalNeumannSegment_(LocalBlockVector &residual,
                             const NeumannContext &neumannVars,
                             int scvIdx,
                             int boundaryFaceIdx,
                             int timeIdx) const
    {
        // temporary vector to store the neumann boundary fluxes
        const BoundaryTypes &bcTypes = neumannVars.elemCtx().boundaryTypes(scvIdx, timeIdx);
        RateVector values;

        // deal with neumann boundaries
        if (bcTypes.hasNeumann()) {
            Valgrind::SetUndefined(values);
            neumannVars.problem().neumann(values,
                                          neumannVars,
                                          boundaryFaceIdx,
                                          timeIdx);
            Valgrind::CheckDefined(values);

            values *=
                neumannVars.fvElemGeom(timeIdx).boundaryFace[boundaryFaceIdx].area
                * neumannVars.elemCtx().volVars(scvIdx, timeIdx).extrusionFactor();

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (bcTypes.isNeumann(eqIdx))
                    residual[scvIdx][eqIdx] += values[eqIdx];
            }
        }
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

            // mass balance within the element. this is the
            // \f$\frac{m}{\partial t}\f$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            asImp_().computeStorage(tmp2,
                                    elemCtx,
                                    scvIdx,
                                    /*timeIdx=*/1);
            asImp_().computeStorage(tmp,
                                    elemCtx,
                                    scvIdx,
                                    /*timeIdx=*/0);

            tmp -= tmp2;
            tmp *=
                elemCtx.fvElemGeom(/*timeIdx=*/0).subContVol[scvIdx].volume
                * extrusionFactor
                / elemCtx.problem().timeManager().timeStepSize();

            storageTerm[scvIdx] += tmp;
            residual[scvIdx] += tmp;

            // subtract the source term from the residual
            asImp_().computeSource(sourceRate, elemCtx, scvIdx, /*timeIdx=*/0);
            sourceRate *= elemCtx.fvElemGeom(/*timeIdx=*/0).subContVol[scvIdx].volume * extrusionFactor;
            residual[scvIdx] -= sourceRate;

            // make sure that only defined quantities were used
            // to calculate the residual.
            Valgrind::CheckDefined(storageTerm[scvIdx]);
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
