/*
  Copyright (C) 2008-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::FvBaseLocalResidual
 */
#ifndef EWOMS_FV_BASE_LOCAL_RESIDUAL_HH
#define EWOMS_FV_BASE_LOCAL_RESIDUAL_HH

#include <opm/material/Valgrind.hpp>

#include <dune/istl/bvector.hh>
#include <dune/grid/common/geometry.hh>

#include <dune/common/fvector.hh>

#include <opm/core/utility/ClassName.hpp>

#include "fvbaseproperties.hh"

#include <cmath>

namespace Ewoms {

/*!
 * \ingroup Discretization
 *
 * \brief Element-wise caculation of the residual matrix for models
 *        based on a finite volume spatial discretization.
 *
 * \copydetails Doxygen::typeTagTParam
 */
template<class TypeTag>
class FvBaseLocalResidual
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dim = GridView::dimension };
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryContext) BoundaryContext;
    typedef typename GET_PROP_TYPE(TypeTag, ConstraintsContext) ConstraintsContext;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

    // copying the local residual class is not a good idea
    FvBaseLocalResidual(const FvBaseLocalResidual &)
    {}

public:
    FvBaseLocalResidual()
    { }

    ~FvBaseLocalResidual()
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
     * \copydetails Doxygen::ecfvScvIdxParam
     */
    const VectorBlock &residual(int dofIdx) const
    { return internalResidual_[dofIdx]; }

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
     * \copydetails Doxygen::ecfvScvIdxParam
     */
    const VectorBlock &storageTerm(int dofIdx) const
    { return internalStorageTerm_[dofIdx]; }

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
     * \copydetails Doxygen::ecfvElemCtxParam
     */
    void eval(const ElementContext &elemCtx)
    {
        int numPrimaryDof = elemCtx.numPrimaryDof(/*timeIdx=*/0);
        internalResidual_.resize(numPrimaryDof);
        internalStorageTerm_.resize(numPrimaryDof);
        asImp_().eval(internalResidual_, internalStorageTerm_, elemCtx);
    }

    /*!
     * \brief Compute the local residual, i.e. the deviation of the
     *        conservation equations from zero.
     *
     * \copydetails Doxygen::residualParam
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::ecfvElemCtxParam
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
        for (int i=0; i < elemCtx.stencil(/*timeIdx=*/0).numPrimaryDof(); i++) {
            for (int j = 0; j < numEq; ++ j)
                assert(std::isfinite(residual[i][j]));
            Valgrind::CheckDefined(residual[i]);
        }
#endif // !defined NDEBUG && HAVE_VALGRIND

        // evaluate the storage and the source terms
        asImp_().evalVolumeTerms_(residual, storage, elemCtx);

        for (int dofIdx=0; dofIdx < elemCtx.stencil(/*timeIdx=*/0).numPrimaryDof(); dofIdx++) {
            storage[dofIdx] /= elemCtx.dofTotalVolume(dofIdx, /*timeIdx=*/0);

#if !defined NDEBUG && HAVE_VALGRIND
            for (int j = 0; j < numEq; ++ j)
                assert(std::isfinite(residual[dofIdx][j]));
            Valgrind::CheckDefined(residual[dofIdx]);
#endif // !defined NDEBUG && HAVE_VALGRIND
        }

        // evaluate the boundary conditions
        asImp_().evalBoundary_(residual, storage, elemCtx, /*timeIdx=*/0);

#if !defined NDEBUG && HAVE_VALGRIND
        for (int i=0; i < elemCtx.stencil(/*timeIdx=*/0).numPrimaryDof(); i++) {
            for (int j = 0; j < numEq; ++ j)
                assert(std::isfinite(residual[i][j]));
            Valgrind::CheckDefined(residual[i]);
        }
#endif // !defined NDEBUG && HAVE_VALGRIND

        // evaluate the constraint DOFs
        asImp_().evalConstraints_(residual, storage, elemCtx, /*timeIdx=*/0);

        // make the residual volume specific (i.e., make it incorrect
        // mass per cubic meter instead of total mass)
        for (int dofIdx=0;
             dofIdx < elemCtx.stencil(/*timeIdx=*/0).numPrimaryDof();
             ++dofIdx)
        {
            assert(elemCtx.dofTotalVolume(dofIdx, /*timeIdx=*/0) > 0);
            residual[dofIdx] /= elemCtx.dofTotalVolume(dofIdx, /*timeIdx=*/0);
            assert(std::isfinite(residual[dofIdx].two_norm()));

#if !defined NDEBUG && HAVE_VALGRIND
            for (int j = 0; j < numEq; ++ j)
                assert(std::isfinite(residual[dofIdx][j]));
            Valgrind::CheckDefined(residual[dofIdx]);
#endif // !defined NDEBUG && HAVE_VALGRIND
        }
    }

    /*!
     * \brief Calculate the amount of all conservation quantities
     *        stored in all element's sub-control volumes for a given
     *        history index.
     *
     * This is used to figure out how much of each conservation
     * quantity is inside the element.
     *
     * \copydetails Doxygen::ecfvElemCtxParam
     * \copydetails Doxygen::timeIdxParam
     */
    void evalStorage(const ElementContext &elemCtx,
                     int timeIdx)
    {
        int numDof = elemCtx.numDof(/*timeIdx=*/0);
        internalStorageTerm_.resize(numDof);
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
     * \copydetails Doxygen::ecfvElemCtxParam
     * \copydetails Doxygen::timeIdxParam
     */
    void evalStorage(LocalBlockVector &storage,
                     const ElementContext &elemCtx,
                     int timeIdx) const
    {
        // calculate the amount of conservation each quantity inside
        // all sub control volumes
        for (int dofIdx=0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); dofIdx++)
        {
            storage[dofIdx] = 0.0;
            asImp_().computeStorage(storage[dofIdx],
                                    elemCtx,
                                    dofIdx,
                                    timeIdx);
            storage[dofIdx] *=
                elemCtx.stencil(timeIdx).subControlVolume(dofIdx).volume()
                * elemCtx.volVars(dofIdx, timeIdx).extrusionFactor();
        }
    }

    /*!
     * \brief Add the flux term to a local residual.
     *
     * \copydetails Doxygen::residualParam
     * \copydetails Doxygen::ecfvElemCtxParam
     * \copydetails Doxygen::timeIdxParam
     */
    void evalFluxes(LocalBlockVector &residual,
                    const ElementContext &elemCtx,
                    int timeIdx) const
    {
        RateVector flux;

        const auto &stencil = elemCtx.stencil(timeIdx);
        // calculate the mass flux over the sub-control volume faces
        for (int scvfIdx = 0; scvfIdx < elemCtx.stencil(timeIdx).numInteriorFaces(); scvfIdx++)
        {
            const auto &face = stencil.interiorFace(scvfIdx);
            int i = face.interiorIndex();
            int j = face.exteriorIndex();

            Valgrind::SetUndefined(flux);
            asImp_().computeFlux(flux, /*context=*/elemCtx, scvfIdx, timeIdx);
            flux *= elemCtx.fluxVars(scvfIdx, timeIdx).extrusionFactor()*face.area();
            Valgrind::CheckDefined(flux);

            // The balance equation for a finite volume is given by
            //
            // dStorage/dt + Flux = Source
            //
            // where the 'Flux' and the 'Source' terms represent the
            // mass per second which leaves the finite
            // volume. Re-arranging this, we get
            //
            // dStorage/dt + Flux - Source = 0
            //
            // Since the mass flux as calculated by computeFlux() goes
            // out of sub-control volume i and into sub-control volume
            // j, we need to add the flux to finite volume i and
            // subtract it from finite volume j
            if (i < stencil.numPrimaryDof())
                residual[i] += flux;
            if (j < stencil.numPrimaryDof())
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
     * \copydetails Doxygen::ecfvScvCtxParams
     */
    void computeStorage(EqVector &storage,
                        const ElementContext &elemCtx,
                        int dofIdx,
                        int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                   "Not implemented: The local residual " << Opm::className<Implementation>()
                   << " does not implement the required method 'computeStorage()'");
    }

    /*!
     * \brief Evaluates the total mass flux of all conservation
     *        quantities over a face of a sub-control volume.
     *
     * \copydetails Doxygen::areaFluxParam
     * \copydetails Doxygen::ecfvScvfCtxParams
     */
    void computeFlux(RateVector &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx,
                     int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                  "Not implemented: The local residual " << Opm::className<Implementation>()
                  << " does not implement the required method 'computeFlux()'");
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \copydoc Doxygen::sourceParam
     * \copydoc Doxygen::ecfvScvCtxParams
     */
    void computeSource(RateVector &source,
                       const ElementContext &elemCtx,
                       int dofIdx,
                       int timeIdx) const
    {
        OPM_THROW(std::logic_error,
                  "Not implemented: The local residual " << Opm::className<Implementation>()
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

        const auto &stencil = elemCtx.stencil(timeIdx);
        BoundaryContext boundaryCtx(elemCtx);

        // Assemble the boundary for all vertices of the current
        // face
        for (int faceIdx = 0; faceIdx < stencil.numBoundaryFaces(); ++faceIdx)
        {
            // add the residual of all vertices of the boundary
            // segment
            evalBoundarySegment_(residual,
                                 boundaryCtx,
                                 faceIdx,
                                 timeIdx);
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

        const auto &stencil = boundaryCtx.stencil(timeIdx);
        int dofIdx = stencil.boundaryFace(boundaryFaceIdx).interiorIndex();
        const auto &insideVolVars = boundaryCtx.elementContext().volVars(dofIdx, timeIdx);
        values *=
            stencil.boundaryFace(boundaryFaceIdx).area()
            * insideVolVars.extrusionFactor();

        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            residual[dofIdx][eqIdx] += values[eqIdx];
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
        for (int dofIdx = 0; dofIdx < constraintsCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
            // ask the problem for the constraint values
            constraints.reset();
            problem.constraints(constraints, elemCtx, dofIdx, timeIdx);

            if (!constraints.isConstraint())
                continue;

            // enforce the constraints
            const PrimaryVariables &priVars = elemCtx.primaryVars(dofIdx, timeIdx);
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (!constraints.isConstraint(eqIdx))
                    continue;

                int pvIdx = constraints.eqToPvIndex(eqIdx);

                assert(0 <= pvIdx && pvIdx < numEq);
                Valgrind::CheckDefined(constraints[pvIdx]);

                residual[dofIdx][eqIdx] = priVars[pvIdx] - constraints[pvIdx];
                storage[dofIdx][eqIdx] = 0.0;
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
        for (int dofIdx=0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); dofIdx++)
        {
            Scalar extrusionFactor =
                elemCtx.volVars(dofIdx, /*timeIdx=*/0).extrusionFactor();
            Scalar scvVolume =
                elemCtx.stencil(/*timeIdx=*/0).subControlVolume(dofIdx).volume() * extrusionFactor;

            // mass balance within the element. this is the
            // \f$\frac{m}{\partial t}\f$ term if using implicit
            // euler as time discretization.
            //
            // TODO (?): we might need a more explicit way for
            // doing the time discretization...
            asImp_().computeStorage(tmp,
                                    elemCtx,
                                    dofIdx,
                                    /*timeIdx=*/0);
            asImp_().computeStorage(tmp2,
                                    elemCtx,
                                    dofIdx,
                                    /*timeIdx=*/1);

            tmp -= tmp2;
            tmp *= scvVolume / elemCtx.problem().timeManager().timeStepSize();

            storage[dofIdx] += tmp;
            residual[dofIdx] += tmp;

            // subtract the source term from the residual
            asImp_().computeSource(sourceRate, elemCtx, dofIdx, /*timeIdx=*/0);
            sourceRate *= scvVolume;
            residual[dofIdx] -= sourceRate;

            // make sure that only defined quantities were used
            // to calculate the residual.
            Valgrind::CheckDefined(storage[dofIdx]);
            Valgrind::CheckDefined(residual[dofIdx]);
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

} // namespace Ewoms

#endif
