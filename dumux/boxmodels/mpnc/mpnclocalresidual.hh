// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
#ifndef DUMUX_MPNC_LOCAL_RESIDUAL_HH
#define DUMUX_MPNC_LOCAL_RESIDUAL_HH

#include "mpncfluxvariables.hh"
#include "diffusion/diffusion.hh"
#include "energy/mpnclocalresidualenergy.hh"
#include "mass/mpnclocalresidualmass.hh"

#include <dumux/boxmodels/common/boxmodel.hh>

#include <dumux/common/math.hh>


namespace Dumux
{
/*!
 * \ingroup MPNCModel
 * \ingroup BoxLocalResidual
 * \brief two-phase, N-component specific details needed to
 *        approximately calculate the local defect in the BOX scheme.
 *
 * This class is used to fill the gaps in BoxLocalResidual for the
 * two-phase, N-component twophase flow.
 */
template<class TypeTag>
class MPNCLocalResidual : public BoxLocalResidual<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, MPNCIndices) Indices;

protected:
    typedef typename GET_PROP_TYPE(TypeTag, LocalResidual) Implementation;
    typedef BoxLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;


    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy),
        enableKineticEnergy = GET_PROP_VALUE(TypeTag, EnableKineticEnergy),

        enableKinetic = GET_PROP_VALUE(TypeTag, EnableKinetic),

        phase0NcpIdx = Indices::phase0NcpIdx
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;




    typedef MPNCLocalResidualEnergy<TypeTag, enableEnergy, enableKineticEnergy> EnergyResid;
    typedef MPNCLocalResidualMass<TypeTag, enableKinetic> MassResid;

    typedef Dune::FieldVector<Scalar, numEq> VectorBlock;
    typedef Dune::BlockVector<VectorBlock> LocalBlockVector;

public:
    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(PrimaryVariables &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int historyIdx) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VolumeVariables &volVars =
            elemCtx.volVars(scvIdx, historyIdx);

        storage = 0;

        // compute mass and energy storage terms
        MassResid::computeStorage(storage, volVars);
        Valgrind::CheckDefined(storage);
        EnergyResid::computeStorage(storage, volVars);
        Valgrind::CheckDefined(storage);
    }

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within all sub-control volumes of an
     *        element.
     */
    void addPhaseStorage(PrimaryVariables &storage,
                         const ElementContext &elemCtx,
                         int phaseIdx) const
    {
        // calculate the phase storage for all sub-control volumes
        for (int scvIdx=0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            PrimaryVariables tmp(0.0);

            // compute mass and energy storage terms in terms of
            // averaged quantities
            MassResid::addPhaseStorage(tmp,
                                       elemCtx.volVars(scvIdx),
                                       phaseIdx);
            EnergyResid::addPhaseStorage(tmp,
                                         elemCtx.volVars(scvIdx),
                                         phaseIdx);

            // multiply with volume of sub-control volume
            tmp *=
                elemCtx.volVars(scvIdx).extrusionFactor() *
                elemCtx.fvElemGeom().subContVol[scvIdx].volume;

            // Add the storage of the current SCV to the total storage
            storage += tmp;
        }
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(PrimaryVariables &source,
                       const ElementContext &elemCtx,
                       int scvIdx,
                       int historyIdx = 0) const
    {
        Valgrind::SetUndefined(source);
        elemCtx.problem().source(source, elemCtx, scvIdx);

        PrimaryVariables tmp(0);
        MassResid::computeSource(tmp, elemCtx, scvIdx);
        source += tmp;
        //EnergyResid::computeSource(tmp, elemCtx, scvIdx);
        Valgrind::CheckDefined(source);
    };

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(PrimaryVariables &flux,
                     const ElementContext &elemCtx,
                     int scvfIdx) const
    {
        flux = 0.0;
        MassResid::computeFlux(flux, elemCtx, scvfIdx);
        Valgrind::CheckDefined(flux);
        /*
         * EnergyResid also called in the MassResid
         * 1) Makes some sense because energy is also carried by mass
         * 2) The component-wise mass flux in each phase is needed.
         */
    }

    /*!
     * \brief Evaluate the local residual.
     */
    using ParentType::eval;
    void eval(LocalBlockVector &residual,
              LocalBlockVector &storageTerm,
              const ElementContext &elemCtx) const
    {
        ParentType::eval(residual,
                         storageTerm,
                         elemCtx);

        // handle the M additional model equations, make sure that
        // the dirichlet boundary condition is conserved
        int numScv = elemCtx.numScv();
        for (int scvIdx = 0; scvIdx < numScv; ++scvIdx) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                if (!elemCtx.boundaryTypes(scvIdx).isDirichlet(phase0NcpIdx + phaseIdx))
                    residual[scvIdx][phase0NcpIdx + phaseIdx] =
                        phaseNcp_(elemCtx, scvIdx, phaseIdx);
            }
        }
    }

protected:
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    /*!
     * \brief Returns the value of the NCP-function for a phase.
     */
    Scalar phaseNcp_(const ElementContext &elemCtx,
                     int scvIdx,
                     int phaseIdx) const
    {
        const auto &fsEval = elemCtx.evalPointVolVars(scvIdx).fluidState();
        const auto &fs = elemCtx.volVars(scvIdx).fluidState();

        Scalar aEval = phaseNotPresentIneq_(fsEval, phaseIdx);
        Scalar bEval = phasePresentIneq_(fsEval, phaseIdx);
        if (aEval > bEval)
            return phasePresentIneq_(fs, phaseIdx);
        return phaseNotPresentIneq_(fs, phaseIdx);
    };

    /*!
     * \brief Returns the value of the inequality where a phase is
     *        present.
     */
    template <class FluidState>
    Scalar phasePresentIneq_(const FluidState &fluidState, int phaseIdx) const
    { return fluidState.saturation(phaseIdx); }

    /*!
     * \brief Returns the value of the inequality where a phase is not
     *        present.
     */
    template <class FluidState>
    Scalar phaseNotPresentIneq_(const FluidState &fluidState, int phaseIdx) const
    {
        // difference of sum of mole fractions in the phase from 100%
        Scalar a = 1;
        for (int i = 0; i < numComponents; ++i)
            a -= fluidState.moleFraction(phaseIdx, i);
        return a;
    }
};

} // end namepace

#endif
