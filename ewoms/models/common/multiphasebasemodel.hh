// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2013 by Andreas Lauser

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
 * \copydoc Ewoms::MultiPhaseBaseModel
 */
#ifndef EWOMS_MULTI_PHASE_BASE_MODEL_HH
#define EWOMS_MULTI_PHASE_BASE_MODEL_HH

#include <ewoms/parallel/mpihelper.hh>
#include "multiphasebaseproperties.hh"

#include <dune/common/unused.hh>

namespace Ewoms {
/*!
 * \ingroup MultiPhaseBaseModel
 * \brief A base class for fully-implicit multi-phase porous-media flow models
 *        which assume multiple fluid phases.
 */
template <class TypeTag>
class MultiPhaseBaseModel : public GET_PROP_TYPE(TypeTag, Discretization)
{
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = FluidSystem::numComponents };

public:
    MultiPhaseBaseModel(Problem &problem)
        : ParentType(problem)
    {}

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        // register runtime parameters of the VTK output modules
        Ewoms::VtkMultiPhaseModule<TypeTag>::registerParameters();
        Ewoms::VtkTemperatureModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::phaseIdxParam
     */
    void globalPhaseStorage(EqVector &storage, int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        storage = 0;
        EqVector tmp;

        ElementContext elemCtx(this->problem_);
        ElementIterator elemIt = this->gridView_.template begin<0>();
        const ElementIterator elemEndIt = this->gridView_.template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            if (elemIt->partitionType() != Dune::InteriorEntity)
                continue; // ignore ghost and overlap elements

            elemCtx.updateStencil(*elemIt);
            elemCtx.updateVolVars(/*timeIdx=*/0);

            const auto &stencil = elemCtx.stencil(/*timeIdx=*/0);

            for (int dofIdx = 0; dofIdx < elemCtx.numDof(/*timeIdx=*/0); ++dofIdx) {
                const auto &scv = stencil.subControlVolume(dofIdx);
                const auto &volVars = elemCtx.volVars(dofIdx, /*timeIdx=*/0);

                tmp = 0;
                this->localResidual().addPhaseStorage(tmp,
                                                      elemCtx,
                                                      dofIdx,
                                                      /*timeIdx=*/0,
                                                      phaseIdx);
                tmp *= scv.volume()*volVars.extrusionFactor();
                storage += tmp;
            }
        };

        storage = this->gridView_.comm().sum(storage);
    }

    void registerVtkModules_()
    {
        ParentType::registerVtkModules_();

        // add the VTK output modules available on all model
        this->vtkOutputModules_.push_back(new Ewoms::VtkMultiPhaseModule<TypeTag>(this->problem_));
        this->vtkOutputModules_.push_back(new Ewoms::VtkTemperatureModule<TypeTag>(this->problem_));
    }

private:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};
} // namespace Ewoms

#include "multiphasebasepropertydefaults.hh"

#endif
