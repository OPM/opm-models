// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::VcfvVtkPrimaryVarsModule
 */
#ifndef EWOMS_VCFV_VTK_PRIMARY_VARS_MODULE_HH
#define EWOMS_VCFV_VTK_PRIMARY_VARS_MODULE_HH

#include "vcfvvtkoutputmodule.hh"

#include <ewoms/common/parametersystem.hh>
#include <opm/core/utility/PropertySystem.hpp>

namespace Opm {
namespace Properties {
// create new type tag for the VTK primary variables output
NEW_TYPE_TAG(VtkPrimaryVars);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePrimaryVars);
NEW_PROP_TAG(VtkWriteProcessRank);

SET_BOOL_PROP(VtkPrimaryVars, VtkWritePrimaryVars, false);
SET_BOOL_PROP(VtkPrimaryVars, VtkWriteProcessRank, false);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvVtk
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class VcfvVtkPrimaryVarsModule : public VcfvVtkOutputModule<TypeTag>
{
    typedef VcfvVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::EqBuffer EqBuffer;

    enum { dim = GridView::dimension };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:
    VcfvVtkPrimaryVarsModule(const Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePrimaryVars, "Include the primary variables in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteProcessRank, "Include the MPI process rank in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (primaryVarsOutput_()) this->resizeEqBuffer_(primaryVars_);
        if (processRankOutput_()) this->resizeScalarBuffer_(processRank_, /*vertexCentered=*/false);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        const auto &elementMapper = elemCtx.problem().elementMapper();
        const auto &vertexMapper = elemCtx.problem().vertexMapper();
        const auto &elem = elemCtx.element();

        auto elemIdx = elementMapper.map(elemCtx.element());
        if (processRankOutput_()) processRank_[elemIdx] = this->problem_.gridView().comm().rank();
        for (int i = 0; i < elemCtx.numScv(); ++i) {

            int I = vertexMapper.map(elem, i, dim);
            const auto &priVars = elemCtx.primaryVars(i, /*timeIdx=*/0);

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (primaryVarsOutput_()) primaryVars_[eqIdx][I] = priVars[eqIdx];
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (primaryVarsOutput_()) this->commitPriVarsBuffer_(writer, "PV_%s", primaryVars_);
        if (processRankOutput_()) this->commitScalarBuffer_(writer, "process rank", processRank_, /*vertexCentered=*/false);
    }

private:
    static bool primaryVarsOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWritePrimaryVars); }
    static bool processRankOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteProcessRank); }

    EqBuffer primaryVars_;
    ScalarBuffer processRank_;
};

} // namespace Ewoms

#endif
