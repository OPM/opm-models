/*
  Copyright (C) 2010-2013 by Andreas Lauser

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
 * \copydoc Ewoms::VtkPrimaryVarsModule
 */
#ifndef EWOMS_VTK_PRIMARY_VARS_MODULE_HH
#define EWOMS_VTK_PRIMARY_VARS_MODULE_HH

#include <ewoms/io/baseoutputmodule.hh>
#include <ewoms/io/vtkmultiwriter.hh>

#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/propertysystem.hh>

namespace Ewoms {
namespace Properties {
// create new type tag for the VTK primary variables output
NEW_TYPE_TAG(VtkPrimaryVars);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePrimaryVars);
NEW_PROP_TAG(VtkWriteProcessRank);
NEW_PROP_TAG(VtkOutputFormat);

SET_BOOL_PROP(VtkPrimaryVars, VtkWritePrimaryVars, false);
SET_BOOL_PROP(VtkPrimaryVars, VtkWriteProcessRank, false);
} // namespace Properties

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class VtkPrimaryVarsModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Ewoms::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::EqBuffer EqBuffer;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:
    VtkPrimaryVarsModule(const Simulator &simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePrimaryVars,
                             "Include the primary variables in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteProcessRank,
                             "Include the MPI process rank in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (primaryVarsOutput_())
            this->resizeEqBuffer_(primaryVars_);
        if (processRankOutput_())
            this->resizeScalarBuffer_(processRank_,
                                      /*bufferType=*/ParentType::ElementBuffer);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        const auto &elementMapper = elemCtx.model().elementMapper();
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        int elemIdx = elementMapper.index(elemCtx.element());
#else
        int elemIdx = elementMapper .map(elemCtx.element());
#endif
        if (processRankOutput_() && !processRank_.empty())
            processRank_[elemIdx] = this->simulator_.gridView().comm().rank();
        for (int i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto &priVars = elemCtx.primaryVars(i, /*timeIdx=*/0);

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (primaryVarsOutput_() && !primaryVars_[eqIdx].empty())
                    primaryVars_[eqIdx][I] = priVars[eqIdx];
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter &baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter) {
            return;
        }

        if (primaryVarsOutput_())
            this->commitPriVarsBuffer_(baseWriter, "PV_%s", primaryVars_);
        if (processRankOutput_())
            this->commitScalarBuffer_(baseWriter,
                                      "process rank",
                                      processRank_,
                                      /*bufferType=*/ParentType::ElementBuffer);
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
