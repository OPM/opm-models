// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Ewoms::VtkPhasePresenceModule
 */
#ifndef EWOMS_VTK_PHASE_PRESENCE_MODULE_HH
#define EWOMS_VTK_PHASE_PRESENCE_MODULE_HH

#include <ewoms/vtk/vtkbaseoutputmodule.hh>

#include <ewoms/common/parametersystem.hh>
#include <opm/core/utility/PropertySystem.hpp>

namespace Opm {
namespace Properties {
// create new type tag for the VTK primary variables output
NEW_TYPE_TAG(VtkPhasePresence);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePhasePresence);

SET_BOOL_PROP(VtkPhasePresence, VtkWritePhasePresence, false);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class VtkPhasePresenceModule : public VtkBaseOutputModule<TypeTag>
{
    typedef VtkBaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;
    typedef typename ParentType::ScalarBuffer ScalarBuffer;

    enum { dim = GridView::dimension };

public:
    VtkPhasePresenceModule(const Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePhasePresence, "Include the phase presence pseudo primary variable in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (phasePresenceOutput_()) this->resizeScalarBuffer_(phasePresence_);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        if (!phasePresenceOutput_())
            return;

        for (int i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            // calculate the phase presence
            int phasePresence = elemCtx.primaryVars(i, /*timeIdx=*/0).phasePresence();
            int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            phasePresence_[I] = phasePresence;
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (phasePresenceOutput_()) this->commitScalarBuffer_(writer, "phase presence", phasePresence_);
    }

private:
    static bool phasePresenceOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWritePhasePresence); }

    ScalarBuffer phasePresence_;
};

} // namespace Ewoms

#endif
