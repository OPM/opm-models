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
 * \copydoc Ewoms::BoxVtkPhasePresenceModule
 */
#ifndef EWOMS_BOX_VTK_PHASE_PRESENCE_MODULE_HH
#define EWOMS_BOX_VTK_PHASE_PRESENCE_MODULE_HH

#include "boxvtkoutputmodule.hh"

#include <ewoms/common/parameters.hh>
#include <ewoms/common/propertysystem.hh>

namespace Ewoms
{
namespace Properties
{
// create new type tag for the VTK primary variables output
NEW_TYPE_TAG(VtkPhasePresence);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePhasePresence);

SET_BOOL_PROP(VtkPhasePresence, VtkWritePhasePresence, false);
}

/*!
 * \ingroup BoxVtk
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class BoxVtkPhasePresenceModule : public BoxVtkOutputModule<TypeTag>
{
    typedef BoxVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;
    typedef typename ParentType::ScalarBuffer ScalarBuffer;

    enum { dim = GridView::dimension };

public:
    BoxVtkPhasePresenceModule(const Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, bool, VtkWritePhasePresence, "Include the phase presence pseudo primary variable in the VTK output files");
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

        const auto &vertexMapper = elemCtx.problem().vertexMapper();
        const auto &elem = elemCtx.element();
        for (int i = 0; i < elemCtx.numScv(); ++i) {
            // calculate the phase presence
            int phasePresence = elemCtx.primaryVars(i, /*timeIdx=*/0).phasePresence();
            int I = vertexMapper.map(elem, i, dim);
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
    { return GET_PARAM(TypeTag, bool, VtkWritePhasePresence); }

    ScalarBuffer phasePresence_;
};

}

#endif
