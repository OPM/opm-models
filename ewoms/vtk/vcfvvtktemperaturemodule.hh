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
 * \copydoc Ewoms::VcfvVtkTemperatureModule
 */
#ifndef EWOMS_VCFV_VTK_TEMPERATURE_MODULE_HH
#define EWOMS_VCFV_VTK_TEMPERATURE_MODULE_HH

#include "vcfvvtkoutputmodule.hh"

#include <ewoms/common/parametersystem.hh>
#include <opm/core/utility/PropertySystem.hpp>

namespace Opm {
namespace Properties {
// create new type tag for the VTK temperature output
NEW_TYPE_TAG(VtkTemperature);

// create the property tags needed for the temperature module
NEW_PROP_TAG(VtkWriteTemperature);
NEW_PROP_TAG(VtkWriteSolidHeatCapacity);
NEW_PROP_TAG(VtkWriteInternalEnergies);
NEW_PROP_TAG(VtkWriteEnthalpies);

// set default values for what quantities to output
SET_BOOL_PROP(VtkTemperature, VtkWriteTemperature, true);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup VcfvVtk
 *
 * \brief VTK output module for the temperature in which assume
 *        thermal equilibrium
 */
template<class TypeTag>
class VcfvVtkTemperatureModule : public VcfvVtkOutputModule<TypeTag>
{
    typedef VcfvVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

public:
    VcfvVtkTemperatureModule(const Problem &problem)
        : ParentType(problem)
    {
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteTemperature, "Include the temperature in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (temperatureOutput_()) this->resizeScalarBuffer_(temperature_);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        const auto &vertexMapper = elemCtx.problem().vertexMapper();
        const auto &elem = elemCtx.element();

        for (int i = 0; i < elemCtx.numScv(); ++i) {
            int I = vertexMapper.map(elem, i, dim);
            const auto &volVars = elemCtx.volVars(i, /*timeIdx=*/0);
            const auto &fs = volVars.fluidState();

            if (temperatureOutput_()) temperature_[I] = fs.temperature(/*phaseIdx=*/0);
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (temperatureOutput_()) this->commitScalarBuffer_(writer, "temperature", temperature_);
    }

private:
    static bool temperatureOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteTemperature); }

    ScalarBuffer temperature_;
};

} // namespace Ewoms

#endif
