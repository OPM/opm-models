/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
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
 *
 * \brief VTK output module for quantities which make sense for models which
 *        assume thermal equilibrium.
 */
#ifndef DUMUX_BOX_VTK_ENERGY_MODULE_HH
#define DUMUX_BOX_VTK_ENERGY_MODULE_HH

#include "boxvtkoutputmodule.hh"

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
// create new type tag for the VTK energy output
NEW_TYPE_TAG(VtkEnergy);

// create the property tags needed for the energy module
NEW_PROP_TAG(VtkWriteSolidHeatCapacity);
NEW_PROP_TAG(VtkWriteInternalEnergies);
NEW_PROP_TAG(VtkWriteEnthalpies);

// set default values for what quantities to output
SET_BOOL_PROP(VtkEnergy, VtkWriteSolidHeatCapacity, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteInternalEnergies, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteEnthalpies, false);
}

/*!
 * \ingroup BoxModels
 *
 * \brief VTK output module for quantities which make sense for models which
 *        assume thermal equilibrium.
 *
 * This module deals with the following quantities:
 * - Specific enthalpy of all fluid phases
 * - Specific internal energy of all fluid phases
 * - Specific heat capacity of the solid phase
 */
template<class TypeTag>
class BoxVtkEnergyModule : public BoxVtkOutputModule<TypeTag>
{
    typedef BoxVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef Dumux::VtkMultiWriter<GridView> VtkMultiWriter;

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;

public:
    BoxVtkEnergyModule(const Problem &problem)
        : ParentType(problem)
    {
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (enthalpyOutput_()) this->resizePhaseBuffer_(enthalpy_);
        if (internalEnergyOutput_()) this->resizePhaseBuffer_(internalEnergy_);

        if (solidHeatCapacityOutput_()) this->resizeScalarBuffer_(solidHeatCapacity_);
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
            const auto &volVars = elemCtx.volVars(i);
            const auto &fs = volVars.fluidState();

            if (solidHeatCapacityOutput_()) solidHeatCapacity_[I] = volVars.heatCapacitySolid();

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (enthalpyOutput_()) enthalpy_[phaseIdx][I] = fs.enthalpy(phaseIdx);
                if (internalEnergyOutput_()) internalEnergy_[phaseIdx][I] = fs.internalEnergy(phaseIdx);
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (solidHeatCapacityOutput_()) this->commitScalarBuffer_(writer, "heatCapacitySolid", solidHeatCapacity_);

        if (enthalpyOutput_()) this->commitPhaseBuffer_(writer, "enthalpy_%s", enthalpy_);
        if (internalEnergyOutput_()) this->commitPhaseBuffer_(writer, "internalEnergy_%s", internalEnergy_);
    }

private:
    static bool solidHeatCapacityOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteSolidHeatCapacity); };

    static bool enthalpyOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteEnthalpies); };

    static bool internalEnergyOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteInternalEnergies); };

    PhaseBuffer enthalpy_;
    PhaseBuffer internalEnergy_;

    ScalarBuffer solidHeatCapacity_;
};

}

#endif
