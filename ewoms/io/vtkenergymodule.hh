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
 * \copydoc Ewoms::VtkEnergyModule
 */
#ifndef EWOMS_VTK_ENERGY_MODULE_HH
#define EWOMS_VTK_ENERGY_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <opm/material/common/MathToolbox.hpp>

namespace Ewoms {
namespace Properties {
// create new type tag for the VTK energy output
NEW_TYPE_TAG(VtkEnergy);

// create the property tags needed for the energy module
NEW_PROP_TAG(VtkWriteSolidHeatCapacity);
NEW_PROP_TAG(VtkWriteHeatConductivity);
NEW_PROP_TAG(VtkWriteInternalEnergies);
NEW_PROP_TAG(VtkWriteEnthalpies);
NEW_PROP_TAG(VtkOutputFormat);

// set default values for what quantities to output
SET_BOOL_PROP(VtkEnergy, VtkWriteSolidHeatCapacity, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteHeatConductivity, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteInternalEnergies, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteEnthalpies, false);
} // namespace Properties
} // namespace Ewoms

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for quantities which make sense for models which
 *        assume thermal equilibrium.
 *
 * This module deals with the following quantities:
 * - Specific enthalpy of all fluid phases
 * - Specific internal energy of all fluid phases
 * - Specific heat capacity of the solid phase
 * - Lumped heat conductivity (solid phase plus all fluid phases)
 */
template <class TypeTag>
class VtkEnergyModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename Opm::MathToolbox<Evaluation> Toolbox;
    typedef Ewoms::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

public:
    VtkEnergyModule(const Simulator &simulator)
        : ParentType(simulator)
    {
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSolidHeatCapacity,
                             "Include the specific heat capacities of rock "
                             "matrix in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteHeatConductivity,
                             "Include the lumped heat conductivity of the "
                             "medium in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteEnthalpies,
                             "Include the specific enthalpy of the phases in "
                             "the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteInternalEnergies,
                             "Include the specific internal energy of the "
                             "phases in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (enthalpyOutput_())
            this->resizePhaseBuffer_(enthalpy_);
        if (internalEnergyOutput_())
            this->resizePhaseBuffer_(internalEnergy_);

        if (solidHeatCapacityOutput_())
            this->resizeScalarBuffer_(solidHeatCapacity_);
        if (heatConductivityOutput_())
            this->resizeScalarBuffer_(heatConductivity_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto &intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto &fs = intQuants.fluidState();

            if (solidHeatCapacityOutput_())
                solidHeatCapacity_[I] = Toolbox::value(intQuants.heatCapacitySolid());
            if (heatConductivityOutput_())
                heatConductivity_[I] = Toolbox::value(intQuants.heatConductivity());

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (enthalpyOutput_())
                    enthalpy_[phaseIdx][I] = Toolbox::value(fs.enthalpy(phaseIdx));
                if (internalEnergyOutput_())
                    internalEnergy_[phaseIdx][I] = Toolbox::value(fs.internalEnergy(phaseIdx));
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

        if (solidHeatCapacityOutput_())
            this->commitScalarBuffer_(baseWriter, "heatCapacitySolid", solidHeatCapacity_);
        if (heatConductivityOutput_())
            this->commitScalarBuffer_(baseWriter, "heatConductivity", heatConductivity_);

        if (enthalpyOutput_())
            this->commitPhaseBuffer_(baseWriter, "enthalpy_%s", enthalpy_);
        if (internalEnergyOutput_())
            this->commitPhaseBuffer_(baseWriter, "internalEnergy_%s", internalEnergy_);
    }

private:
    static bool solidHeatCapacityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSolidHeatCapacity); }

    static bool heatConductivityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteHeatConductivity); }

    static bool enthalpyOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteEnthalpies); }

    static bool internalEnergyOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteInternalEnergies); }

    PhaseBuffer enthalpy_;
    PhaseBuffer internalEnergy_;

    ScalarBuffer heatConductivity_;
    ScalarBuffer solidHeatCapacity_;
};

} // namespace Ewoms

#endif
