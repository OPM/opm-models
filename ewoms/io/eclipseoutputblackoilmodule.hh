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
 * \copydoc Ewoms::EclipseOutputBlackOilModule
 */
#ifndef EWOMS_ECLIPSE_OUTPUT_BLACK_OIL_MODULE_HH
#define EWOMS_ECLIPSE_OUTPUT_BLACK_OIL_MODULE_HH

#include <ewoms/io/baseoutputmodule.hh>
#include <ewoms/io/eclipsewriter.hh>

#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/common/fvector.hh>

#include <type_traits>

namespace Opm {
namespace Properties {
// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(EclipseOutputBlackOil);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(EclipseOutputWriteSaturations);
NEW_PROP_TAG(EclipseOutputWritePressures);
NEW_PROP_TAG(EclipseOutputWriteGasDissolutionFactor);
NEW_PROP_TAG(EclipseOutputWriteGasFormationVolumeFactor);
NEW_PROP_TAG(EclipseOutputWriteOilFormationVolumeFactor);
NEW_PROP_TAG(EclipseOutputWriteOilSaturationPressure);

// set default values for what quantities to output
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteSaturations, true);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWritePressures, true);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteGasDissolutionFactor, true);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteGasFormationVolumeFactor, true);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteOilFormationVolumeFactor, true);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteOilSaturationPressure, true);
}} // namespace Opm, Properties

namespace Ewoms {

// forward declaration
template <class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup EclipseOutput
 *
 * \brief Output module for the results black oil model writing in
 *        Eclipse binary format.
 */
template <class TypeTag>
class EclipseOutputBlackOilModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Ewoms::EclipseWriter<TypeTag> EclipseWriter;

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;

public:
    EclipseOutputBlackOilModule(const Simulator &simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteSaturations,
                             "Include the saturations of all fluid phases in the "
                             "Eclipse output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWritePressures,
                             "Include the absolute pressures of all fluid phases in the "
                             "Eclipse output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteGasDissolutionFactor,
                             "Include the gas dissolution factor in the "
                             "Eclipse output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteGasFormationVolumeFactor,
                             "Include the gas formation volume factor in the "
                             "Eclipse output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteOilFormationVolumeFactor,
                             "Include the oil formation volume factor of saturated oil "
                             "in the Eclipse output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteOilSaturationPressure,
                             "Include the saturation pressure of oil in the "
                             "Eclipse output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        auto bufferType = ParentType::ElementBuffer;
        if (saturationsOutput_()) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                this->resizeScalarBuffer_(saturation_[phaseIdx], bufferType);
        }
        if (pressuresOutput_()) {
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                this->resizeScalarBuffer_(pressure_[phaseIdx], bufferType);
        }
        if (gasDissolutionFactorOutput_())
            this->resizeScalarBuffer_(gasDissolutionFactor_, bufferType);
        if (gasFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(gasFormationVolumeFactor_, bufferType);
        if (saturatedOilFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(saturatedOilFormationVolumeFactor_, bufferType);
        if (oilSaturationPressureOutput_())
            this->resizeScalarBuffer_(oilSaturationPressure_, bufferType);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        for (int i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            const auto &fs = elemCtx.intensiveQuantities(/*spaceIdx=*/i, /*timeIdx=*/0).fluidState();
            int I = elemCtx.globalSpaceIndex(/*spaceIdx=*/i, /*timeIdx=*/0);
            Scalar po = fs.pressure(oilPhaseIdx);
            Scalar XoG = fs.massFraction(oilPhaseIdx, gasCompIdx);

            if (saturationsOutput_()) {
                for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    saturation_[phaseIdx][I] = fs.saturation(phaseIdx);
                    Valgrind::CheckDefined(saturation_[phaseIdx][I]);
                }
            }
            if (pressuresOutput_()) {
                for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    pressure_[phaseIdx][I] = fs.pressure(phaseIdx);
                    Valgrind::CheckDefined(pressure_[phaseIdx][I]);
                }
            }
            if (gasDissolutionFactorOutput_()) {
                gasDissolutionFactor_[I] = FluidSystem::gasDissolutionFactor(po);
                Valgrind::CheckDefined(gasDissolutionFactor_[I]);
            }
            if (gasFormationVolumeFactorOutput_()) {
                gasFormationVolumeFactor_[I] = FluidSystem::gasFormationVolumeFactor(po);
                Valgrind::CheckDefined(gasFormationVolumeFactor_[I]);
            }
            if (saturatedOilFormationVolumeFactorOutput_()) {
                saturatedOilFormationVolumeFactor_[I] = FluidSystem::saturatedOilFormationVolumeFactor(po);
                Valgrind::CheckDefined(saturatedOilFormationVolumeFactor_[I]);
            }
            if (oilSaturationPressureOutput_()) {
                oilSaturationPressure_[I] = FluidSystem::oilSaturationPressure(XoG);
                Valgrind::CheckDefined(oilSaturationPressure_[I]);
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter &writer)
    {
        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        if (!dynamic_cast<EclipseWriter*>(&writer))
            return; // this module only consideres eclipse writers...

        typename ParentType::BufferType bufferType = ParentType::ElementBuffer;
        if (saturationsOutput_()) {
            this->commitScalarBuffer_(writer, "SOIL", saturation_[oilPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "SGAS", saturation_[gasPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "SWAT", saturation_[waterPhaseIdx], bufferType);
        }
        if (pressuresOutput_()) {
            this->commitScalarBuffer_(writer, "PRESSURE", pressure_[oilPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "PGAS", pressure_[gasPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "PWAT", pressure_[waterPhaseIdx], bufferType);
        }
        if (gasDissolutionFactorOutput_())
            this->commitScalarBuffer_(writer, "RS", gasDissolutionFactor_, bufferType);
        if (gasFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(writer, "BG", gasFormationVolumeFactor_, bufferType);
        if (saturatedOilFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(writer, "BP", saturatedOilFormationVolumeFactor_, bufferType);
        if (oilSaturationPressureOutput_())
            this->commitScalarBuffer_(writer, "PSAT", oilSaturationPressure_);
    }

private:
    static bool saturationsOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteSaturations); }

    static bool pressuresOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWritePressures); }

    static bool gasDissolutionFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteGasDissolutionFactor); }

    static bool gasFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteGasFormationVolumeFactor); }

    static bool saturatedOilFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteOilFormationVolumeFactor); }

    static bool oilSaturationPressureOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteOilSaturationPressure); }

    ScalarBuffer saturation_[numPhases];
    ScalarBuffer pressure_[numPhases];
    ScalarBuffer gasDissolutionFactor_;
    ScalarBuffer gasFormationVolumeFactor_;
    ScalarBuffer saturatedOilFormationVolumeFactor_;
    ScalarBuffer oilSaturationPressure_;
};
} // namespace Ewoms

#endif
