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

        for (int dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto &fs = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState();
            int globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            int regionIdx = elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex();
            Scalar po = fs.pressure(oilPhaseIdx);
            Scalar XoG = fs.massFraction(oilPhaseIdx, gasCompIdx);

            if (saturationsOutput_()) {
                for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    saturation_[phaseIdx][globalDofIdx] = fs.saturation(phaseIdx);
                    Valgrind::CheckDefined(saturation_[phaseIdx][globalDofIdx]);
                }
            }
            if (pressuresOutput_()) {
                for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    pressure_[phaseIdx][globalDofIdx] = fs.pressure(phaseIdx) / 1e5;
                    Valgrind::CheckDefined(pressure_[phaseIdx][globalDofIdx]);
                }
            }
            if (gasDissolutionFactorOutput_()) {
                gasDissolutionFactor_[globalDofIdx] =
                    FluidSystem::gasDissolutionFactor(po, regionIdx);
                Valgrind::CheckDefined(gasDissolutionFactor_[globalDofIdx]);
            }
            if (gasFormationVolumeFactorOutput_()) {
                gasFormationVolumeFactor_[globalDofIdx] =
                    FluidSystem::gasFormationVolumeFactor(po, regionIdx);
                Valgrind::CheckDefined(gasFormationVolumeFactor_[globalDofIdx]);
            }
            if (saturatedOilFormationVolumeFactorOutput_()) {
                saturatedOilFormationVolumeFactor_[globalDofIdx] =
                    FluidSystem::saturatedOilFormationVolumeFactor(po, regionIdx);
                Valgrind::CheckDefined(saturatedOilFormationVolumeFactor_[globalDofIdx]);
            }
            if (oilSaturationPressureOutput_()) {
                oilSaturationPressure_[globalDofIdx] =
                    FluidSystem::oilSaturationPressure(XoG, regionIdx);
                Valgrind::CheckDefined(oilSaturationPressure_[globalDofIdx]);
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
        if (pressuresOutput_()) {
            this->commitScalarBuffer_(writer, "PRESSURE", pressure_[oilPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "PGAS", pressure_[gasPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "PWAT", pressure_[waterPhaseIdx], bufferType);
        }
        if (saturationsOutput_()) {
            this->commitScalarBuffer_(writer, "SWAT", saturation_[waterPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "SGAS", saturation_[gasPhaseIdx], bufferType);
            // the oil saturation is _NOT_ written to disk. Instead, it is calculated by
            // the visualization tool. Wondering why is probably a waste of time...
        }
        if (gasDissolutionFactorOutput_())
            this->commitScalarBuffer_(writer, "RS", gasDissolutionFactor_, bufferType);
        if (gasFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(writer, "BG", gasFormationVolumeFactor_, bufferType);
        if (saturatedOilFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(writer, "BOSAT", saturatedOilFormationVolumeFactor_, bufferType);
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
