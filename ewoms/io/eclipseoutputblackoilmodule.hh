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

#include <cstdio>

namespace Opm {
namespace Properties {
// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(EclipseOutputBlackOil);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(EclipseOutputWriteGasFormationFactor);
NEW_PROP_TAG(EclipseOutputWriteGasFormationVolumeFactor);
NEW_PROP_TAG(EclipseOutputWriteOilFormationVolumeFactor);
NEW_PROP_TAG(EclipseOutputWriteOilSaturationPressure);

// set default values for what quantities to output
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteGasFormationFactor, false);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteGasFormationVolumeFactor, false);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteOilFormationVolumeFactor, false);
SET_BOOL_PROP(EclipseOutputBlackOil, EclipseOutputWriteOilSaturationPressure, false);
}} // namespace Opm, Properties

namespace Ewoms {
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
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Ewoms::EclipseWriter<TypeTag> EclipseWriter;

    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
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
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteGasFormationFactor,
                             "Include the gas formation factor in the VTK "
                             "output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteGasFormationVolumeFactor,
                             "Include the gas formation volume factor in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteOilFormationVolumeFactor,
                             "Include the oil formation volume factor in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclipseOutputWriteOilSaturationPressure,
                             "Include the saturation pressure of oil in the "
                             "VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
#warning TODO
#if 0
        auto bufferType = ParentType::ElementBuffer;
        if (gasDissolutionFactorOutput_())
            this->resizeScalarBuffer_(gasDissolutionFactor_, bufferType);
        if (gasFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(gasFormationVolumeFactor_, bufferType);
        if (oilFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(oilFormationVolumeFactor_, bufferType);
        if (oilSaturationPressureOutput_())
            this->resizeScalarBuffer_(oilSaturationPressure_, bufferType);
#endif
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
#warning TODO
#if 0
        for (int i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            const auto &fs
                = elemCtx.volVars(/*spaceIdx=*/i, /*timeIdx=*/0).fluidState();
            int I = elemCtx.globalSpaceIndex(/*spaceIdx=*/i, /*timeIdx=*/0);
            Scalar po = fs.pressure(oilPhaseIdx);
            Scalar X_oG = fs.massFraction(oilPhaseIdx, gasCompIdx);

            if (gasDissolutionFactorOutput_())
                gasDissolutionFactor_[I] = FluidSystem::gasDissolutionFactor(po);
            if (gasFormationVolumeFactorOutput_())
                gasFormationVolumeFactor_[I]
                    = FluidSystem::gasFormationVolumeFactor(po);
            if (oilFormationVolumeFactorOutput_())
                oilFormationVolumeFactor_[I]
                    = FluidSystem::oilFormationVolumeFactor(po);
            if (oilSaturationPressureOutput_())
                oilSaturationPressure_[I]
                    = FluidSystem::oilSaturationPressure(X_oG);
        }
#endif
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter &writer)
    {
        if (!dynamic_cast<EclipseWriter*>(&writer))
            return; // this module only consideres eclipse writers...

#warning TODO
/*
        if (gasDissolutionFactorOutput_())
            this->commitScalarBuffer_(writer, "RS", gasDissolutionFactor_, bufferType);
        if (gasFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(writer, "BG", gasFormationVolumeFactor_, bufferType);
        if (oilFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(writer, "BP", oilFormationVolumeFactor_, bufferType);
        if (oilSaturationPressureOutput_())
            this->commitScalarBuffer_(writer, "PRESSURE", oilSaturationPressure_);
*/
    }

private:
#warning TODO
/*
    static bool gasDissolutionFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteGasFormationFactor); }

    static bool gasFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteGasFormationVolumeFactor); }

    static bool oilFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteOilFormationVolumeFactor); }

    static bool oilSaturationPressureOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclipseOutputWriteOilSaturationPressure); }
*/
    ScalarBuffer gasDissolutionFactor_;
    ScalarBuffer gasFormationVolumeFactor_;
    ScalarBuffer oilFormationVolumeFactor_;
    ScalarBuffer oilSaturationPressure_;
};
} // namespace Ewoms

#endif
