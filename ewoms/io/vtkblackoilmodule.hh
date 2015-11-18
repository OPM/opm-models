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
 * \copydoc Ewoms::VtkBlackOilModule
 */
#ifndef EWOMS_VTK_BLACK_OIL_MODULE_HH
#define EWOMS_VTK_BLACK_OIL_MODULE_HH

#include <opm/material/localad/Math.hpp>

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/models/blackoil/blackoilproperties.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

namespace Ewoms {
namespace Properties {
// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(VtkBlackOil);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(VtkWriteGasDissolutionFactor);
NEW_PROP_TAG(VtkWriteOilVaporizationFactor);
NEW_PROP_TAG(VtkWriteOilFormationVolumeFactor);
NEW_PROP_TAG(VtkWriteGasFormationVolumeFactor);
NEW_PROP_TAG(VtkWriteWaterFormationVolumeFactor);
NEW_PROP_TAG(VtkWriteOilSaturationPressure);
NEW_PROP_TAG(VtkWriteGasSaturationPressure);
NEW_PROP_TAG(VtkWriteSaturatedOilGasDissolutionFactor);
NEW_PROP_TAG(VtkWriteSaturatedGasOilVaporizationFactor);
NEW_PROP_TAG(VtkWriteSaturatedOilFormationVolumeFactor);
NEW_PROP_TAG(VtkWriteSaturatedGasFormationVolumeFactor);

// set default values for what quantities to output
SET_BOOL_PROP(VtkBlackOil, VtkWriteGasDissolutionFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteOilVaporizationFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteOilFormationVolumeFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteGasFormationVolumeFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteWaterFormationVolumeFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteOilSaturationPressure, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteGasSaturationPressure, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteSaturatedOilGasDissolutionFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteSaturatedGasOilVaporizationFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteSaturatedOilFormationVolumeFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteSaturatedGasFormationVolumeFactor, false);
} // namespace Properties
} // namespace Ewoms

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the black oil model's parameters.
 */
template <class TypeTag>
class VtkBlackOilModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Ewoms::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;

public:
    VtkBlackOilModule(const Simulator &simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteGasDissolutionFactor,
                             "Include the gas dissolution factor (R_s) of the observed oil "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteOilVaporizationFactor,
                             "Include the oil vaporization factor (R_v) of the observed gas "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteOilFormationVolumeFactor,
                             "Include the oil formation volume factor (B_o) of gas saturated "
                             "oil in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteGasFormationVolumeFactor,
                             "Include the gas formation volume factor (B_g) in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteWaterFormationVolumeFactor,
                             "Include the water formation volume factor (B_w) in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteOilSaturationPressure,
                             "Include the saturation pressure of oil (p_o,sat) in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteGasSaturationPressure,
                             "Include the saturation pressure of gas (p_g,sat) in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSaturatedOilGasDissolutionFactor,
                             "Include the gas dissolution factor (R_s,sat) of gas saturated "
                             "oil in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSaturatedGasOilVaporizationFactor,
                             "Include the oil vaporization factor (R_v,sat) of oil saturated "
                             "gas in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSaturatedOilFormationVolumeFactor,
                             "Include the formation volume factor (B_o,sat) of saturated oil in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSaturatedGasFormationVolumeFactor,
                             "Include the formation volume factor (B_g,sat) of saturated gas in the "
                             "VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (gasDissolutionFactorOutput_())
            this->resizeScalarBuffer_(gasDissolutionFactor_);
        if (oilVaporizationFactorOutput_())
            this->resizeScalarBuffer_(oilVaporizationFactor_);
        if (oilFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(oilFormationVolumeFactor_);
        if (gasFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(gasFormationVolumeFactor_);
        if (waterFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(waterFormationVolumeFactor_);
        if (oilSaturationPressureOutput_())
            this->resizeScalarBuffer_(oilSaturationPressure_);
        if (gasSaturationPressureOutput_())
            this->resizeScalarBuffer_(gasSaturationPressure_);
        if (saturatedOilGasDissolutionFactorOutput_())
            this->resizeScalarBuffer_(saturatedOilGasDissolutionFactor_);
        if (saturatedGasOilVaporizationFactorOutput_())
            this->resizeScalarBuffer_(saturatedGasOilVaporizationFactor_);
        if (saturatedOilFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(saturatedOilFormationVolumeFactor_);
        if (saturatedGasFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(saturatedGasFormationVolumeFactor_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        typedef Opm::MathToolbox<Evaluation> Toolbox;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto &fs = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState();
            int globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            Scalar po = Toolbox::value(fs.pressure(oilPhaseIdx));
            Scalar pg = Toolbox::value(fs.pressure(oilPhaseIdx));
            Scalar pw = Toolbox::value(fs.pressure(oilPhaseIdx));
            Scalar T = Toolbox::value(fs.temperature(oilPhaseIdx));
            Scalar X_oG = Toolbox::value(fs.massFraction(oilPhaseIdx, gasCompIdx));
            Scalar X_gO = Toolbox::value(fs.massFraction(gasPhaseIdx, oilCompIdx));
            int regionIdx = elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex();
            Scalar rhooRef = FluidSystem::referenceDensity(oilPhaseIdx, regionIdx);
            Scalar rhogRef = FluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
            Scalar Rs = FluidSystem::convertXoGToRs(X_oG, regionIdx);
            Scalar Rv = FluidSystem::convertXgOToRv(X_gO, regionIdx);

            if (gasDissolutionFactorOutput_())
                gasDissolutionFactor_[globalDofIdx] = Rs;
            if (oilVaporizationFactorOutput_())
                oilVaporizationFactor_[globalDofIdx] = Rv;
            if (oilFormationVolumeFactorOutput_())
                oilFormationVolumeFactor_[globalDofIdx] =
                    FluidSystem::oilFormationVolumeFactor(T, po, Rs, regionIdx);
            if (gasFormationVolumeFactorOutput_())
                gasFormationVolumeFactor_[globalDofIdx] =
                    FluidSystem::gasFormationVolumeFactor(T, pg, Rv, regionIdx);
            if (waterFormationVolumeFactorOutput_())
                waterFormationVolumeFactor_[globalDofIdx] =
                    FluidSystem::waterFormationVolumeFactor(T, pw, regionIdx);
            if (oilSaturationPressureOutput_())
                oilSaturationPressure_[globalDofIdx] =
                    FluidSystem::oilSaturationPressure(T, Rs, regionIdx);
            if (gasSaturationPressureOutput_())
                gasSaturationPressure_[globalDofIdx] =
                    FluidSystem::gasSaturationPressure(T, Rv, regionIdx);
            if (saturatedOilGasDissolutionFactorOutput_())
                saturatedOilGasDissolutionFactor_[globalDofIdx] =
                    FluidSystem::template gasDissolutionFactor<Scalar>(T, po, regionIdx);
            if (saturatedGasOilVaporizationFactorOutput_())
                saturatedGasOilVaporizationFactor_[globalDofIdx] =
                    FluidSystem::template oilVaporizationFactor<Scalar>(T, pg, regionIdx);
            if (saturatedOilFormationVolumeFactorOutput_())
                saturatedOilFormationVolumeFactor_[globalDofIdx] =
                    FluidSystem::saturatedOilFormationVolumeFactor(T, po, regionIdx);
            if (saturatedGasFormationVolumeFactorOutput_())
                saturatedGasFormationVolumeFactor_[globalDofIdx] =
                    FluidSystem::saturatedOilFormationVolumeFactor(T, pg, regionIdx);
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter &baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter)
            return;

        if (gasDissolutionFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "R_s", gasDissolutionFactor_);
        if (oilVaporizationFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "R_v", oilVaporizationFactor_);
        if (oilFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "B_o", oilFormationVolumeFactor_);
        if (gasFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "B_g", gasFormationVolumeFactor_);
        if (waterFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "B_w", waterFormationVolumeFactor_);
        if (oilSaturationPressureOutput_())
            this->commitScalarBuffer_(baseWriter, "p_o,sat", oilSaturationPressure_);
        if (gasSaturationPressureOutput_())
            this->commitScalarBuffer_(baseWriter, "p_g,sat", gasSaturationPressure_);
        if (saturatedOilGasDissolutionFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "R_s,sat", saturatedOilGasDissolutionFactor_);
        if (saturatedGasOilVaporizationFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "R_v,sat", saturatedGasOilVaporizationFactor_);
        if (saturatedOilFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "B_o,sat", saturatedOilFormationVolumeFactor_);
        if (saturatedGasFormationVolumeFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "B_g,sat", saturatedGasFormationVolumeFactor_);
    }

private:
    static bool gasDissolutionFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteGasDissolutionFactor); }

    static bool oilVaporizationFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteOilVaporizationFactor); }

    static bool oilFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteOilFormationVolumeFactor); }

    static bool gasFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteGasFormationVolumeFactor); }

    static bool waterFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteWaterFormationVolumeFactor); }

    static bool oilSaturationPressureOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteOilSaturationPressure); }

    static bool gasSaturationPressureOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteGasSaturationPressure); }

    static bool saturatedOilGasDissolutionFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSaturatedOilGasDissolutionFactor); }

    static bool saturatedGasOilVaporizationFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSaturatedGasOilVaporizationFactor); }

    static bool saturatedOilFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSaturatedOilFormationVolumeFactor); }

    static bool saturatedGasFormationVolumeFactorOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSaturatedGasFormationVolumeFactor); }


    ScalarBuffer gasDissolutionFactor_;
    ScalarBuffer oilVaporizationFactor_;
    ScalarBuffer oilFormationVolumeFactor_;
    ScalarBuffer gasFormationVolumeFactor_;
    ScalarBuffer waterFormationVolumeFactor_;
    ScalarBuffer oilSaturationPressure_;
    ScalarBuffer gasSaturationPressure_;

    ScalarBuffer saturatedOilGasDissolutionFactor_;
    ScalarBuffer saturatedGasOilVaporizationFactor_;
    ScalarBuffer saturatedOilFormationVolumeFactor_;
    ScalarBuffer saturatedGasFormationVolumeFactor_;
};
} // namespace Ewoms

#endif
