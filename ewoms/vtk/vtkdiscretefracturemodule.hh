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
 * \copydoc Ewoms::VtkDiscreteFractureModule
 */
#ifndef EWOMS_VTK_DISCRETE_FRACTURE_MODULE_HH
#define EWOMS_VTK_DISCRETE_FRACTURE_MODULE_HH

#include <ewoms/vtk/vtkbaseoutputmodule.hh>

#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

namespace Opm {
namespace Properties {
// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(VtkDiscreteFracture);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(VtkWriteFractureSaturations);
NEW_PROP_TAG(VtkWriteFractureMobilities);
NEW_PROP_TAG(VtkWriteFractureRelativePermeabilities);
NEW_PROP_TAG(VtkWriteFracturePorosity);
NEW_PROP_TAG(VtkWriteFractureIntrinsicPermeabilities);
NEW_PROP_TAG(VtkWriteFractureFilterVelocities);
NEW_PROP_TAG(VtkWriteFractureVolumeFraction);
NEW_PROP_TAG(DiscVtkBaseOutputModule);

// set default values for what quantities to output
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureSaturations, true);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureMobilities, false);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureRelativePermeabilities, true);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFracturePorosity, true);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureIntrinsicPermeabilities,
              false);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureFilterVelocities, false);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureVolumeFraction, true);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for quantities which make sense for all
 *        models which deal with discrete fractures in porous media.
 *
 * This module deals with the following quantities:
 * - Saturations of all fluid phases in the fracture
 * - Mobilities of all fluid phases in the fracture
 * - Relative permeabilities of all fluid phases in the fracture
 * - Porosity of the medium in the fracture
 * - Norm of the intrinsic permeability of the medium in the fracture
 */
template <class TypeTag>
class VtkDiscreteFractureModule : public VtkBaseOutputModule<TypeTag>
{
    typedef VtkBaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, DiscVtkBaseOutputModule) DiscVtkBaseOutputModule;

    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;

    typedef Dune::FieldVector<Scalar, dim> VelocityVector;
    typedef Dune::BlockVector<VelocityVector> VelocityField;
    typedef std::array<VelocityField, numPhases> PhaseVectorField;

public:
    VtkDiscreteFractureModule(const Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureSaturations,
                             "Include the phase saturations in the VTK output "
                             "files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureMobilities,
                             "Include the phase mobilities in the VTK output "
                             "files");
        EWOMS_REGISTER_PARAM(TypeTag, bool,
                             VtkWriteFractureRelativePermeabilities,
                             "Include the phase relative permeabilities in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFracturePorosity,
                             "Include the porosity in the VTK output files");
        EWOMS_REGISTER_PARAM(
            TypeTag, bool, VtkWriteFractureIntrinsicPermeabilities,
            "Include the intrinsic permeability in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureFilterVelocities,
                             "Include in the filter velocities of the phases "
                             "the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureVolumeFraction,
                             "Add the fraction of the total volume which is "
                             "occupied by fractures in the VTK output");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (saturationOutput_())
            this->resizePhaseBuffer_(saturation_);
        if (mobilityOutput_())
            this->resizePhaseBuffer_(mobility_);
        if (relativePermeabilityOutput_())
            this->resizePhaseBuffer_(relativePermeability_);

        if (porosityOutput_())
            this->resizeScalarBuffer_(porosity_);
        if (intrinsicPermeabilityOutput_())
            this->resizeScalarBuffer_(intrinsicPermeability_);
        if (volumeFractionOutput_())
            this->resizeScalarBuffer_(volumeFraction_);

        if (velocityOutput_()) {
            Scalar nElem = this->problem_.gridView().size(/*codim=*/0);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                velocity_[phaseIdx].resize(nElem);
                velocity_[phaseIdx] = 0;
            }
            this->resizePhaseBuffer_(velocityWeight_);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        const auto &fractureMapper = elemCtx.problem().fractureMapper();

        for (int i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            if (!fractureMapper.isFractureVertex(I))
                continue;

            const auto &volVars = elemCtx.volVars(i, /*timeIdx=*/0);
            const auto &fs = volVars.fractureFluidState();

            if (porosityOutput_())
                porosity_[I] = volVars.fracturePorosity();
            if (intrinsicPermeabilityOutput_()) {
                const auto &K = volVars.fractureIntrinsicPermeability();
                intrinsicPermeability_[I] = K[0][0];
            }

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (saturationOutput_())
                    saturation_[phaseIdx][I] = fs.saturation(phaseIdx);
                if (mobilityOutput_())
                    mobility_[phaseIdx][I] = volVars.fractureMobility(phaseIdx);
                if (relativePermeabilityOutput_())
                    relativePermeability_[phaseIdx][I]
                        = volVars.fractureRelativePermeability(phaseIdx);
                if (volumeFractionOutput_())
                    volumeFraction_[I] += volVars.fractureVolume();
            }
        }

        if (velocityOutput_()) {
            // calculate velocities if requested by the problem
            for (int scvfIdx = 0; scvfIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++ scvfIdx) {
                const auto &fluxVars = elemCtx.fluxVars(scvfIdx, /*timeIdx=*/0);

                int i = fluxVars.interiorIndex();
                int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                int j = fluxVars.exteriorIndex();
                int J = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                if (!fractureMapper.isFractureEdge(I, J))
                    continue;

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight
                        = std::max(1e-16, std::abs(fluxVars.fractureVolumeFlux(
                                              phaseIdx)));
                    Valgrind::CheckDefined(fluxVars.extrusionFactor());
                    assert(fluxVars.extrusionFactor() > 0);
                    weight *= fluxVars.extrusionFactor();

                    Dune::FieldVector<Scalar, dim> v(
                        fluxVars.fractureFilterVelocity(phaseIdx));
                    v *= weight;

                    velocity_[phaseIdx][I] += v;
                    velocity_[phaseIdx][J] += v;

                    velocityWeight_[phaseIdx][I] += weight;
                    velocityWeight_[phaseIdx][J] += weight;
                } // end for all phases
            }     // end for all faces
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (saturationOutput_())
            this->commitPhaseBuffer_(writer, "fractureSaturation_%s",
                                     saturation_);
        if (mobilityOutput_())
            this->commitPhaseBuffer_(writer, "fractureMobility_%s", mobility_);
        if (relativePermeabilityOutput_())
            this->commitPhaseBuffer_(writer, "fractureRelativePerm_%s",
                                     relativePermeability_);

        if (porosityOutput_())
            this->commitScalarBuffer_(writer, "fracturePorosity", porosity_);
        if (intrinsicPermeabilityOutput_())
            this->commitScalarBuffer_(writer, "fractureIntrinsicPerm",
                                      intrinsicPermeability_);
        if (volumeFractionOutput_()) {
            // divide the fracture volume by the total volume of the finite
            // volumes
            for (unsigned I = 0; I < volumeFraction_.size(); ++I)
                volumeFraction_[I] /= this->problem_.model().dofTotalVolume(I);
            this->commitScalarBuffer_(writer, "fractureVolumeFraction",
                                      volumeFraction_);
        }

        if (velocityOutput_()) {
            int nDof = this->problem_.model().numDof();

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (int i = 0; i < nDof; ++i)
                    velocity_[phaseIdx][i]
                        /= std::max<Scalar>(1e-20, velocityWeight_[phaseIdx][i]);
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "fractureFilterVelocity_%s", FluidSystem::phaseName(phaseIdx));

                DiscVtkBaseOutputModule::attachDofData_(writer, velocity_[phaseIdx], name, dim);
            }
        }
    }

private:
    static bool saturationOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureSaturations); }

    static bool mobilityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureMobilities); }

    static bool relativePermeabilityOutput_()
    {
        return EWOMS_GET_PARAM(TypeTag, bool,
                               VtkWriteFractureRelativePermeabilities);
    }

    static bool porosityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFracturePorosity); }

    static bool intrinsicPermeabilityOutput_()
    {
        return EWOMS_GET_PARAM(TypeTag, bool,
                               VtkWriteFractureIntrinsicPermeabilities);
    }

    static bool volumeFractionOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureVolumeFraction); }

    static bool velocityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureFilterVelocities); }

    PhaseBuffer pressure_;
    PhaseBuffer density_;
    PhaseBuffer saturation_;
    PhaseBuffer mobility_;
    PhaseBuffer relativePermeability_;
    PhaseBuffer viscosity_;
    PhaseBuffer averageMolarMass_;

    ScalarBuffer porosity_;
    ScalarBuffer volumeFraction_;
    ScalarBuffer intrinsicPermeability_;

    PhaseVectorField velocity_;
    PhaseBuffer velocityWeight_;

    PhaseVectorField potentialGradient_;
    PhaseBuffer potentialWeight_;
};

} // namespace Ewoms

#endif
