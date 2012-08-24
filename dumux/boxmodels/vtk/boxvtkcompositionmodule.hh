// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \brief VTK output module for the fluid composition for models assuming
 *        chemical equilibrium.
 */
#ifndef DUMUX_BOX_VTK_COMPOSITION_MODULE_HH
#define DUMUX_BOX_VTK_COMPOSITION_MODULE_HH

#include "boxvtkoutputmodule.hh"

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{
namespace Properties
{
// create new type tag for the VTK composition output
NEW_TYPE_TAG(VtkComposition);

// create the property tags needed for the composition module
NEW_PROP_TAG(VtkWriteMassFractions);
NEW_PROP_TAG(VtkWriteMoleFractions);
NEW_PROP_TAG(VtkWriteTotalMassFractions);
NEW_PROP_TAG(VtkWriteTotalMoleFractions);
NEW_PROP_TAG(VtkWriteMolarities);
NEW_PROP_TAG(VtkWriteFugacities);
NEW_PROP_TAG(VtkWriteFugacityCoeffs);

// set default values for what quantities to output
SET_BOOL_PROP(VtkComposition, VtkWriteMassFractions, false);
SET_BOOL_PROP(VtkComposition, VtkWriteMoleFractions, true);
SET_BOOL_PROP(VtkComposition, VtkWriteTotalMassFractions, false);
SET_BOOL_PROP(VtkComposition, VtkWriteTotalMoleFractions, false);
SET_BOOL_PROP(VtkComposition, VtkWriteMolarities, false);
SET_BOOL_PROP(VtkComposition, VtkWriteFugacities, false);
SET_BOOL_PROP(VtkComposition, VtkWriteFugacityCoeffs, false);
}

/*!
 * \ingroup BoxModels
 *
 * \brief VTK output module for the fluid composition
 *
 * This module deals with the following quantities:
 * - Mole fraction of a component in a fluid phase
 * - Mass fraction of a component in a fluid phase
 * - Molarity (i.e. molar concentration) of a component in a fluid phase
 * - Fugacity of all components
 * - FugacityCoefficient of all components in all phases
 */
template<class TypeTag>
class BoxVtkCompositionModule : public BoxVtkOutputModule<TypeTag>
{
    typedef BoxVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };

    typedef Dumux::VtkMultiWriter<GridView> VtkMultiWriter;

    typedef typename ParentType::ComponentBuffer ComponentBuffer;
    typedef typename ParentType::PhaseComponentBuffer PhaseComponentBuffer;

public:
    BoxVtkCompositionModule(const Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (moleFracOutput_()) this->resizePhaseComponentBuffer_(moleFrac_);
        if (massFracOutput_()) this->resizePhaseComponentBuffer_(massFrac_);
        if (totalMassFracOutput_()) this->resizeComponentBuffer_(totalMassFrac_);
        if (totalMoleFracOutput_()) this->resizeComponentBuffer_(totalMoleFrac_);
        if (molarityOutput_()) this->resizePhaseComponentBuffer_(molarity_);

        if (fugacityOutput_()) this->resizeComponentBuffer_(fugacity_);
        if (fugacityCoeffOutput_()) this->resizePhaseComponentBuffer_(fugacityCoeff_);
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

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (moleFracOutput_()) moleFrac_[phaseIdx][compIdx][I] = fs.moleFraction(phaseIdx, compIdx);
                    if (massFracOutput_()) massFrac_[phaseIdx][compIdx][I] = fs.massFraction(phaseIdx, compIdx);
                    if (molarityOutput_()) molarity_[phaseIdx][compIdx][I] = fs.molarity(phaseIdx, compIdx);

                    if (fugacityCoeffOutput_()) fugacityCoeff_[phaseIdx][compIdx][I] = fs.fugacityCoefficient(phaseIdx, compIdx);
                }
            }

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (totalMassFracOutput_()) {
                    Scalar compMass = 0;
                    Scalar totalMass = 0;
                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        totalMass += fs.density(phaseIdx)*fs.saturation(phaseIdx);
                        compMass += fs.density(phaseIdx)*fs.saturation(phaseIdx)*fs.massFraction(phaseIdx, compIdx);
                    }
                    totalMassFrac_[compIdx][I] = compMass/totalMass;
                }
                if (totalMoleFracOutput_()) {
                    Scalar compMoles = 0;
                    Scalar totalMoles = 0;
                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        totalMoles += fs.molarDensity(phaseIdx)*fs.saturation(phaseIdx);
                        compMoles += fs.molarDensity(phaseIdx)*fs.saturation(phaseIdx)*fs.moleFraction(phaseIdx, compIdx);
                    }
                    totalMoleFrac_[compIdx][I] = compMoles/totalMoles;
                }
                if (fugacityOutput_()) fugacity_[compIdx][I] = volVars.fluidState().fugacity(/*phaseIdx=*/0, compIdx);
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (moleFracOutput_()) this->commitPhaseComponentBuffer_(writer, "moleFrac_%s^%s", moleFrac_);
        if (massFracOutput_()) this->commitPhaseComponentBuffer_(writer, "massFrac_%s^%s", massFrac_);
        if (molarityOutput_()) this->commitPhaseComponentBuffer_(writer, "molarity_%s^%s", molarity_);
        if (totalMassFracOutput_()) this->commitComponentBuffer_(writer, "totalMassFrac^%s", totalMassFrac_);
        if (totalMoleFracOutput_()) this->commitComponentBuffer_(writer, "totalMoleFrac^%s", totalMoleFrac_);

        if (fugacityOutput_()) this->commitComponentBuffer_(writer, "fugacity^%s", fugacity_);
        if (fugacityCoeffOutput_()) this->commitPhaseComponentBuffer_(writer, "fugacityCoeff_%s^%s", fugacityCoeff_);
    }

private:
    static bool massFracOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteMassFractions); }

    static bool moleFracOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteMoleFractions); }

    static bool totalMassFracOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteTotalMassFractions); }

    static bool totalMoleFracOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteTotalMoleFractions); }

    static bool molarityOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteMolarities); }

    static bool fugacityOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteFugacities); }

    static bool fugacityCoeffOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteFugacityCoeffs); }

    PhaseComponentBuffer moleFrac_;
    PhaseComponentBuffer massFrac_;
    PhaseComponentBuffer molarity_;
    ComponentBuffer totalMassFrac_;
    ComponentBuffer totalMoleFrac_;

    ComponentBuffer fugacity_;
    PhaseComponentBuffer fugacityCoeff_;
};

}

#endif
