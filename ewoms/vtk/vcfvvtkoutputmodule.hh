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
 * \copydoc Ewoms::VcfvVtkOutputModule
 */
#ifndef EWOMS_VCFV_VTK_OUTPUT_MODULE_HH
#define EWOMS_VCFV_VTK_OUTPUT_MODULE_HH

#include <ewoms/io/vtkmultiwriter.hh>
#include <ewoms/common/parametersystem.hh>
#include <opm/core/utility/PropertySystem.hpp>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

#include <vector>
#include <sstream>
#include <string>
#include <array>

#include <cstdio>

namespace Opm {
namespace Properties {
// forward definition of property tags
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
NEW_PROP_TAG(NumEq);

NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ElementContext);
NEW_PROP_TAG(VtkMultiWriter);
NEW_PROP_TAG(FluidSystem);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \brief A the base class for VTK writer modules.
 *
 * This class also provides some convenience methods for buffer
 * management and is the base class for all other VTK writer modules.
 */
template<class TypeTag>
class VcfvVtkOutputModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { dim = GridView::dimension };

public:
    typedef std::vector<Dune::FieldVector<Scalar, 1> > ScalarBuffer;
    typedef std::array<ScalarBuffer, numEq> EqBuffer;
    typedef std::array<ScalarBuffer, numPhases> PhaseBuffer;
    typedef std::array<ScalarBuffer, numComponents> ComponentBuffer;
    typedef std::array<ComponentBuffer, numPhases> PhaseComponentBuffer;

    VcfvVtkOutputModule(const Problem &problem)
        : problem_(problem)
    {
    }

    virtual ~VcfvVtkOutputModule()
    {}

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    virtual void allocBuffers(VtkMultiWriter &writer) = 0;

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    virtual void processElement(const ElementContext &elemCtx) = 0;

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    virtual void commitBuffers(VtkMultiWriter &writer) = 0;

protected:
    /*!
     * \brief Allocate the space for a buffer storing a scalar quantity
     */
    void resizeScalarBuffer_(ScalarBuffer &buffer,
                             bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        buffer.resize(n);
        std::fill(buffer.begin(), buffer.end(), 0.0);
    }

    /*!
     * \brief Allocate the space for a buffer storing a equation specific
     *        quantity
     */
    void resizeEqBuffer_(EqBuffer &buffer,
                         bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int i = 0; i < numEq; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a phase-specific
     *        quantity
     */
    void resizePhaseBuffer_(PhaseBuffer &buffer,
                            bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int i = 0; i < numPhases; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a component
     *        specific quantity
     */
    void resizeComponentBuffer_(ComponentBuffer &buffer,
                                bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int i = 0; i < numComponents; ++i) {
            buffer[i].resize(n);
            std::fill(buffer[i].begin(), buffer[i].end(), 0.0);
        }
    }

    /*!
     * \brief Allocate the space for a buffer storing a phase and
     *        component specific buffer
     */
    void resizePhaseComponentBuffer_(PhaseComponentBuffer &buffer,
                                     bool vertexCentered = true)
    {
        Scalar n;
        if (vertexCentered)
            n = problem_.gridView().size(dim);
        else
            n = problem_.gridView().size(0);

        for (int i = 0; i < numPhases; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                buffer[i][j].resize(n);
                std::fill(buffer[i][j].begin(), buffer[i][j].end(), 0.0);
            }
        }
    }

    /*!
     * \brief Add a phase-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitScalarBuffer_(MultiWriter &writer,
                             const char *name,
                             ScalarBuffer &buffer,
                             bool vertexCentered = true)
    {
        if (vertexCentered)
            writer.attachVertexData(buffer, name, 1);
        else
            writer.attachCellData(buffer, name, 1);
    }

    /*!
     * \brief Add a buffer with as many variables as PDEs to the VTK result file.
     */
    template <class MultiWriter>
    void commitPriVarsBuffer_(MultiWriter &writer,
                              const char *pattern,
                              EqBuffer &buffer,
                              bool vertexCentered = true)
    {
        char name[512];
        for (int i = 0; i < numEq; ++i) {
            std::string eqName = problem_.model().primaryVarName(i);
            snprintf(name, 512, pattern, eqName.c_str());

            if (vertexCentered)
                writer.attachVertexData(buffer[i], name, 1);
            else
                writer.attachCellData(buffer[i], name, 1);
        }
    }

    /*!
     * \brief Add a buffer with as many variables as PDEs to the VTK result file.
     */
    template <class MultiWriter>
    void commitEqBuffer_(MultiWriter &writer,
                         const char *pattern,
                         EqBuffer &buffer,
                         bool vertexCentered = true)
    {
        char name[512];
        for (int i = 0; i < numEq; ++i) {
            std::ostringstream oss;
            oss << i;
            snprintf(name, 512, pattern, oss.str().c_str());

            if (vertexCentered)
                writer.attachVertexData(buffer[i], name, 1);
            else
                writer.attachCellData(buffer[i], name, 1);
        }
    }

    /*!
     * \brief Add a phase-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitPhaseBuffer_(MultiWriter &writer,
                            const char *pattern,
                            PhaseBuffer &buffer,
                            bool vertexCentered = true)
    {
        char name[512];
        for (int i = 0; i < numPhases; ++i) {
            snprintf(name, 512, pattern, FluidSystem::phaseName(i));

            if (vertexCentered)
                writer.attachVertexData(buffer[i], name, 1);
            else
                writer.attachCellData(buffer[i], name, 1);
        }
    }

    /*!
     * \brief Add a component-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitComponentBuffer_(MultiWriter &writer,
                                const char *pattern,
                                ComponentBuffer &buffer,
                                bool vertexCentered = true)
    {
        char name[512];
        for (int i = 0; i < numComponents; ++i) {
            snprintf(name, 512, pattern, FluidSystem::componentName(i));

            if (vertexCentered)
                writer.attachVertexData(buffer[i], name, 1);
            else
                writer.attachCellData(buffer[i], name, 1);
        }
    }

    /*!
     * \brief Add a phase and component specific quantities to the output.
     */
    template <class MultiWriter>
    void commitPhaseComponentBuffer_(MultiWriter &writer,
                                     const char *pattern,
                                     PhaseComponentBuffer &buffer,
                                     bool vertexCentered = true)
    {
        char name[512];
        for (int i= 0; i < numPhases; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                snprintf(name, 512, pattern,
                         FluidSystem::phaseName(i),
                         FluidSystem::componentName(j));

                if (vertexCentered)
                    writer.attachVertexData(buffer[i][j], name, 1);
                else
                    writer.attachCellData(buffer[i][j], name, 1);
            }
        }
    }

    const Problem &problem_;
};

} // namespace Ewoms

#endif
