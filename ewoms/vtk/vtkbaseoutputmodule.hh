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
 * \copydoc Ewoms::VtkBaseOutputModule
 */
#ifndef EWOMS_VTK_BASE_OUTPUT_MODULE_HH
#define EWOMS_VTK_BASE_OUTPUT_MODULE_HH

#include <ewoms/io/vtkmultiwriter.hh>
#include <ewoms/common/parametersystem.hh>

#include <opm/core/utility/PropertySystem.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

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

NEW_PROP_TAG(Model);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ElementContext);
NEW_PROP_TAG(VtkMultiWriter);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(DiscVtkBaseOutputModule);
}} // namespace Properties, Opm

namespace Ewoms {
/*!
 * \brief A the base class for VTK writer modules.
 *
 * This class also provides some convenience methods for buffer
 * management and is the base class for all other VTK writer modules.
 */
template<class TypeTag>
class VtkBaseOutputModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, DiscVtkBaseOutputModule) DiscVtkBaseOutputModule;

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

    VtkBaseOutputModule(const Problem &problem)
        : problem_(problem)
    {
    }

    virtual ~VtkBaseOutputModule()
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
    enum BufferType {
        DofBuffer, //!< Buffer contains data associated with the degrees of freedom
        VertexBuffer, //!< Buffer contains data associated with the grid's vertices
        ElementBuffer //!< Buffer contains data associated with the grid's elements
    };

    /*!
     * \brief Allocate the space for a buffer storing a scalar quantity
     */
    void resizeScalarBuffer_(ScalarBuffer &buffer,
                             BufferType bufferType = DofBuffer)
    {
        Scalar n;
        if (bufferType == VertexBuffer)
            n = problem_.gridView().size(dim);
        else if (bufferType == ElementBuffer)
            n = problem_.gridView().size(0);
        else if (bufferType == DofBuffer)
            n = problem_.model().numDof();
        else
            OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");

        buffer.resize(n);
        std::fill(buffer.begin(), buffer.end(), 0.0);
    }

    /*!
     * \brief Allocate the space for a buffer storing a equation specific
     *        quantity
     */
    void resizeEqBuffer_(EqBuffer &buffer,
                         BufferType bufferType = DofBuffer)
    {
        Scalar n;
        if (bufferType == VertexBuffer)
            n = problem_.gridView().size(dim);
        else if (bufferType == ElementBuffer)
            n = problem_.gridView().size(0);
        else if (bufferType == DofBuffer)
            n = problem_.model().numDof();
        else
            OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");

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
                            BufferType bufferType = DofBuffer)
    {
        Scalar n;
        if (bufferType == VertexBuffer)
            n = problem_.gridView().size(dim);
        else if (bufferType == ElementBuffer)
            n = problem_.gridView().size(0);
        else if (bufferType == DofBuffer)
            n = problem_.model().numDof();
        else
            OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");

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
                                BufferType bufferType = DofBuffer)
    {
        Scalar n;
        if (bufferType == VertexBuffer)
            n = problem_.gridView().size(dim);
        else if (bufferType == ElementBuffer)
            n = problem_.gridView().size(0);
        else if (bufferType == DofBuffer)
            n = problem_.model().numDof();
        else
            OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");

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
                                     BufferType bufferType = DofBuffer)
    {
        Scalar n;
        if (bufferType == VertexBuffer)
            n = problem_.gridView().size(dim);
        else if (bufferType == ElementBuffer)
            n = problem_.gridView().size(0);
        else if (bufferType == DofBuffer)
            n = problem_.model().numDof();
        else
            OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");

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
                             BufferType bufferType = DofBuffer)
    {
        if (bufferType == DofBuffer)
            DiscVtkBaseOutputModule::attachDofData_(writer, buffer, name, 1);
        else if (bufferType == VertexBuffer)
            attachVertexData_(writer, buffer, name, 1);
        else if (bufferType == ElementBuffer)
            attachElementData_(writer, buffer, name, 1);
        else
            OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");
    }

    /*!
     * \brief Add a buffer with as many variables as PDEs to the VTK result file.
     */
    template <class MultiWriter>
    void commitPriVarsBuffer_(MultiWriter &writer,
                              const char *pattern,
                              EqBuffer &buffer,
                              BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (int i = 0; i < numEq; ++i) {
            std::string eqName = problem_.model().primaryVarName(i);
            snprintf(name, 512, pattern, eqName.c_str());

            if (bufferType == DofBuffer)
                DiscVtkBaseOutputModule::attachDofData_(writer, buffer[i], name, 1);
            else if (bufferType == VertexBuffer)
                attachVertexData_(writer, buffer[i], name, 1);
            else if (bufferType == ElementBuffer)
                attachElementData_(writer, buffer[i], name, 1);
            else
                OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a buffer with as many variables as PDEs to the VTK result file.
     */
    template <class MultiWriter>
    void commitEqBuffer_(MultiWriter &writer,
                         const char *pattern,
                         EqBuffer &buffer,
                         BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (int i = 0; i < numEq; ++i) {
            std::ostringstream oss;
            oss << i;
            snprintf(name, 512, pattern, oss.str().c_str());

            if (bufferType == DofBuffer)
                DiscVtkBaseOutputModule::attachDofData_(writer, buffer[i], name, 1);
            else if (bufferType == VertexBuffer)
                attachVertexData_(writer, buffer[i], name, 1);
            else if (bufferType == ElementBuffer)
                attachElementData_(writer, buffer[i], name, 1);
            else
                OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a phase-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitPhaseBuffer_(MultiWriter &writer,
                            const char *pattern,
                            PhaseBuffer &buffer,
                            BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (int i = 0; i < numPhases; ++i) {
            snprintf(name, 512, pattern, FluidSystem::phaseName(i));

            if (bufferType == DofBuffer)
                DiscVtkBaseOutputModule::attachDofData_(writer, buffer[i], name, 1);
            else if (bufferType == VertexBuffer)
                attachVertexData_(writer, buffer[i], name, 1);
            else if (bufferType == ElementBuffer)
                attachElementData_(writer, buffer[i], name, 1);
            else
                OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a component-specific buffer to the VTK result file.
     */
    template <class MultiWriter>
    void commitComponentBuffer_(MultiWriter &writer,
                                const char *pattern,
                                ComponentBuffer &buffer,
                                BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (int i = 0; i < numComponents; ++i) {
            snprintf(name, 512, pattern, FluidSystem::componentName(i));

            if (bufferType == DofBuffer)
                DiscVtkBaseOutputModule::attachDofData_(writer, buffer[i], name, 1);
            else if (bufferType == VertexBuffer)
                attachVertexData_(writer, buffer[i], name, 1);
            else if (bufferType == ElementBuffer)
                attachElementData_(writer, buffer[i], name, 1);
            else
                OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");
        }
    }

    /*!
     * \brief Add a phase and component specific quantities to the output.
     */
    template <class MultiWriter>
    void commitPhaseComponentBuffer_(MultiWriter &writer,
                                     const char *pattern,
                                     PhaseComponentBuffer &buffer,
                                     BufferType bufferType = DofBuffer)
    {
        char name[512];
        for (int i= 0; i < numPhases; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                snprintf(name, 512, pattern,
                         FluidSystem::phaseName(i),
                         FluidSystem::componentName(j));

                if (bufferType == DofBuffer)
                    DiscVtkBaseOutputModule::attachDofData_(writer, buffer[i][j], name, 1);
                else if (bufferType == VertexBuffer)
                    attachVertexData_(writer, buffer[i][j], name, 1);
                else if (bufferType == ElementBuffer)
                    attachElementData_(writer, buffer[i][j], name, 1);
                else
                    OPM_THROW(std::logic_error, "bufferType must be one of Dof, Vertex or Element");
            }
        }
    }

    template <class MultiWriter, class Buffer>
    void attachElementData_(MultiWriter &writer, Buffer &buffer, const char *name, int numComponents)
    { writer.attachElementData(buffer, name, numComponents); }

    template <class MultiWriter, class Buffer>
    void attachVertexData_(MultiWriter &writer, Buffer &buffer, const char *name, int numComponents)
    { writer.attachVertexData(buffer, name, numComponents); }

    const Problem &problem_;
};

} // namespace Ewoms

#endif
