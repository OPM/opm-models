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
 * \brief VTK output module for the primary variables and boundary conditions
 */
#ifndef DUMUX_BOX_VTK_PRIMARY_VARS_MODULE_HH
#define DUMUX_BOX_VTK_PRIMARY_VARS_MODULE_HH

#include "boxvtkoutputmodule.hh"

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
// create new type tag for the VTK primary variables output
NEW_TYPE_TAG(VtkPrimaryVars);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePrimaryVars);
NEW_PROP_TAG(VtkWriteBoundaryTypes);
NEW_PROP_TAG(VtkWriteNeumann);
NEW_PROP_TAG(VtkWriteDirichlet);

SET_BOOL_PROP(VtkPrimaryVars, VtkWritePrimaryVars, false);
SET_BOOL_PROP(VtkPrimaryVars, VtkWriteBoundaryTypes, false);
SET_BOOL_PROP(VtkPrimaryVars, VtkWriteNeumann, false);
SET_BOOL_PROP(VtkPrimaryVars, VtkWriteDirichlet, false);
};

/*!
 * \ingroup BoxModels
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class BoxVtkPrimaryVarsModule : public BoxVtkOutputModule<TypeTag>
{
    typedef BoxVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { dim = GridView::dimension };

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef Dumux::VtkMultiWriter<GridView> VtkMultiWriter;

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::EqBuffer EqBuffer;

public:
    BoxVtkPrimaryVarsModule(const Problem &problem)
        : ParentType(problem)
    {
    }


    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (boundaryTypesOutput_()) this->resizeScalarBuffer_(boundaryTypes_);
        if (neumannOutput_()) this->resizeEqBuffer_(neumann_);
        if (dirichletOutput_()) this->resizeEqBuffer_(dirichlet_);
        if (primaryVarsOutput_()) this->resizeEqBuffer_(primaryVars_);
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
            const auto &priVars = elemCtx.primaryVars(i, /*timeIdx=*/0);

            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (primaryVarsOutput_()) primaryVars_[eqIdx][I] = priVars[eqIdx];

                if (boundaryTypesOutput_()) {
                    // calculate a single value for the boundary type: use one
                    // bit for each equation and set it to 1 if the equation
                    // is used for a dirichlet condition
                    int tmp = 0;
                    for (int j = 0; j < numEq; ++j) {
                        if (elemCtx.boundaryTypes(i, /*timeIdx=*/0).isDirichlet(j))
                            tmp += (1 << j);
                    }
                    boundaryTypes_[I] = tmp;
                }

                if (neumannOutput_()) {
                    DUNE_THROW(Dune::NotImplemented,
                               "VTK output Neumann fluxes is not implemented yet");
                };

                if (dirichletOutput_()) {
                    DUNE_THROW(Dune::NotImplemented,
                               "VTK output Dirichlet conditions is not implemented yet");
                };
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (primaryVarsOutput_()) this->commitPriVarsBuffer_(writer, "PV_%s", primaryVars_);
        if (boundaryTypesOutput_()) this->commitScalarBuffer_(writer, "BC_%s", boundaryTypes_);
        if (neumannOutput_()) this->commitEqBuffer_(writer, "neumann_%s", neumann_);
        if (dirichletOutput_()) this->commitEqBuffer_(writer, "dirichlet_%s", dirichlet_);
    }

private:
    static bool primaryVarsOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WritePrimaryVars); };

    static bool boundaryTypesOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteBoundaryTypes); };

    static bool neumannOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteNeumann); };

    static bool dirichletOutput_()
    { return GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteDirichlet); };

    EqBuffer primaryVars_;
    ScalarBuffer boundaryTypes_;
    EqBuffer neumann_;
    EqBuffer dirichlet_;
};

}

#endif
