// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
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
 * \brief spatial parameters for the RichardsLensProblem
 */
#ifndef DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_LENS_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/Mp/2padapter.hh>

#include <dumux/material/fluidmatrixinteractions/Mp/2padapter.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class RichardsLensSpatialParameters;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(RichardsLensSpatialParameters);

// Set the spatial parameters
SET_TYPE_PROP(RichardsLensSpatialParameters, SpatialParameters, Dumux::RichardsLensSpatialParameters<TypeTag>);

// Set the material Law
SET_PROP(RichardsLensSpatialParameters, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;

    // ... and one for absolute saturations.
    typedef EffToAbsLaw<EffectiveLaw> TwoPMaterialLaw;
    
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };

public:
    typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};
}

/*!
 * \ingroup RichardsProblems
 * \brief The spatial parameters for the RichardsLensProblem
 */
template<class TypeTag>
class RichardsLensSpatialParameters : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP(TypeTag, ParameterTree) Params;
    typedef typename Grid::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dimWorld,dimWorld> Tensor;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    //! The parameters of the material law to be used
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief Constructor
     *
     * \param gridView The DUNE GridView representing the spatial
     *                 domain of the problem.
     */
    RichardsLensSpatialParameters(const GridView& gridView)
        : ParentType(gridView)
    {
        lensPorosity_ = Params::tree().template get<double>("Soil.FinePorosity");
        outerPorosity_ = Params::tree().template get<double>("Soil.CoarsePorosity");

        for(int i = 0; i < dim; i++)
        {
            lensK_[i][i] = Params::tree().template get<double>("Soil.FinePermeability");
            outerK_[i][i] = Params::tree().template get<double>("Soil.CoarsePermeability");
        }

        // residual saturations
        lensMaterialParams_.setSwr(Params::tree().template get<double>("Soil.FineResidualSaturationWetting"));
        lensMaterialParams_.setSnr(Params::tree().template get<double>("Soil.FineResidualSaturationNonWetting"));
        outerMaterialParams_.setSwr(Params::tree().template get<double>("Soil.CoarseResidualSaturationWetting"));
        outerMaterialParams_.setSnr(Params::tree().template get<double>("Soil.CoarseResidualSaturationNonWetting"));

        // parameters for the Van Genuchten law
        // alpha and n
//        lensMaterialParams_.setVgAlpha(0.00045);
//        lensMaterialParams_.setVgN(7.3);
//        outerMaterialParams_.setVgAlpha(0.0037);
//        outerMaterialParams_.setVgN(4.7);

        // parameters for the Brooks-Corey law
        lensMaterialParams_.setPe(Params::tree().template get<double>("Soil.FineBrooksCoreyEntryPressure"));
        lensMaterialParams_.setLambda(Params::tree().template get<double>("Soil.FineBrooksCoreyLambda"));
        outerMaterialParams_.setPe(Params::tree().template get<double>("Soil.CoarseBrooksCoreyEntryPressure"));
        outerMaterialParams_.setLambda(Params::tree().template get<double>("Soil.CoarseBrooksCoreyLambda"));

        // parameters for the linear law
        // minimum and maximum pressures
//        lensMaterialParams_.setEntryPC(0);
//        outerMaterialParams_.setEntryPC(0);
//        lensMaterialParams_.setMaxPC(0);
//        outerMaterialParams_.setMaxPC(0);
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    const Tensor &intrinsicPermeability(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &globalPos = context.pos(spaceIdx, timeIdx);
        if (isInLens_(globalPos))
            return lensPorosity_;
        return outerPorosity_;
    }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * \param element An arbitrary DUNE Codim<0> entity of the grid view
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &globalPos = context.pos(spaceIdx, timeIdx);
        if (isInLens_(globalPos))
            return lensMaterialParams_;
        return outerMaterialParams_;
    }

    /*!
     * \brief Set the bounding box of the low-permeability lens
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param lensLowerLeft the lower left corner of the lens
     * \param lensUpperRight the upper right corner of the lens
     */
    void setLensCoords(const GlobalPosition& lensLowerLeft,
                       const GlobalPosition& lensUpperRight)
    {
        lensLowerLeft_ = lensLowerLeft;
        lensUpperRight_ = lensUpperRight;
    }

private:
    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dim; ++i) {
            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Tensor lensK_;
    Tensor outerK_;
    Scalar lensPorosity_;
    Scalar outerPorosity_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams outerMaterialParams_;
};

} // end namespace

#endif

