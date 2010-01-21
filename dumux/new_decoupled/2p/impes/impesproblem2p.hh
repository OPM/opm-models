// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for all problems which use the box scheme
 */
#ifndef DUMUX_IMPESPROBLEM_2P_HH
#define DUMUX_IMPESPROBLEM_2P_HH

#include <dumux/new_decoupled/common/impes.hh>
#include <dumux/new_decoupled/common/impesproblem.hh>
#include <dumux/new_decoupled/2p/variableclass2p.hh>
#include <dumux/new_decoupled/2p/2pproperties.hh>


namespace Dune
{
/*!
 * \ingroup Decoupled
 * \brief  Base class for all 2-phase problems which use an impes algorithm
 *
 * \todo Please doc me more!
 */
template<class TypeTag, class Implementation>
class IMPESProblem2P : public IMPESProblem<TypeTag, Implementation>
{
    typedef IMPESProblem<TypeTag, Implementation> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid                         Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;

    // material properties
//    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem))       FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase))    WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonwettingPhase)) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters))            SpatialParameters;


    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld>      GlobalPosition;

public:
    IMPESProblem2P(const GridView &gridView)
        : ParentType(gridView),
        gravity_(0),spatialParameters_(gridView)
    {
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = - 9.81;
    }

    IMPESProblem2P(const GridView &gridView, WettingPhase& wPhase, NonwettingPhase& nPhase)
        : ParentType(gridView),
        gravity_(0), wPhase_(wPhase), nPhase_(nPhase),
        spatialParameters_(gridView)
    {
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = - 9.81;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { return asImp_()->temperature(); };

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    /*!
     * \brief Fluid properties of the wetting phase.
     */
    const WettingPhase &wettingPhase() const
    { return wPhase_; }

    /*!
     * \brief Fluid properties of the non-wetting phase.
     */
    const NonwettingPhase &nonwettingPhase() const
    { return nPhase_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParameters &spatialParameters()
    { return spatialParameters_; }

    /*!
     * \copydoc spatialParameters()
     */
    const SpatialParameters &spatialParameters() const
    { return spatialParameters_; }


    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }

    GlobalPosition  gravity_;

    // fluids and material properties
    WettingPhase    wPhase_;
    NonwettingPhase nPhase_;
    SpatialParameters  spatialParameters_;
};

}

#endif
