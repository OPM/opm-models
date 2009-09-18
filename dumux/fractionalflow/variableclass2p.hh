// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
#ifndef DUNE_VARIABLECLASS2P_NEW_HH
#define DUNE_VARIABLECLASS2P_NEW_HH

#include <dune/istl/bvector.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations
 * @author Markus Wolff
 */

namespace Dune
{
/*!
 * \ingroup fracflow
 * \ingroup diffusion
 * \ingroup transport
 */
//! Class including the variables and data of discretized data of the constitutive relations.
/*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 */
template<class GridView, class Scalar>
class VariableClass
{
private:
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wetting = 0, nonWetting = 1
    };

typedef    typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > ScalarVectorType;//!<type for vector of scalars
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,2> > PhasePropVectorType;//!<type for vector of phase properties
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,2> > FluidPropVectorType;//!<type for vector of fluid properties
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, 2>, 2*dim> > PotType;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars

private:
    GridView& gridViewDiffusion_;
    GridView& gridViewTransport_;
    const IndexSet& indexSetDiffusion_;
    const IndexSet& indexSetTransport_;
    const int gridSizeDiffusion_;
    const int gridSizeTransport_;

    bool multiscale_;
    const int codim_;

    Scalar time_;

    ScalarVectorType saturation_;
    ScalarVectorType pressure_;
    PhasePropVectorType mobility_;//store lambda for efficiency reasons
    PhasePropVectorType fracFlowFunc_;
    ScalarVectorType capillaryPressure_;
    VelType velocity_;
    VelType velocitySecondPhase_;

    PotType potential_;

    FluidPropVectorType density_;
    FluidPropVectorType viscosity_;

    ScalarVectorType volumecorrection_;

//    ScalarVectorType elementVolumes_;
//    ScalarVectorType Srn_;

public:
    //! Constructs a VariableClass object
    /**
     *  @param gridViewDiff a DUNE gridview object corresponding to the diffusion equation
     *  @param gridViewTrans a DUNE gridview object corresponding to the transport equation
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(GridView& gridViewDiff, GridView& gridViewTrans, Scalar& initialSat = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridViewDiff), gridViewTransport_(gridViewTrans),
    indexSetDiffusion_(gridViewDiff.indexSet()),indexSetTransport_(gridViewTrans.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(0)),gridSizeTransport_(indexSetTransport_.size(0)), multiscale_(true), codim_(0), time_(0)
    {
        initializeGlobalVariablesDiffPart(initialVel);
        initializeGlobalVariablesTransPart(initialVel);

//        analyzeMassInitialize();
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(GridView& gridView, Scalar& initialSat = *(new Scalar(1)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridView), gridViewTransport_(gridView),
    indexSetDiffusion_(gridView.indexSet()),indexSetTransport_(gridView.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(0)),gridSizeTransport_(indexSetTransport_.size(0)), multiscale_(false), codim_(0), time_(0)
    {
        initializeGlobalVariablesDiffPart(initialVel);
        initializeGlobalVariablesTransPart(initialSat);

//        analyzeMassInitialize();
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param codim codimension of the entity of which data has to be strored
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(GridView& gridView, int codim, Scalar& initialSat = *(new Scalar(1)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridView), gridViewTransport_(gridView),
    indexSetDiffusion_(gridView.indexSet()),indexSetTransport_(gridView.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(codim)),gridSizeTransport_(indexSetTransport_.size(codim)), multiscale_(false), codim_(codim), time_(0)
    {
        initializeGlobalVariablesDiffPart(initialVel);
        initializeGlobalVariablesTransPart(initialSat);

//        analyzeMassInitialize();
    }
private:
    void initializeGlobalVariablesDiffPart(Dune::FieldVector<Scalar, dim>& initialVel)
    {
        //resize to grid size
        pressure_.resize(gridSizeDiffusion_);
        velocity_.resize(gridSizeDiffusion_);//depends on pressure
        velocitySecondPhase_.resize(gridSizeDiffusion_);//depends on pressure
        potential_.resize(gridSizeDiffusion_);//depends on pressure
        density_.resize(gridSizeDiffusion_);//depends on pressure
        viscosity_.resize(gridSizeDiffusion_);//depends on pressure

        //initialise variables
        pressure_ = 0;
        velocity_ = initialVel;
        velocitySecondPhase_ = initialVel;
        initializePotentials(initialVel);
        density_=0;
        viscosity_=0;
    }
    void initializeGlobalVariablesTransPart(int initialSat)
    {
        //resize to grid size
        saturation_.resize(gridSizeTransport_);
        mobility_.resize(gridSizeTransport_);//lambda is dependent on saturation! ->choose same size
        fracFlowFunc_.resize(gridSizeTransport_);//depends on saturation
        capillaryPressure_.resize(gridSizeTransport_);//depends on saturation
        volumecorrection_.resize(gridSizeTransport_);//dS/dt for correction of pressure equation

        //initialise variables
        saturation_ = initialSat;
        mobility_ = 0;
        fracFlowFunc_ = 0;
        capillaryPressure_ = 0;
        volumecorrection_=0;
    }
    void initializePotentials (Dune::FieldVector<Scalar, dim>& initialVel)
    {
        if (initialVel.two_norm())
        {
            // compute update vector
            ElementIterator eItEnd = gridViewTransport_.template end<0>();
            for (ElementIterator eIt = gridViewTransport_.template begin<0>(); eIt != eItEnd; ++eIt)
            {
                // cell index
                int globalIdxI = indexSetTransport_.index(*eIt);

                // run through all intersections with neighbors and boundary
                IntersectionIterator
                isItEnd = gridViewTransport_.template iend(*eIt);
                for (IntersectionIterator
                        isIt = gridViewTransport_.template ibegin(*eIt); isIt
                        !=isItEnd; ++isIt)
                {
                    // local number of facet
                    int indexInInside = isIt->indexInInside();

                    // get geometry type of face
                    Dune::GeometryType faceGT = isIt->geometryInInside().type();

                    // center in face's reference element
                    const Dune::FieldVector<Scalar,dim-1>&
                    faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                    Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

                    potential_[globalIdxI][indexInInside][wetting] = initialVel*unitOuterNormal;
                    potential_[globalIdxI][indexInInside][nonWetting] = initialVel*unitOuterNormal;
                }
            }
        }
        else
        {
            potential_ = FieldVector<Scalar,2>(0);
        }
        return;
    }
//    void analyzeMassInitialize()
//    {
//        elementVolumes_.resize(gridSizeTransport_);
//        Srn_.resize(gridSizeTransport_);
//        Srn_=0;
//        ElementIterator eItEnd = gridViewTransport_.template end<0>();
//        for (ElementIterator eIt = gridViewTransport_.template begin<0>(); eIt != eItEnd; ++eIt)
//        elementVolumes_[indexSetTransport_.index(*eIt)] = (*eIt).geometry().volume();
//    }
//    void analyzeMass()
//    {
//        Scalar totalMass = 0;
//        Scalar trappedMass = 0;
//        for (int i = 0; i < gridSizeTransport_; i++)
//        {
//            totalMass += saturation_[i]*elementVolumes_[i]*0.15*density_[i][1];//fixed porosity 0.15
//            trappedMass += 0.15*elementVolumes_[i]*Srn_[i]*density_[i][1];
//
//            //            if (Srn_[i] > (saturation_[i] + 1e-4))
//            //                DUNE_THROW(MathError, "Srn = " << Srn_[i] << " is greater than Sn = " << (saturation_[i]));
//        }
//        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
//        std::cout.precision(3);
//        std::cout << time_ << ": Mass non-wetting Phase: " << totalMass << " kg "
//        << "Residually trapped: "<< trappedMass/totalMass*100.0 << "%." << std::endl;
//    }
    //Write saturation and pressure into file
    void vtkoutMultiLevel(const char* name, int k) const
    {
        if (codim_ == 0)
        {
            if (multiscale_)
            {
                Dune::VTKWriter<GridView> vtkwriterpressure(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-press%05d", name, k);
                vtkwriterpressure.addCellData(pressure_, "pressure");
                vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

                Dune::VTKWriter<GridView> vtkwritersaturation(gridViewTransport_);
                sprintf(fname, "%s-%05d", name, k);
                vtkwritersaturation.addCellData(saturation_, "saturation");
                vtkwritersaturation.write(fname, VTKOptions::ascii);
            }
            else
            {
                VTKWriter<GridView> vtkwriter(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-%05d", name, k);
                vtkwriter.addCellData(pressure_, "pressure");
                vtkwriter.addCellData(saturation_, "saturation");
                ScalarVectorType densityW(gridSizeDiffusion_);
                ScalarVectorType densityNW(gridSizeDiffusion_);
                ScalarVectorType viscosityW(gridSizeDiffusion_);
                ScalarVectorType viscosityNW(gridSizeDiffusion_);

                for (int i=0;i< gridSizeDiffusion_;i++)
                {
                    densityW[i] = density_[i][wetting];
                    densityNW[i] = density_[i][nonWetting];
                    viscosityW[i] = viscosity_[i][wetting];
                    viscosityNW[i] = viscosity_[i][nonWetting];
                }
                vtkwriter.addCellData(densityW, "wetting phase density");
                vtkwriter.addCellData(densityNW, "nonwetting phase density");
                vtkwriter.addCellData(viscosityW, "wetting phase viscosity");
                vtkwriter.addCellData(viscosityNW, "nonwetting phase viscosity");

                vtkwriter.write(fname, VTKOptions::ascii);
            }
        }
        if (codim_ == dim)
        {
            if (multiscale_)
            {
                Dune::VTKWriter<GridView> vtkwriterpressure(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-press%05d", name, k);
                vtkwriterpressure.addVertexData(pressure_, "pressure");
                vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

                Dune::VTKWriter<GridView> vtkwritersaturation(gridViewTransport_);
                sprintf(fname, "%s-%05d", name, k);
                vtkwritersaturation.addVertexData(saturation_, "saturation");
                vtkwritersaturation.write(fname, VTKOptions::ascii);
            }
            else
            {
                VTKWriter<GridView> vtkwriter(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-%05d", name, k);
                vtkwriter.addVertexData(pressure_, "pressure");
                vtkwriter.addVertexData(saturation_, "saturation");

                vtkwriter.write(fname, VTKOptions::ascii);
            }
        }

        return;
    }
public:
    //! Return saturation vector
    ScalarVectorType& saturation()
    {
        return saturation_;
    }

    //! Return pressure vector
    ScalarVectorType& pressure()
    {
        return pressure_;
    }

    //! Return velocity vector
    VelType& velocity()
    {
        return velocity_;
    }

    //! Return velocity vector
    VelType& velocitySecondPhase()
    {
        return velocitySecondPhase_;
    }

    //! Return vector of wetting phase potential gradients
    Scalar& potentialWetting(int Idx1,int Idx2)
    {
        return potential_[Idx1][Idx2][wetting];
    }

    //! Return vector of non-wetting phase potential gradients
    Scalar& potentialNonWetting(int Idx1,int Idx2)
    {
        return potential_[Idx1][Idx2][nonWetting];
    }

    //! Return vector of wetting phase mobilities
    Scalar& mobilityWetting(int Idx)
    {
        return mobility_[Idx][wetting];
    }

    //! Return vector of non-wetting phase mobilities
    Scalar& mobilityNonWetting(int Idx)
    {
        return mobility_[Idx][nonWetting];
    }

    //! Return vector of wetting phase fractional flow functions
    Scalar& fracFlowFuncWetting(int Idx)
    {
        return fracFlowFunc_[Idx][wetting];
    }

    //! Return vector of non-wetting phase fractional flow functions
    Scalar& fracFlowFuncNonWetting(int Idx)
    {
        return fracFlowFunc_[Idx][nonWetting];
    }

    //! Return capillary pressure vector
    Scalar& capillaryPressure(int Idx)
    {
        return capillaryPressure_[Idx][0];
    }

    //! Return density vector
    Scalar& densityWetting(int Idx)
    {
        return density_[Idx][wetting];
    }
    //! Return density vector
    Scalar& densityNonWetting(int Idx)
    {
        return density_[Idx][nonWetting];
    }

    //! Return density vector
    Scalar& viscosityWetting(int Idx)
    {
        return viscosity_[Idx][wetting];
    }

    //! Return density vector
    Scalar& viscosityNonWetting(int Idx)
    {
        return viscosity_[Idx][nonWetting];
    }

    Scalar& volumecorrection(int Idx)
    {
        return volumecorrection_[Idx][0];
    }

    //! Return current time
    Scalar& time()
    {
        return time_;
    }

    void updateTime(Scalar dt)
    {
        time_ += dt;
    }

//    void storeSrn(Scalar Srn, int index)
//    {
//        Srn_[index]=Srn;
//    }

    //! Get index of element (codim 0 entity) corresponding to the grid of the discretized diffusion equation.
    /*! Get index of element (codim 0 entity) corresponding to the grid of the discretized diffusion equation.
     * @param element codim 0 entity
     * \return element index
     */
    int indexDiffusion(const Element& element)
    {
        return indexSetDiffusion_.index(element);
    }

    //! Get index of element (codim 0 entity) corresponding to the grid of the discretized transport equation.
    /*! Get index of element (codim 0 entity) corresponding to the grid of the discretized transport equation.
     * @param element codim 0 entity
     * \return element index
     */
    int indexTransport(const Element& element)
    {
        return indexSetTransport_.index(element);
    }

    //!Return the number of data elements of the discretized diffusion equation
    int gridSizeDiffusion()
    {
        return gridSizeDiffusion_;
    }

    //!Return the number of data elements of the discretized transport equation
    int gridSizeTransport()
    {
        return gridSizeTransport_;
    }

    //!Return gridView on the grid of the discretized diffusion equation
    GridView& gridViewDiffusion()
    {
        return gridViewDiffusion_;
    }

    //!Return gridView on the grid of the discretized transport equation
    GridView& gridViewTransport()
    {
        return gridViewTransport_;
    }

    //! Get saturation
    /*! evaluate saturation at given element
     @param  element      entity of codim 0
     \return     value of saturation
     */
    const Dune::FieldVector<Scalar,1>& satElement(const Element& element) const
    {
        return saturation_[indexSetTransport_.index(element)];;
    }

    //! Get pressure
    /*! evaluate pressure at given element
     @param  element      entity of codim 0
     \return     value of pressure
     */
    const Dune::FieldVector<Scalar,1>& pressElement(const Element& element) const
    {
        return pressure_[indexSetDiffusion_.index(element)];
    }

    //! Get velocity at given element face
    /*! evaluate velocity at given location
     @param  element      entity of codim 0
     @param  indexInInside     index in reference element
     \return     vector of velocity
     */
    const Dune::FieldVector<Scalar,dim>& vTotalElementFace(const Element& element,
            const int indexInInside) const
    {
        int elemId = indexSetTransport_.index(element);

        return (velocity_[elemId][indexInInside]);
    }

    //! \brief Write data files
    /*!
     *  \param name file name
     *  \param k format parameter
     */
    void vtkout(const char* name, int k) const
    {
//        analyzeMass();
        vtkoutMultiLevel(name, k);
    }
};
}
#endif
