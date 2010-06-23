// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2007-2010 by Yufei Cao                                    *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@mathematik.uni-stuttgart.de                   *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#ifndef DUMUX_MPFAOVELOCITIES2P_UPWIND_HH
#define DUMUX_MPFAOVELOCITIES2P_UPWIND_HH

#include "fvmpfaopressure2p_new.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @author Yufei Cao
 */

namespace Dumux
{

template<class TypeTag> class FVMPFAOVelocities2P: public FVMPFAOPressure2P<TypeTag>
{
    typedef FVMPFAOVelocities2P<TypeTag> ThisType;
    typedef FVMPFAOPressure2P<TypeTag> ParentType;typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) ReferenceElements;
    typedef typename ReferenceElements::Container ReferenceElementContainer;
    typedef typename ReferenceElements::ContainerFaces ReferenceElementFaceContainer;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridTypeIndices)) GridTypeIndices;

    typedef Dumux::InteractionVolume<TypeTag> InteractionVolume;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    enum
    {
        globalCorner = 2, globalEdge = 3, neumannNeumann = 0, dirichletDirichlet = 1, dirichletNeumann = 2, neumannDirichlet = 3
    };

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

public:
    FVMPFAOVelocities2P(Problem& problem) :
        ParentType(problem)
    {
    }

    void calculateVelocity();

private:
    static const Scalar threshold_ = 1e-15;
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation)); //!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
}; // end of template

template<class TypeTag>
void FVMPFAOVelocities2P<TypeTag>::calculateVelocity()
{
    // reset velocity
    this->problem().variables().velocity() = FieldVector(0);
    this->problem().variables().velocitySecondPhase() = FieldVector(0);

    // run through all elements
    VertexIterator vItEnd = this->problem().gridView().template end<dim> ();
    for (VertexIterator vIt = this->problem().gridView().template begin<dim> (); vIt != vItEnd; ++vIt)
    {
        int globalVertIdx = this->problem().variables().index(*vIt);

        std::vector<int> bcTypeFace(2 * dim);
        bcTypeFace[0] = this->interactionVolumes_[globalVertIdx].getBoundaryType(0);
        bcTypeFace[1] = this->interactionVolumes_[globalVertIdx].getBoundaryType(1);
        bcTypeFace[2] = this->interactionVolumes_[globalVertIdx].getBoundaryType(2);
        bcTypeFace[3] = this->interactionVolumes_[globalVertIdx].getBoundaryType(3);

        if (bcTypeFace[0] == InteractionVolume::inside && bcTypeFace[1] == InteractionVolume::inside && bcTypeFace[2]
                == InteractionVolume::inside && bcTypeFace[3] == InteractionVolume::inside)
        {

            ElementPointer& elementPointer1 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
            ElementPointer& elementPointer2 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(1);
            ElementPointer& elementPointer3 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(2);
            ElementPointer& elementPointer4 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(3);

            // cell index
            int globalIdx1 = this->problem().variables().index(*elementPointer1);
            int globalIdx2 = this->problem().variables().index(*elementPointer2);
            int globalIdx3 = this->problem().variables().index(*elementPointer3);
            int globalIdx4 = this->problem().variables().index(*elementPointer4);

            // get global coordinate of cell centers
            const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
            const GlobalPosition& globalPos2 = elementPointer2->geometry().center();
            const GlobalPosition& globalPos3 = elementPointer3->geometry().center();
            const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

            //get the viscosities
            Scalar viscosityW = this->problem().variables().viscosityWetting(globalIdx1);
            Scalar viscosityNW = this->problem().variables().viscosityNonwetting(globalIdx1);

            // get pressure values
            Scalar press1 = this->problem().variables().pressure()[globalIdx1];
            Scalar press2 = this->problem().variables().pressure()[globalIdx2];
            Scalar press3 = this->problem().variables().pressure()[globalIdx3];
            Scalar press4 = this->problem().variables().pressure()[globalIdx4];

            //calculate potential gradients
            Scalar potentialW12 = 0;
            Scalar potentialNW12 = 0;
            Scalar potentialW14 = 0;
            Scalar potentialNW14 = 0;
            Scalar potentialW34 = 0;
            Scalar potentialNW34 = 0;
            Scalar potentialW32 = 0;
            Scalar potentialNW32 = 0;

            Scalar dist12 = (globalPos1 - globalPos2).two_norm();
            Scalar dist23 = (globalPos2 - globalPos3).two_norm();
            Scalar dist34 = (globalPos3 - globalPos4).two_norm();
            Scalar dist41 = (globalPos4 - globalPos1).two_norm();

            potentialW12 = (press1 - press2) / dist12;
            potentialNW12 = (press1 - press2) / dist12;
            potentialW14 = (press1 - press4) / dist41;
            potentialNW14 = (press1 - press4) / dist41;
            potentialW34 = (press3 - press4) / dist34;
            potentialNW34 = (press3 - press4) / dist34;
            potentialW32 = (press3 - press2) / dist23;
            potentialNW32 = (press3 - press2) / dist23;

            //store potentials for further calculations (saturation, ...)
            this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 0))
                    = potentialW12;
            this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 0))
                    = potentialNW12;
            this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 1))
                    = potentialW14;
            this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 1))
                    = potentialNW14;
            this->problem().variables().potentialWetting(globalIdx2, this->interactionVolumes_[globalIdx2].getIndexOnElement(1, 0))
                    = -potentialW32;
            this->problem().variables().potentialNonwetting(globalIdx2, this->interactionVolumes_[globalIdx2].getIndexOnElement(1, 0))
                    = -potentialNW32;
            this->problem().variables().potentialWetting(globalIdx2, this->interactionVolumes_[globalIdx2].getIndexOnElement(1, 1))
                    = -potentialW12;
            this->problem().variables().potentialNonwetting(globalIdx2, this->interactionVolumes_[globalIdx2].getIndexOnElement(1, 1))
                    = -potentialNW12;
            this->problem().variables().potentialWetting(globalIdx3, this->interactionVolumes_[globalIdx3].getIndexOnElement(2, 0))
                    = potentialW34;
            this->problem().variables().potentialNonwetting(globalIdx3, this->interactionVolumes_[globalIdx3].getIndexOnElement(2, 0))
                    = potentialNW34;
            this->problem().variables().potentialWetting(globalIdx3, this->interactionVolumes_[globalIdx3].getIndexOnElement(2, 1))
                    = potentialW32;
            this->problem().variables().potentialNonwetting(globalIdx3, this->interactionVolumes_[globalIdx3].getIndexOnElement(2, 1))
                    = potentialNW32;
            this->problem().variables().potentialWetting(globalIdx4, this->interactionVolumes_[globalIdx4].getIndexOnElement(3, 0))
                    = -potentialW14;
            this->problem().variables().potentialNonwetting(globalIdx4, this->interactionVolumes_[globalIdx4].getIndexOnElement(3, 0))
                    = -potentialNW14;
            this->problem().variables().potentialWetting(globalIdx4, this->interactionVolumes_[globalIdx4].getIndexOnElement(3, 1))
                    = -potentialW34;
            this->problem().variables().potentialNonwetting(globalIdx4, this->interactionVolumes_[globalIdx4].getIndexOnElement(3, 1))
                    = -potentialNW34;

            //compute mobilities of face 1
            Dune::FieldVector<Scalar, numPhases> lambda12(0.0);
            Dune::FieldVector<Scalar, numPhases> lambda21(0.0);
            if (potentialW12 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx1];
                lambda12[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda21[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda12[wPhaseIdx] /= viscosityW;
                lambda21[wPhaseIdx] /= viscosityW;
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx2];
                lambda12[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda21[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda12[wPhaseIdx] /= viscosityW;
                lambda21[wPhaseIdx] /= viscosityW;
            }
            if (potentialNW12 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx1];
                lambda12[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda21[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda12[nPhaseIdx] /= viscosityNW;
                lambda21[nPhaseIdx] /= viscosityNW;
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx2];
                lambda12[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda21[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda12[nPhaseIdx] /= viscosityNW;
                lambda21[nPhaseIdx] /= viscosityNW;
            }

            //compute mobilities of face 4
            Dune::FieldVector<Scalar, numPhases> lambda14(0.0);
            Dune::FieldVector<Scalar, numPhases> lambda41(0.0);
            if (potentialW14 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx1];
                lambda14[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda41[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
                lambda14[wPhaseIdx] /= viscosityW;
                lambda41[wPhaseIdx] /= viscosityW;
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx4];
                lambda14[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda41[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
                lambda14[wPhaseIdx] /= viscosityW;
                lambda41[wPhaseIdx] /= viscosityW;
            }
            if (potentialNW14 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx1];
                lambda14[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda41[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx4];
                lambda14[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1), sat);
                lambda41[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
                lambda14[nPhaseIdx] /= viscosityNW;
                lambda41[nPhaseIdx] /= viscosityNW;
            }

            //compute mobilities of face 2
            Dune::FieldVector<Scalar, numPhases> lambda23(0.0);
            Dune::FieldVector<Scalar, numPhases> lambda32(0.0);
            if (potentialW32 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx3];
                lambda32[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda23[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda32[wPhaseIdx] /= viscosityW;
                lambda23[wPhaseIdx] /= viscosityW;
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx2];
                lambda32[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda23[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda32[wPhaseIdx] /= viscosityW;
                lambda23[wPhaseIdx] /= viscosityW;
            }
            if (potentialNW32 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx3];
                lambda32[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda23[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda32[nPhaseIdx] /= viscosityNW;
                lambda23[nPhaseIdx] /= viscosityNW;
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx2];
                lambda32[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda23[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2), sat);
                lambda32[nPhaseIdx] /= viscosityNW;
                lambda23[nPhaseIdx] /= viscosityNW;
            }

            //compute mobilities of face 3
            Dune::FieldVector<Scalar, numPhases> lambda34(0.0);
            Dune::FieldVector<Scalar, numPhases> lambda43(0.0);
            if (potentialW34 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx3];
                lambda34[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda43[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
                lambda34[wPhaseIdx] /= viscosityW;
                lambda43[wPhaseIdx] /= viscosityW;
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx4];
                lambda34[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda43[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
                lambda34[wPhaseIdx] /= viscosityW;
                lambda43[wPhaseIdx] /= viscosityW;
            }
            if (potentialNW34 >= 0)
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx3];
                lambda34[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda43[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
                lambda34[nPhaseIdx] /= viscosityNW;
                lambda43[nPhaseIdx] /= viscosityNW;
            }
            else
            {
                Scalar sat = this->problem().variables().saturation()[globalIdx4];
                lambda34[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos3, *elementPointer3), sat);
                lambda43[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4), sat);
                lambda34[nPhaseIdx] /= viscosityNW;
                lambda43[nPhaseIdx] /= viscosityNW;
            }

            //compute total mobility of cell 1
            Scalar lambda1(this->problem().variables().mobilityWetting(globalIdx1));
            lambda1 += this->problem().variables().mobilityNonwetting(globalIdx1);

            Scalar lambda2(this->problem().variables().mobilityWetting(globalIdx2));
            lambda2 += this->problem().variables().mobilityNonwetting(globalIdx2);

            Scalar lambda3(this->problem().variables().mobilityWetting(globalIdx3));
            lambda3 += this->problem().variables().mobilityNonwetting(globalIdx3);

            Scalar lambda4(this->problem().variables().mobilityWetting(globalIdx4));
            lambda4 += this->problem().variables().mobilityNonwetting(globalIdx4);

            Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
            Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
            Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
            Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
            Scalar gn12nu23 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 0);
            Scalar gn12nu21 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 1);
            Scalar gn23nu23 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 0);
            Scalar gn23nu21 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 0, 1);
            Scalar gn43nu32 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 0, 1);
            Scalar gn43nu34 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 0, 0);
            Scalar gn23nu32 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 1, 1);
            Scalar gn23nu34 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda3, 2, 1, 0);
            Scalar gn43nu41 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 0);
            Scalar gn43nu43 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 1, 1);
            Scalar gn14nu41 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 0);
            Scalar gn14nu43 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 1);

            // compute transmissibility matrix T = CA^{-1}B+F
            Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> C(0), F(0), A(0), B(0);

            // evaluate matrix C, F, A, B
            C[0][0] = -gn12nu12;
            C[0][3] = -gn12nu14;
            C[1][0] = gn23nu21;
            C[1][1] = -gn23nu23;
            C[2][1] = gn43nu32;
            C[2][2] = gn43nu34;
            C[3][2] = -gn14nu43;
            C[3][3] = gn14nu41;

            F[0][0] = gn12nu12 + gn12nu14;
            F[1][1] = -gn23nu21 + gn23nu23;
            F[2][2] = -gn43nu34 - gn43nu32;
            F[3][3] = gn14nu43 - gn14nu41;

            A[0][0] = gn12nu12 + gn12nu21;
            A[0][1] = -gn12nu23;
            A[0][3] = gn12nu14;
            A[1][0] = -gn23nu21;
            A[1][1] = gn23nu23 + gn23nu32;
            A[1][2] = gn23nu34;
            A[2][1] = -gn43nu32;
            A[2][2] = -gn43nu34 - gn43nu43;
            A[2][3] = gn43nu41;
            A[3][0] = -gn14nu12;
            A[3][2] = gn14nu43;
            A[3][3] = -gn14nu41 - gn14nu14;

            //                        std::cout << A << "\n";

            B[0][0] = gn12nu12 + gn12nu14;
            B[0][1] = gn12nu21 - gn12nu23;
            B[1][1] = -gn23nu21 + gn23nu23;
            B[1][2] = gn23nu34 + gn23nu32;
            B[2][2] = -gn43nu34 - gn43nu32;
            B[2][3] = -gn43nu43 + gn43nu41;
            B[3][0] = -gn14nu12 - gn14nu14;
            B[3][3] = gn14nu43 - gn14nu41;

            // pressures flux calculation
            Dune::FieldVector<Scalar, 2 * dim> p(0);
            p[0] = press1;
            p[1] = press2;
            p[2] = press3;
            p[3] = press4;

            //flux vector
            Dune::FieldVector<Scalar, 2 * dim> flux(0);

            // compute T
            A.invert();
            F += C.rightmultiply(B.leftmultiply(A));
            Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> T(F);

            T.mv(p, flux);

            // evaluate parts of velocity
            FieldVector vel12 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 0);
            FieldVector vel14 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 1);
            FieldVector vel23 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(1, 0);
            FieldVector vel21 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(1, 1);
            FieldVector vel34 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(2, 0);
            FieldVector vel32 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(2, 1);
            FieldVector vel41 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(3, 0);
            FieldVector vel43 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(3, 1);

            vel12 *= flux[0] / 2;//divide by 2 because the flux is related to the half face!
            vel14 *= flux[3] / 2;
            vel23 *= flux[1] / 2;
            vel21 *= flux[0] / 2;
            vel34 *= flux[2] / 2;
            vel32 *= flux[1] / 2;
            vel41 *= flux[3] / 2;
            vel43 *= flux[2] / 2;

            for (int i = 0; i < numPhases; i++)
            {
                // evaluate parts of velocity
                FieldVector phaseVel12(vel12);
                FieldVector phaseVel14(vel14);
                FieldVector phaseVel23(vel23);
                FieldVector phaseVel21(vel21);
                FieldVector phaseVel34(vel34);
                FieldVector phaseVel32(vel32);
                FieldVector phaseVel41(vel41);
                FieldVector phaseVel43(vel43);

                switch (velocityType_)
                {
                case vw:
                {

                    switch (i)
                    {
                    case 0:
                        this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 0)]
                                += phaseVel12 *= (lambda12[i] / (lambda12[wPhaseIdx] + lambda12[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 1)]
                                += phaseVel14 *= (lambda14[i] / (lambda14[wPhaseIdx] + lambda14[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(1, 0)]
                                += phaseVel23 *= (lambda23[i] / (lambda23[wPhaseIdx] + lambda23[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(1, 1)]
                                += phaseVel21 *= (lambda21[i] / (lambda21[wPhaseIdx] + lambda21[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(2, 0)]
                                += phaseVel34 *= (lambda34[i] / (lambda34[wPhaseIdx] + lambda34[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(2, 1)]
                                += phaseVel32 *= (lambda32[i] / (lambda32[wPhaseIdx] + lambda32[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(3, 0)]
                                += phaseVel41 *= (lambda41[i] / (lambda41[wPhaseIdx] + lambda41[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(3, 1)]
                                += phaseVel43 *= (lambda43[i] / (lambda43[wPhaseIdx] + lambda43[nPhaseIdx]));
                        break;
                    case 1:
                        this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                0, 0)] += phaseVel12 *= (lambda12[i] / (lambda12[wPhaseIdx] + lambda12[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                0, 1)] += phaseVel14 *= (lambda14[i] / (lambda14[wPhaseIdx] + lambda14[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                1, 0)] += phaseVel23 *= (lambda23[i] / (lambda23[wPhaseIdx] + lambda23[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                1, 1)] += phaseVel21 *= (lambda21[i] / (lambda21[wPhaseIdx] + lambda21[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                2, 0)] += phaseVel34 *= (lambda34[i] / (lambda34[wPhaseIdx] + lambda34[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                2, 1)] += phaseVel32 *= (lambda32[i] / (lambda32[wPhaseIdx] + lambda32[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                3, 0)] += phaseVel41 *= (lambda41[i] / (lambda41[wPhaseIdx] + lambda41[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                3, 1)] += phaseVel43 *= (lambda43[i] / (lambda43[wPhaseIdx] + lambda43[nPhaseIdx]));
                        break;
                    }
                    break;
                }
                case vn:
                {
                    switch (i)
                    {
                    case 1:
                        this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 0)]
                                += phaseVel12 *= (lambda12[i] / (lambda12[wPhaseIdx] + lambda12[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 1)]
                                += phaseVel14 *= (lambda14[i] / (lambda14[wPhaseIdx] + lambda14[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(1, 0)]
                                += phaseVel23 *= (lambda23[i] / (lambda23[wPhaseIdx] + lambda23[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(1, 1)]
                                += phaseVel21 *= (lambda21[i] / (lambda21[wPhaseIdx] + lambda21[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(2, 0)]
                                += phaseVel34 *= (lambda34[i] / (lambda34[wPhaseIdx] + lambda34[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(2, 1)]
                                += phaseVel32 *= (lambda32[i] / (lambda32[wPhaseIdx] + lambda32[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(3, 0)]
                                += phaseVel41 *= (lambda41[i] / (lambda41[wPhaseIdx] + lambda41[nPhaseIdx]));
                        this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(3, 1)]
                                += phaseVel43 *= (lambda43[i] / (lambda43[wPhaseIdx] + lambda43[nPhaseIdx]));
                        break;
                    case 0:
                        this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                0, 0)] += phaseVel12 *= (lambda12[i] / (lambda12[wPhaseIdx] + lambda12[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                0, 1)] += phaseVel14 *= (lambda14[i] / (lambda14[wPhaseIdx] + lambda14[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                1, 0)] += phaseVel23 *= (lambda23[i] / (lambda23[wPhaseIdx] + lambda23[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                1, 1)] += phaseVel21 *= (lambda21[i] / (lambda21[wPhaseIdx] + lambda21[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                2, 0)] += phaseVel34 *= (lambda34[i] / (lambda34[wPhaseIdx] + lambda34[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                2, 1)] += phaseVel32 *= (lambda32[i] / (lambda32[wPhaseIdx] + lambda32[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                3, 0)] += phaseVel41 *= (lambda41[i] / (lambda41[wPhaseIdx] + lambda41[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                3, 1)] += phaseVel43 *= (lambda43[i] / (lambda43[wPhaseIdx] + lambda43[nPhaseIdx]));
                        break;
                    }
                    break;
                }
                case vt:
                {
                    if (i == 0)
                    {
                        this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                0, 0)] += phaseVel12 *= (lambda12[i] / (lambda12[wPhaseIdx] + lambda12[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                0, 1)] += phaseVel14 *= (lambda14[i] / (lambda14[wPhaseIdx] + lambda14[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                1, 0)] += phaseVel23 *= (lambda23[i] / (lambda23[wPhaseIdx] + lambda23[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                1, 1)] += phaseVel21 *= (lambda21[i] / (lambda21[wPhaseIdx] + lambda21[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                2, 0)] += phaseVel34 *= (lambda34[i] / (lambda34[wPhaseIdx] + lambda34[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                2, 1)] += phaseVel32 *= (lambda32[i] / (lambda32[wPhaseIdx] + lambda32[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                3, 0)] += phaseVel41 *= (lambda41[i] / (lambda41[wPhaseIdx] + lambda41[nPhaseIdx]));
                        this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                3, 1)] += phaseVel43 *= (lambda43[i] / (lambda43[wPhaseIdx] + lambda43[nPhaseIdx]));
                    }
                    this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 0)]
                            += phaseVel12 *= (lambda12[i] / (lambda12[wPhaseIdx] + lambda12[nPhaseIdx]));
                    this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 1)]
                            += phaseVel14 *= (lambda14[i] / (lambda14[wPhaseIdx] + lambda14[nPhaseIdx]));
                    this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(1, 0)]
                            += phaseVel23 *= (lambda23[i] / (lambda23[wPhaseIdx] + lambda23[nPhaseIdx]));
                    this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(1, 1)]
                            += phaseVel21 *= (lambda21[i] / (lambda21[wPhaseIdx] + lambda21[nPhaseIdx]));
                    this->problem().variables().velocity()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(2, 0)]
                            += phaseVel34 *= (lambda34[i] / (lambda34[wPhaseIdx] + lambda34[nPhaseIdx]));
                    this->problem().variables().velocity()[globalIdx3][this->interactionVolumes_[globalVertIdx].getIndexOnElement(2, 1)]
                            += phaseVel32 *= (lambda32[i] / (lambda32[wPhaseIdx] + lambda32[nPhaseIdx]));
                    this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(3, 0)]
                            += phaseVel41 *= (lambda41[i] / (lambda41[wPhaseIdx] + lambda41[nPhaseIdx]));
                    this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(3, 1)]
                            += phaseVel43 *= (lambda43[i] / (lambda43[wPhaseIdx] + lambda43[nPhaseIdx]));

                    break;
                }
                }
            }
        }

        // at least one face on boundary!
        else
        {
            std::vector<int> interactionVolFaces(0);
            for (int faceIdx = 0; faceIdx < 2 * dim; faceIdx++)
            {
                if (bcTypeFace[faceIdx] != InteractionVolume::outside)
                {
                    interactionVolFaces.push_back(faceIdx);
                }
            }

            int numInteractionVolFaces = interactionVolFaces.size();

            switch (numInteractionVolFaces)
            {
            case globalCorner:
            {
                ElementPointer& elementPointer1 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(0);

                // cell index
                int globalIdx1 = this->problem().variables().index(*elementPointer1);

                // get global coordinate of cell centers
                const GlobalPosition& globalPos1 = elementPointer1->geometry().center();

                //get the densities
                Scalar densityW = this->problem().variables().densityWetting(globalIdx1);
                Scalar densityNW = this->problem().variables().densityNonwetting(globalIdx1);

                //get the viscosities
                Scalar viscosityW = this->problem().variables().viscosityWetting(globalIdx1);
                Scalar viscosityNW = this->problem().variables().viscosityNonwetting(globalIdx1);

                //compute total mobility of cell 1
                Scalar sat = this->problem().variables().saturation()[globalIdx1];
                Scalar satBound = 0;

                //check boundary sat at face 1
                if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                {
                    satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                }
                else
                {
                    satBound = sat;
                }
                Scalar lambda12 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1),
                        satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                        *elementPointer1), satBound) / viscosityNW;

                //check boundary sat at face 4
                if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                {
                    satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                }
                else
                {
                    satBound = sat;
                }
                Scalar lambda14 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1),
                        satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                        *elementPointer1), satBound) / viscosityNW;

                std::vector<FieldVector> vel12(numPhases, this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 0));
                std::vector<FieldVector> vel14(numPhases, this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 1));

                //neumann - neumann
                if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::neumann)
                {
                    for (int i = 0; i < numPhases; i++)
                    {
                        // get neumann boundary value
                        std::vector<Scalar> J1(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                                BoundaryConditions::neumann> (interactionVolFaces[0]));
                        J1[wPhaseIdx] /= densityW;
                        J1[nPhaseIdx] /= densityNW;

                        std::vector<Scalar> J4(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                                BoundaryConditions::neumann> (interactionVolFaces[1]));
                        J4[wPhaseIdx] /= densityW;
                        J4[nPhaseIdx] /= densityNW;

                        vel12[i] *= J1[i] / 2;
                        vel14[i] *= J4[i] / 2;
                    }
                }
                //dirichlet - dirichlet
                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::dirichlet)
                {
                    // get pressure values
                    Scalar press1 = this->problem().variables().pressure()[globalIdx1];

                    Scalar g1 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[0]);
                    Scalar g4 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);

                    //calculate potential gradients
                    Scalar potentialW12 = (press1 - g1);
                    Scalar potentialNW12 = (press1 - g1);
                    Scalar potentialW14 = (press1 - g4);
                    Scalar potentialNW14 = (press1 - g4);

                    //store potentials for further calculations (saturation, ...)
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            0)) = potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 0)) = potentialNW12;
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            1)) = potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 1)) = potentialNW14;

                    //compute mobilities of face 1
                    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                    if (potentialW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 4
                    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                    if (potentialW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }

                    Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 1);
                    Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 0);
                    Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 1);
                    Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 0);

                    for (int i = 0; i < numPhases; i++)
                    {
                        vel12[i] *= ((gn12nu12 + gn12nu14) * press1 - gn12nu12 * g1 - gn12nu14 * g4) / 2;//divide by 2 because the flux is related to the half face!
                        vel14[i] *= ((gn14nu12 + gn14nu14) * press1 - gn14nu12 * g1 - gn14nu14 * g4) / 2;

                        vel12[i] *= lambda12Upw[i] / (lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx]);
                        vel14[i] *= lambda12Upw[i] / (lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx]);
                    }
                }
                //neumann - dirichlet
                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::dirichlet)
                {
                    // get neumann boundary value
                    std::vector<Scalar> J1(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                            BoundaryConditions::neumann> (interactionVolFaces[0]));
                    J1[wPhaseIdx] /= densityW;
                    J1[nPhaseIdx] /= densityNW;

                    // get dirichlet boundary value
                    Scalar g4 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);

                    // get pressure values
                    Scalar press1 = this->problem().variables().pressure()[globalIdx1];

                    //calculate potential gradients
                    Scalar potentialW12 = (J1[0]);
                    Scalar potentialNW12 = (J1[1]);
                    Scalar potentialW14 = (press1 - g4);
                    Scalar potentialNW14 = (press1 - g4);

                    //                    std::cout<<"potentials 12 = "<<potentialW12<<", "<<potentialNW12<<
                    //                            "\npotentials 14 = "<<potentialW14<<", "<<potentialNW14<<"\n";

                    //store potentials for further calculations (saturation, ...)
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            0)) = potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 0)) = potentialNW12;
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            1)) = potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 1)) = potentialNW14;

                    //compute mobilities of face 1
                    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                    if (potentialW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 4
                    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                    if (potentialW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //                                                std::cout<<"lambda12 = "<<lambda12<<
                    //                                                        "\nlambda14 = "<<lambda14<<"\n";

                    Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 1);
                    Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 0);
                    Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 1);
                    Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 0);

                    //                        vel12[i] *= this->interactionVolumes_[globalVertIdx].getFaceArea(0, 0);
                    //std::cout<<this->interactionVolumes_[globalVertIdx].getFaceArea(0, 0)<<"\n";

                    for (int i = 0; i < numPhases; i++)
                    {
                        vel12[i] *= J1[i] / 2;

                        vel14[i] *= ((gn14nu14 - gn14nu12 * gn12nu14 / gn12nu12) * press1 + (gn14nu12 * gn12nu14 / gn12nu12 - gn14nu14)
                                * g4 + gn14nu12 / gn12nu12 * J1[i] * this->interactionVolumes_[globalVertIdx].getFaceArea(0, 0)) / 2;//divide by 2 because the flux is related to the half face!

                        vel14[i] *= lambda12Upw[i] / (lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx]);
                    }
                }
                //dirichlet - neumann
                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::neumann)
                {
                    // get dirichlet boundary value
                    Scalar g1 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[0]);

                    // get pressure values
                    Scalar press1 = this->problem().variables().pressure()[globalIdx1];

                    // get neumann boundary value
                    std::vector<Scalar> J4(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                            BoundaryConditions::neumann> (interactionVolFaces[1]));
                    J4[wPhaseIdx] /= densityW;
                    J4[nPhaseIdx] /= densityNW;

                    //calculate potential gradients
                    Scalar potentialW12 = (press1 - g1);
                    Scalar potentialNW12 = (press1 - g1);
                    Scalar potentialW14 = (J4[0]);
                    Scalar potentialNW14 = (J4[1]);

                    //                    std::cout<<"potentials 12 = "<<potentialW12<<", "<<potentialNW12<<
                    //                            "\npotentials 14 = "<<potentialW14<<", "<<potentialNW14<<"\n";

                    //store potentials for further calculations (saturation, ...)
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            0)) = potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 0)) = potentialNW12;
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            1)) = potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 1)) = potentialNW14;

                    //compute mobilities of face 1
                    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                    if (potentialW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 4
                    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                    if (potentialW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }

//                                                                    std::cout<<"lambda12 = "<<lambda12<<
//                                                                            "\nlambda14 = "<<lambda14<<"\n";

                    Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 1);
                    Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 0);
                    Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 1);
                    Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 0);

                    for (int i = 0; i < numPhases; i++)
                    {
                        vel12[i] *= ((gn12nu12 - gn12nu14 * gn14nu12 / gn14nu14) * press1 + (gn12nu14 * gn14nu12 / gn14nu14 - gn12nu12)
                                * g1 + gn12nu14 / gn14nu14 * J4[i] * this->interactionVolumes_[globalVertIdx].getFaceArea(0, 1)) / 2;//divide by 2 because the flux is related to the half face!


                        //                        vel14[i] *= this->interactionVolumes_[globalVertIdx].getFaceArea(0, 1);
                        vel14[i] *= J4[i] / 2;

                        vel12[i] *= lambda12Upw[i] / (lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx]);
                    }
                }
                else
                {
                    std::cout << interactionVolFaces[0] << ", " << interactionVolFaces[1] << ", " << interactionVolFaces[2] << ", "
                            << interactionVolFaces[3] << "\n";
                    DUNE_THROW(Dune::NotImplemented, "Boundary combination not supported in MPFA implementation");
                }

                for (int i = 0; i < numPhases; i++)
                {
                    switch (velocityType_)
                    {
                    case vw:
                    {
                        switch (i)
                        {
                        case 0:
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += vel12[i];
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += vel14[i];
                            break;
                        case 1:
                            this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += vel12[i];
                            this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += vel14[i];
                            break;
                        }
                        break;
                    }
                    case vn:
                    {
                        switch (i)
                        {
                        case 1:
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += vel12[i];
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += vel14[i];
                            break;
                        case 0:
                            this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += vel12[i];
                            this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += vel14[i];
                            break;
                        }
                        break;
                    }
                    case vt:
                    {
                        if (i == 0)
                        {
                            this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += vel12[i];
                            this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += vel14[i];
                        }

                        this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 0)]
                                += vel12[i];
                        this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(0, 1)]
                                += vel14[i];

                        break;
                    }
                    }
                }

                break;
            }
            case globalEdge:
            {
                //numbering assures that either faces 2 and 4 or faces 1 and 3 are on the boundary
                //neumann - neumann

                if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::neumann)
                {
                    ElementPointer& elementPointer1 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer4 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(3);

                    // cell index
                    int globalIdx1 = this->problem().variables().index(*elementPointer1);
                    int globalIdx4 = this->problem().variables().index(*elementPointer4);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

                    // get pressure values
                    Scalar press1 = this->problem().variables().pressure()[globalIdx1];
                    Scalar press4 = this->problem().variables().pressure()[globalIdx4];

                    //get the densities
                    Scalar densityW = this->problem().variables().densityWetting(globalIdx1);
                    Scalar densityNW = this->problem().variables().densityNonwetting(globalIdx1);

                    //get the viscosities
                    Scalar viscosityW = this->problem().variables().viscosityWetting(globalIdx1);
                    Scalar viscosityNW = this->problem().variables().viscosityNonwetting(globalIdx1);

                    // get neumann boundary value
                    std::vector<Scalar> J1(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                            BoundaryConditions::neumann> (interactionVolFaces[0]));
                    J1[wPhaseIdx] /= densityW;
                    J1[nPhaseIdx] /= densityNW;

                    std::vector<Scalar> J4(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                            BoundaryConditions::neumann> (interactionVolFaces[1]));
                    J4[wPhaseIdx] /= densityW;
                    J4[nPhaseIdx] /= densityNW;

                    //calculate potential gradients
                    Scalar potentialW12 = (J1[0]);
                    Scalar potentialNW12 = (J1[1]);
                    Scalar potentialW14 = (press1 - press4);
                    Scalar potentialNW14 = (press1 - press4);
                    Scalar potentialW43 = (J4[0]);
                    Scalar potentialNW43 = (J4[1]);

                    //                    std::cout<<"potentialW12 = "<<potentialW12<<", potentialNW12 = "<<potentialNW12<<
                    //                            "\npotentialW14 = "<<potentialW14<<", potentialNW14 = "<<potentialNW14<<
                    //                            "\npotentialW43 = "<<potentialW43<<", potentialNW43 = "<<potentialNW43<<"\n";

                    //store potentials for further calculations (saturation, ...)
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            0)) = potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 0)) = potentialNW12;
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            1)) = potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 1)) = potentialNW14;
                    this->problem().variables().potentialWetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(3,
                            0)) = -potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            3, 0)) = -potentialNW14;
                    this->problem().variables().potentialWetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(3,
                            1)) = potentialW43;
                    this->problem().variables().potentialNonwetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            3, 1)) = potentialNW43;

                    //compute mobilities of face 1
                    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                    if (potentialW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 3
                    Dune::FieldVector<Scalar, numPhases> lambda43Upw(0.0);

                    if (potentialW43 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda43Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(2))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(2);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx4];
                        }
                        lambda43Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW43 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda43Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(2))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(2);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx4];
                        }
                        lambda43Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 4
                    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                    Dune::FieldVector<Scalar, numPhases> lambda41Upw(0.0);
                    if (potentialW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                        lambda41Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                        lambda41Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                        lambda41Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                        lambda41Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute total mobility of cell 1
                    Scalar lambda1(this->problem().variables().mobilityWetting(globalIdx1));
                    lambda1 += this->problem().variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda4(this->problem().variables().mobilityWetting(globalIdx4));
                    lambda4 += this->problem().variables().mobilityNonwetting(globalIdx4);

                    //compute total mobility at the boundaries
                    Scalar satBound = 0;
                    //check boundary sat at face 1
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[0]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[0]);
                    }
                    else
                    {
                        satBound = this->problem().variables().saturation()[globalIdx1];
                    }
                    Scalar lambda12 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                            *elementPointer1), satBound) / viscosityNW;

                    //check boundary sat at face 4
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[1]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[1]);
                    }
                    else
                    {
                        satBound = this->problem().variables().saturation()[globalIdx4];
                    }
                    Scalar lambda43 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                            *elementPointer4), satBound) / viscosityNW;

                    // evaluate parts of velocity
                    FieldVector vel14 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 1);
                    FieldVector vel41 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(3, 0);

                    Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 1);
                    Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 0);
                    Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
                    Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
                    Scalar gn43nu41 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda43, 3, 1, 0);
                    Scalar gn43nu43 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda43, 3, 1, 1);
                    Scalar gn14nu41 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 0);
                    Scalar gn14nu43 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 1);

                    // compute transmissibility matrix T = CA^{-1}B+F
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 1> C(0), A(0);
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> F(0), B(0);
                    Dune::FieldVector<Scalar, 2 * dim - 1> fB(0);

                    // evaluate matrix C, F, A, B
                    C[0][0] = -gn12nu12;
                    C[0][2] = -gn12nu14;
                    C[1][1] = -gn43nu43;
                    C[1][2] = gn43nu41;
                    C[2][1] = -gn14nu43;
                    C[2][2] = gn14nu41;

                    F[0][0] = gn12nu12 + gn12nu14;
                    F[1][1] = gn43nu43 - gn43nu41;
                    F[2][1] = gn14nu43 - gn14nu41;

                    A[0][0] = gn12nu12;
                    A[0][2] = gn12nu14;
                    A[1][1] = gn43nu43;
                    A[1][2] = -gn43nu41;
                    A[2][0] = -gn14nu12;
                    A[2][1] = gn14nu43;
                    A[2][2] = -gn14nu41 - gn14nu14;

                    B[0][0] = gn12nu12 + gn12nu14;
                    B[1][1] = gn43nu43 - gn43nu41;
                    B[2][0] = -gn14nu12 - gn14nu14;
                    B[2][1] = gn14nu43 - gn14nu41;

                    fB[0] = J1[wPhaseIdx] + J1[nPhaseIdx];
                    fB[1] = J4[wPhaseIdx] + J4[nPhaseIdx];

                    // use the pressure values to compute the fluxes
                    Dune::FieldVector<Scalar, 2 * dim - 2> p(0);
                    p[0] = press1;
                    p[1] = press4;

                    Dune::FieldVector<Scalar, 2 * dim - 1> flux(0);

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);
                    // compute T
                    A.invert();

                    //calc RHS r
                    C.rightmultiplyany(A).mv(fB, r);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> AinvB(A.rightmultiplyany(B));
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(AinvB.leftmultiplyany(C));
                    T += F;

                    T.mv(p, flux);

                    flux += r;

                    vel14 *= flux[2] / 2;//divide by 2 because the flux is related to the half face!
                    vel41 *= flux[2] / 2;

                    for (int i = 0; i < numPhases; i++)
                    {
                        FieldVector phaseVel12 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 0);
                        FieldVector phaseVel43 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(3, 1);
                        phaseVel12 *= J1[i] / 2;
                        phaseVel43 *= J4[i] / 2;

                        // evaluate parts of velocity
                        FieldVector phaseVel14(vel14);
                        FieldVector phaseVel41(vel41);

                        phaseVel14  *= (lambda14Upw[i] / (lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx]));
                        phaseVel41  *= (lambda14Upw[i] / (lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx]));

                        switch (velocityType_)
                        {
                        case vw:
                        {
                            switch (i)
                            {
                            case 0:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            case 1:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            }
                            break;
                        }
                        case vn:
                        {
                            switch (i)
                            {
                            case 1:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            case 0:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            }
                            break;
                        }
                        case vt:
                        {
                            if (i == 0)
                            {
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                            }

                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += phaseVel12;
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += phaseVel14;
                            this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    3, 0)] += phaseVel41;
                            this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    3, 1)] += phaseVel43;

                            break;
                        }
                        }
                    }

                }
                //neumann - neumann

                else if (bcTypeFace[interactionVolFaces[1]] == BoundaryConditions::neumann && bcTypeFace[interactionVolFaces[2]]
                        == BoundaryConditions::neumann)
                {
                    ElementPointer& elementPointer1 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer2 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(1);

                    // cell index
                    int globalIdx1 = this->problem().variables().index(*elementPointer1);
                    int globalIdx2 = this->problem().variables().index(*elementPointer2);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos2 = elementPointer2->geometry().center();

                    // get pressure values
                    Scalar press1 = this->problem().variables().pressure()[globalIdx1];
                    Scalar press2 = this->problem().variables().pressure()[globalIdx2];

                    //get the viscosities
                    Scalar viscosityW = this->problem().variables().viscosityWetting(globalIdx1);
                    Scalar viscosityNW = this->problem().variables().viscosityNonwetting(globalIdx1);

                    //get the densities
                    Scalar densityW = this->problem().variables().densityWetting(globalIdx1);
                    Scalar densityNW = this->problem().variables().densityNonwetting(globalIdx1);

                    // get neumann boundary value
                    std::vector<Scalar> J1(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                            BoundaryConditions::neumann> (interactionVolFaces[1]));
                    J1[wPhaseIdx] /= densityW;
                    J1[nPhaseIdx] /= densityNW;

                    std::vector<Scalar> J2(this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<
                            BoundaryConditions::neumann> (interactionVolFaces[2]));
                    J2[wPhaseIdx] /= densityW;
                    J2[nPhaseIdx] /= densityNW;

                    //calculate potential gradients
                    Scalar potentialW12 = (press1 - press2);
                    Scalar potentialNW12 = (press1 - press2);
                    Scalar potentialW14 = (J1[0]);
                    Scalar potentialNW14 = (J1[1]);
                    Scalar potentialW23 = (J2[0]);
                    Scalar potentialNW23 = (J2[1]);

                    //                    std::cout<<"potentialW12 = "<<potentialW12<<", potentialNW12 = "<<potentialNW12<<
                    //                            "\npotentialW14 = "<<potentialW14<<", potentialNW14 = "<<potentialNW14<<
                    //                            "\npotentialW23 = "<<potentialW23<<", potentialNW23 = "<<potentialNW23<<"\n";

                    //store potentials for further calculations (saturation, ...)
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            0)) = potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 0)) = potentialNW12;
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            1)) = potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 1)) = potentialNW14;
                    this->problem().variables().potentialWetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(1,
                            0)) = potentialW23;
                    this->problem().variables().potentialNonwetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            1, 0)) = potentialNW23;
                    this->problem().variables().potentialWetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(1,
                            1)) = -potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            1, 1)) = -potentialNW12;

                    //compute mobilities of face 1
                    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                    Dune::FieldVector<Scalar, numPhases> lambda21Upw(0.0);
                    if (potentialW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                        lambda21Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                        lambda21Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                        lambda21Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                        lambda21Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 2
                    Dune::FieldVector<Scalar, numPhases> lambda23Upw(0.0);
                    if (potentialW23 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda23Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(1))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(1);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx2];
                        }
                        lambda23Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW23 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda23Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(1))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(1);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx2];
                        }
                        lambda23Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 4
                    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                    if (potentialW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute total mobility of cell 1
                    Scalar lambda1(this->problem().variables().mobilityWetting(globalIdx1));
                    lambda1 += this->problem().variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda2(this->problem().variables().mobilityWetting(globalIdx2));
                    lambda2 += this->problem().variables().mobilityNonwetting(globalIdx2);

                    //compute total mobility at the boundaries
                    Scalar satBound = 0;
                    //check boundary sat at face 1
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[1]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[1]);
                    }
                    else
                    {
                        satBound = this->problem().variables().saturation()[globalIdx1];
                    }
                    Scalar lambda14 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                            *elementPointer1), satBound) / viscosityNW;

                    //check boundary sat at face 4
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[2]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[2]);
                    }
                    else
                    {
                        satBound = this->problem().variables().saturation()[globalIdx2];
                    }
                    Scalar lambda23 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                            *elementPointer2), satBound) / viscosityNW;

                    // initialize velocities with normal divided by face area
                    FieldVector vel12 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 0);
                    FieldVector vel21 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(1, 1);

                    Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
                    Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
                    Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 1);
                    Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 0);
                    Scalar gn12nu23 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 0);
                    Scalar gn12nu21 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 1);
                    Scalar gn23nu23 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda23, 1, 0, 0);
                    Scalar gn23nu21 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda23, 1, 0, 1);

                    // compute transmissibility matrix T = CA^{-1}B+F
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 1> C(0), A(0);
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> F(0), B(0);
                    Dune::FieldVector<Scalar, 2 * dim - 1> fB(0);

                    // evaluate matrix C, F, A, B
                    C[0][0] = -gn12nu12;
                    C[0][2] = -gn12nu14;
                    C[1][0] = gn23nu21;
                    C[1][1] = -gn23nu23;
                    C[2][0] = -gn14nu12;
                    C[2][2] = -gn14nu14;

                    F[0][0] = gn12nu12 + gn12nu14;
                    F[1][1] = -gn23nu21 + gn23nu23;
                    F[2][0] = gn14nu12 + gn14nu14;

                    A[0][0] = gn12nu12 + gn12nu21;
                    A[0][1] = -gn12nu23;
                    A[0][2] = gn12nu14;
                    A[1][0] = -gn23nu21;
                    A[1][1] = gn23nu23;
                    A[2][0] = gn14nu12;
                    A[2][2] = gn14nu14;

                    //                                        std::cout << A << "\n";

                    B[0][0] = gn12nu12 + gn12nu14;
                    B[0][1] = gn12nu21 - gn12nu23;
                    B[1][1] = -gn23nu21 + gn23nu23;
                    B[2][0] = gn14nu12 + gn14nu14;

                    // get neumann boundary value
                    fB[1] = J2[wPhaseIdx] + J2[nPhaseIdx];
                    fB[2] = J1[wPhaseIdx] + J1[nPhaseIdx];

                    // use the pressure values to compute the fluxes
                    Dune::FieldVector<Scalar, 2 * dim - 2> p(0);
                    p[0] = press1;
                    p[1] = press2;

                    Dune::FieldVector<Scalar, 2 * dim - 1> flux(0);

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);
                    // compute T
                    A.invert();

                    //calc RHS r
                    C.rightmultiplyany(A).mv(fB, r);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> AinvB(A.rightmultiplyany(B));
                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(AinvB.leftmultiplyany(C));
                    T += F;

                    T.mv(p, flux);

                    flux += r;

                    vel12 *= flux[0] / 2;//divide by 2 because the flux is related to the half face!
                    vel21 *= flux[0] / 2;

                    //                        std::cout<<"flux2 = "<<flux<<"\n";
                    for (int i = 0; i < numPhases; i++)
                    {
                        FieldVector phaseVel14 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 1);
                        FieldVector phaseVel23 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(1, 0);
                        phaseVel14 *= J1[i] / 2;
                        phaseVel23 *= J2[i] / 2;

                        // evaluate parts of velocity
                        FieldVector phaseVel12(vel12);
                        FieldVector phaseVel21(vel21);

                        phaseVel12 *= lambda12Upw[i] / (lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx]);
                        phaseVel21 *= lambda21Upw[i] / (lambda21Upw[wPhaseIdx] + lambda21Upw[nPhaseIdx]);

                        switch (velocityType_)
                        {
                        case vw:
                        {
                            switch (i)
                            {
                            case 0:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            case 1:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            }
                            break;
                        }
                        case vn:
                        {
                            switch (i)
                            {
                            case 1:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            case 0:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            }
                            break;
                        }
                        case vt:
                        {
                            if (i == 0)
                            {
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                            }

                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += phaseVel12;
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += phaseVel14;
                            this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    1, 0)] += phaseVel23;
                            this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    1, 1)] += phaseVel21;

                            break;
                        }
                        }
                    }
                }
                //dirichlet- dirichlet

                else if (bcTypeFace[interactionVolFaces[0]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[1]]
                        == BoundaryConditions::dirichlet)
                {
                    ElementPointer& elementPointer1 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer4 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(3);

                    // cell index
                    int globalIdx1 = this->problem().variables().index(*elementPointer1);
                    int globalIdx4 = this->problem().variables().index(*elementPointer4);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos4 = elementPointer4->geometry().center();

                    //get the viscosities
                    Scalar viscosityW = this->problem().variables().viscosityWetting(globalIdx1);
                    Scalar viscosityNW = this->problem().variables().viscosityNonwetting(globalIdx1);

                    // get pressure values
                    Scalar press1 = this->problem().variables().pressure()[globalIdx1];
                    Scalar press4 = this->problem().variables().pressure()[globalIdx4];

                    Scalar g1 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[0]);
                    Scalar g4 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);

                    //calculate potential gradients
                    Scalar potentialW12 = (press1 - g1);
                    Scalar potentialNW12 = (press1 - g1);
                    Scalar potentialW14 = (press1 - press4);
                    Scalar potentialNW14 = (press1 - press4);
                    Scalar potentialW43 = (press4 - g4);
                    Scalar potentialNW43 = (press4 - g4);

                    //store potentials for further calculations (saturation, ...)
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            0)) = potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 0)) = potentialNW12;
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            1)) = potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 1)) = potentialNW14;
                    this->problem().variables().potentialWetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(3,
                            0)) = -potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            3, 0)) = -potentialNW14;
                    this->problem().variables().potentialWetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(3,
                            1)) = potentialW43;
                    this->problem().variables().potentialNonwetting(globalIdx4, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            3, 1)) = potentialNW43;

                    //compute mobilities of face 1
                    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                    if (potentialW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(0))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(0);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 3
                    Dune::FieldVector<Scalar, numPhases> lambda43Upw(0.0);

                    if (potentialW43 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda43Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(2))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(2);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx4];
                        }
                        lambda43Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW43 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda43Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(2))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(2);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx4];
                        }
                        lambda43Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda43Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 4
                    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                    Dune::FieldVector<Scalar, numPhases> lambda41Upw(0.0);
                    if (potentialW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                        lambda41Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                        lambda41Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                        lambda41Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx4];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda41Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                                *elementPointer4), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                        lambda41Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute total mobility of cell 1
                    Scalar lambda1(this->problem().variables().mobilityWetting(globalIdx1));
                    lambda1 += this->problem().variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda4(this->problem().variables().mobilityWetting(globalIdx4));
                    lambda4 += this->problem().variables().mobilityNonwetting(globalIdx4);

                    //compute total mobility at the boundaries
                    Scalar sat = this->problem().variables().saturation()[globalIdx1];
                    Scalar satBound = 0;
                    //check boundary sat at face 1
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[0]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[0]);
                    }
                    else
                    {
                        satBound = sat;
                    }
                    Scalar lambda12 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                            *elementPointer1), satBound) / viscosityNW;

                    //check boundary sat at face 4
                    sat = this->problem().variables().saturation()[globalIdx4];
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[1]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[1]);
                    }
                    else
                    {
                        satBound = sat;
                    }
                    Scalar lambda43 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos4, *elementPointer4),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos4,
                            *elementPointer4), satBound) / viscosityNW;

                    //                                        std::cout<<"lambda121 = "<<lambda12<<"\nlambda14"<<lambda14<<"\nlambda41 = "<<lambda41<<"\nlambda34 = "<<lambda43<<"\n";

                    Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 1);
                    Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda12, 0, 0, 0);
                    Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 1);
                    Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 1, 0);
                    Scalar gn43nu41 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda43, 3, 1, 0);
                    Scalar gn43nu43 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda43, 3, 1, 1);
                    Scalar gn14nu41 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 0);
                    Scalar gn14nu43 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda4, 3, 0, 1);

                    // evaluate parts of velocity
                    FieldVector vel12 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 0);
                    FieldVector vel14 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 1);
                    FieldVector vel41 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(3, 0);
                    FieldVector vel43 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(3, 1);

                    Dune::FieldVector<Scalar, 2 * dim - 1> flux(0);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(0);

                    Scalar b = 1 / (gn14nu41 + gn14nu14);

                    T[0][0] = gn12nu12 + gn12nu14 - gn12nu14 * b * (gn14nu12 + gn14nu14);
                    T[0][1] = -gn12nu14 * b * (gn14nu41 - gn14nu43);
                    T[1][0] = gn43nu41 * b * (gn14nu12 + gn14nu14);
                    T[1][1] = gn43nu43 - gn43nu41 + gn43nu41 * b * (gn14nu41 - gn14nu43);
                    T[2][0] = gn14nu14 * b * (gn14nu12 + gn14nu14);
                    T[2][1] = gn14nu43 - gn14nu41 + gn14nu41 * b * (gn14nu41 - gn14nu43);

                    //                    std::cout << "Dirichlet-dirichlet T: " << T << "\n";

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);

                    r[0] = -gn12nu12 * g1 - gn12nu14 * b * (gn14nu43 * g4 - gn14nu12 * g1);
                    r[1] = -gn43nu43 * g4 + gn43nu41 * b * (gn14nu43 * g4 - gn14nu12 * g1);
                    r[2] = -gn14nu43 * g4 + gn14nu41 * b * (gn14nu43 * g4 - gn14nu12 * g1);

                    //                    std::cout << "Dirichlet-dirichlet r: " << r << "\n";

                    // use the pressure values to compute the fluxes
                    Dune::FieldVector<Scalar, 2 * dim - 2> p(0);
                    p[0] = press1;
                    p[1] = press4;

                    T.mv(p, flux);

                    flux += r;

                    vel12 *= flux[0] / 2;//divide by 2 because the flux is related to the half face!
                    vel14 *= flux[2] / 2;
                    vel41 *= flux[2] / 2;
                    vel43 *= flux[1] / 2;

                    for (int i = 0; i < numPhases; i++)
                    {
                        // evaluate parts of velocity
                        FieldVector phaseVel12(vel12);
                        FieldVector phaseVel14(vel14);
                        FieldVector phaseVel41(vel41);
                        FieldVector phaseVel43(vel43);

                        phaseVel12 *= lambda12Upw[i] / (lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx]);
                        phaseVel14 *= lambda14Upw[i] / (lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx]);
                        phaseVel41 *= lambda41Upw[i] / (lambda41Upw[wPhaseIdx] + lambda41Upw[nPhaseIdx]);
                        phaseVel43 *= lambda43Upw[i] / (lambda43Upw[wPhaseIdx] + lambda43Upw[nPhaseIdx]);

                        switch (velocityType_)
                        {
                        case vw:
                        {
                            switch (i)
                            {
                            case 0:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            case 1:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            }
                            break;
                        }
                        case vn:
                        {
                            switch (i)
                            {
                            case 1:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            case 0:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                                break;
                            }
                            break;
                        }
                        case vt:
                        {
                            if (i == 0)
                            {
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 0)] += phaseVel41;
                                this->problem().variables().velocitySecondPhase()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        3, 1)] += phaseVel43;
                            }

                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += phaseVel12;
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += phaseVel14;
                            this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    3, 0)] += phaseVel41;
                            this->problem().variables().velocity()[globalIdx4][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    3, 1)] += phaseVel43;

                            break;

                            break;
                        }
                        }
                    }
                }
                //dirichlet - dirichlet

                else if (bcTypeFace[interactionVolFaces[1]] == BoundaryConditions::dirichlet && bcTypeFace[interactionVolFaces[2]]
                        == BoundaryConditions::dirichlet)
                {
                    ElementPointer& elementPointer1 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(0);
                    ElementPointer& elementPointer2 = this->interactionVolumes_[globalVertIdx].getSubVolumeElement(1);

                    // cell index
                    int globalIdx1 = this->problem().variables().index(*elementPointer1);
                    int globalIdx2 = this->problem().variables().index(*elementPointer2);

                    // get global coordinate of cell centers
                    const GlobalPosition& globalPos1 = elementPointer1->geometry().center();
                    const GlobalPosition& globalPos2 = elementPointer2->geometry().center();

                    //get the viscosities
                    Scalar viscosityW = this->problem().variables().viscosityWetting(globalIdx1);
                    Scalar viscosityNW = this->problem().variables().viscosityNonwetting(globalIdx1);

                    // get pressure values
                    Scalar press1 = this->problem().variables().pressure()[globalIdx1];
                    Scalar press2 = this->problem().variables().pressure()[globalIdx2];

                    Scalar g1 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[1]);
                    Scalar g2 = this->interactionVolumes_[globalVertIdx].template getBoundaryCondition<BoundaryConditions::dirichlet> (
                            interactionVolFaces[2]);

                    //calculate potential gradients
                    Scalar potentialW12 = (press1 - press2);
                    Scalar potentialNW12 = (press1 - press2);
                    Scalar potentialW14 = (press1 - g1);
                    Scalar potentialNW14 = (press1 - g1);
                    Scalar potentialW23 = (press2 - g2);
                    Scalar potentialNW23 = (press2 - g2);

                    //store potentials for further calculations (saturation, ...)
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            0)) = potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 0)) = potentialNW12;
                    this->problem().variables().potentialWetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(0,
                            1)) = potentialW14;
                    this->problem().variables().potentialNonwetting(globalIdx1, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            0, 1)) = potentialNW14;
                    this->problem().variables().potentialWetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(1,
                            0)) = potentialW23;
                    this->problem().variables().potentialNonwetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            1, 0)) = potentialNW23;
                    this->problem().variables().potentialWetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(1,
                            1)) = -potentialW12;
                    this->problem().variables().potentialNonwetting(globalIdx2, this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                            1, 1)) = -potentialNW12;

                    //compute mobilities of face 1
                    Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                    Dune::FieldVector<Scalar, numPhases> lambda21Upw(0.0);
                    if (potentialW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                        lambda21Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda12Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[wPhaseIdx] /= viscosityW;
                        lambda21Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW12 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                        lambda21Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda12Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda21Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda12Upw[nPhaseIdx] /= viscosityNW;
                        lambda21Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 2
                    Dune::FieldVector<Scalar, numPhases> lambda23Upw(0.0);
                    if (potentialW23 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda23Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(1))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(1);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx2];
                        }
                        lambda23Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW23 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx2];
                        lambda23Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(1))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(1);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx2];
                        }
                        lambda23Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                                *elementPointer2), sat);
                        lambda23Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute mobilities of face 4
                    Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                    if (potentialW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[wPhaseIdx] = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[wPhaseIdx] /= viscosityW;
                    }
                    if (potentialNW14 >= 0)
                    {
                        Scalar sat = this->problem().variables().saturation()[globalIdx1];
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }
                    else
                    {
                        Scalar sat = 0;
                        if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(3))
                        {
                            sat = this->interactionVolumes_[globalVertIdx].getDirichletSat(3);
                        }
                        else
                        {
                            sat = this->problem().variables().saturation()[globalIdx1];
                        }
                        lambda14Upw[nPhaseIdx] = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                                *elementPointer1), sat);
                        lambda14Upw[nPhaseIdx] /= viscosityNW;
                    }

                    //compute total mobility of cell 1
                    Scalar lambda1(this->problem().variables().mobilityWetting(globalIdx1));
                    lambda1 += this->problem().variables().mobilityNonwetting(globalIdx1);

                    Scalar lambda2(this->problem().variables().mobilityWetting(globalIdx2));
                    lambda2 += this->problem().variables().mobilityNonwetting(globalIdx2);

                    //compute total mobility at the boundaries
                    Scalar sat = this->problem().variables().saturation()[globalIdx1];
                    Scalar satBound = 0;
                    //check boundary sat at face 1
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[1]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[1]);
                    }
                    else
                    {
                        satBound = sat;
                    }
                    Scalar lambda14 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos1, *elementPointer1),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos1,
                            *elementPointer1), satBound) / viscosityNW;

                    //check boundary sat at face 4
                    sat = this->problem().variables().saturation()[globalIdx2];
                    if (this->interactionVolumes_[globalVertIdx].isDirichletSatBound(interactionVolFaces[2]))
                    {
                        satBound = this->interactionVolumes_[globalVertIdx].getDirichletSat(interactionVolFaces[2]);
                    }
                    else
                    {
                        satBound = sat;
                    }
                    Scalar lambda23 = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos2, *elementPointer2),
                            satBound) / viscosityW + MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos2,
                            *elementPointer2), satBound) / viscosityNW;

                    //                                        std::cout<<"lambda122 = "<<lambda12<<"\nlambda21"<<lambda21<<"\nlambda14 = "<<lambda14<<"\nlambda23 = "<<lambda23<<"\n";

                    Scalar gn12nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 1);
                    Scalar gn12nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda1, 0, 0, 0);
                    Scalar gn14nu14 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 1);
                    Scalar gn14nu12 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda14, 0, 1, 0);
                    Scalar gn12nu23 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 0);
                    Scalar gn12nu21 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda2, 1, 1, 1);
                    Scalar gn23nu23 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda23, 1, 0, 0);
                    Scalar gn23nu21 = this->interactionVolumes_[globalVertIdx].getNTKKrNu_by_dF(lambda23, 1, 0, 1);

                    // evaluate parts of velocity
                    FieldVector vel12 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 0);
                    FieldVector vel14 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(0, 1);
                    FieldVector vel23 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(1, 0);
                    FieldVector vel21 = this->interactionVolumes_[globalVertIdx].getUnitOuterNormalByFace(1, 1);

                    Dune::FieldVector<Scalar, 2 * dim - 1> flux(0);

                    Dune::FieldMatrix<Scalar, 2 * dim - 1, 2 * dim - 2> T(0);

                    Scalar b = 1 / (gn12nu12 + gn12nu21);

                    T[0][0] = gn12nu12 + gn12nu14 - gn12nu12 * b * (gn12nu14 + gn12nu12);
                    T[0][1] = -gn12nu12 * b * (gn12nu21 - gn12nu23);
                    T[1][0] = gn23nu21 * b * (gn12nu14 + gn12nu12);
                    T[1][1] = -gn23nu21 + gn23nu23 + gn23nu21 * b * (gn12nu21 - gn12nu23);
                    T[2][0] = gn14nu12 + gn14nu14 - gn14nu12 * b * (gn12nu14 + gn12nu12);
                    T[2][1] = -gn14nu12 * b * (gn12nu21 - gn12nu23);

                    //                    std::cout << "Dirichlet-dirichlet T: " << T << "\n";

                    Dune::FieldVector<Scalar, 2 * dim - 1> r(0);

                    r[0] = -gn12nu14 * g1 - gn12nu12 * b * (gn12nu23 * g2 - gn12nu14 * g1);
                    r[1] = -gn23nu23 * g2 + gn23nu21 * b * (gn12nu23 * g2 - gn12nu14 * g1);
                    r[2] = -gn14nu14 * g1 - gn14nu12 * b * (gn12nu23 * g2 - gn12nu14 * g1);

                    //                    std::cout << "Dirichlet-dirichlet r: " << r << "\n";

                    // use the pressure values to compute the fluxes
                    Dune::FieldVector<Scalar, 2 * dim - 2> p(0);
                    p[0] = press1;
                    p[1] = press2;

                    T.mv(p, flux);

                    flux += r;

                    vel12 *= flux[0] / 2;//divide by 2 because the flux is related to the half face!
                    vel14 *= flux[2] / 2;
                    vel23 *= flux[1] / 2;
                    vel21 *= flux[0] / 2;

                    for (int i = 0; i < numPhases; i++)
                    {
                        // evaluate parts of velocity
                        FieldVector phaseVel12(vel12);
                        FieldVector phaseVel14(vel14);
                        FieldVector phaseVel23(vel23);
                        FieldVector phaseVel21(vel21);

                        phaseVel12 *= lambda12Upw[i] / (lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx]);
                        phaseVel14 *= lambda14Upw[i] / (lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx]);
                        phaseVel23 *= lambda23Upw[i] / (lambda23Upw[wPhaseIdx] + lambda23Upw[nPhaseIdx]);
                        phaseVel21 *= lambda21Upw[i] / (lambda21Upw[wPhaseIdx] + lambda21Upw[nPhaseIdx]);

                        switch (velocityType_)
                        {
                        case vw:
                        {
                            switch (i)
                            {
                            case 0:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            case 1:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            }
                            break;
                        }
                        case vn:
                        {
                            switch (i)
                            {
                            case 1:
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            case 0:
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                                break;
                            }
                            break;
                        }
                        case vt:
                        {
                            if (i == 0)
                            {
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 0)] += phaseVel12;
                                this->problem().variables().velocitySecondPhase()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        0, 1)] += phaseVel14;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 0)] += phaseVel23;
                                this->problem().variables().velocitySecondPhase()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                        1, 1)] += phaseVel21;
                            }

                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 0)] += phaseVel12;
                            this->problem().variables().velocity()[globalIdx1][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    0, 1)] += phaseVel14;
                            this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    1, 0)] += phaseVel23;
                            this->problem().variables().velocity()[globalIdx2][this->interactionVolumes_[globalVertIdx].getIndexOnElement(
                                    1, 1)] += phaseVel21;

                            break;
                        }
                        }
                    }
                }
                else
                {
                    DUNE_THROW(Dune::NotImplemented, "Boundary combination not supported in MPFA implementation");
                }

                break;
            }
            }

        } // end boundaries

    } // end vertex iterator
//    printvector(std::cout, this->problem().variables().velocity(), "velocity", "row", 4, 1, 3);
//    printvector(std::cout, this->problem().variables().velocitySecondPhase(), "velocity second phase", "row", 4, 1, 3);
    return;
} // end method calcTotalVelocity

}
// end of Dune namespace
#endif
