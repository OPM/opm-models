/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: andreas.lauser _at_ iws.uni-stuttgart.de                         *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_NEW_2P2C_BOX_MODEL_HH
#define DUMUX_NEW_2P2C_BOX_MODEL_HH

#include <dumux/new_models/2p2c/2p2cboxjacobianbase.hh>

namespace Dune
{

/*!
 * \brief The local jacobian operator for the isothermal two-phase,
 *        two-component model.
 *
 * This is basically just a wrapper for TwoPTwoCBoxJacobianBase so
 * that it can be instantiated.
 */
template<class TypeTag>
class TwoPTwoCBoxJacobian : public TwoPTwoCBoxJacobianBase<TypeTag,
                                                           // implementation
                                                           TwoPTwoCBoxJacobian<TypeTag> >
{
    typedef TwoPTwoCBoxJacobian<TypeTag>                   ThisType;
    typedef TwoPTwoCBoxJacobianBase<TypeTag, ThisType>     ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

public:
    TwoPTwoCBoxJacobian(Problem &problem)
        : ParentType(problem)
    {
    };
};

/**
 * \brief Isothermal two-phase two-component model.
 *
 * This implements an isothermal two phase two component
 * model. Depending on which traits are used the primary variables are
 * either $p_w$ and $S_n;X$ or $p_n$ or $S_w;X$. By default they are
 * $p_w$ and $S_n$
 */
template<class TypeTag >
class TwoPTwoCBoxModel
    : public BoxScheme<TypeTag,
                       // Implementation of the box scheme
                       TwoPTwoCBoxModel<TypeTag> >
{
    typedef TwoPTwoCBoxModel<TypeTag>                             ThisType;
    typedef BoxScheme<TypeTag, ThisType>                          ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))        Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))       Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))      GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) LocalJacobian;


    enum {
        dim = GridView::dimension
    };
    typedef typename GridView::template Codim<0>::Entity     Element;
    typedef typename GridView::template Codim<dim>::Entity   Vertex;

public:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCTraits)) TwoPTwoCTraits;

    TwoPTwoCBoxModel(Problem &prob)
        : ParentType(prob, twoPTwoCLocalJacobian_),
          twoPTwoCLocalJacobian_(prob)
    {
    }

    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailedTry()
    {
        ParentType::updateFailedTry();

        twoPTwoCLocalJacobian_.setSwitched(false);
        twoPTwoCLocalJacobian_.resetPhaseState();
        twoPTwoCLocalJacobian_.updateStaticData(this->curSolFunction(),
                                                this->prevSolFunction());
    };

    /*!
     * \brief Called by the BoxScheme's update method.
     */
    void updateSuccessful()
    {
        ParentType::updateSuccessful();

        twoPTwoCLocalJacobian_.updateOldPhaseState();
        twoPTwoCLocalJacobian_.setSwitched(false);
    }

    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        twoPTwoCLocalJacobian_.addVtkFields(writer, this->curSolFunction());
    }

    /*!
     * \brief Calculate the masses in the system for
     *        the current timestep.
     */
    void calculateMass(Dune::FieldVector<Scalar, 4> &mass)
    {
        twoPTwoCLocalJacobian_.calculateMass(this->curSolFunction(), mass);
    }

    /*!
     * \brief Returns true if there was a primary variable switch
     *        after the last time step.
     */
    bool switched() const
    { return twoPTwoCLocalJacobian_.switched(); }

    /*!
     * \brief Write the current solution to a restart file.
     */
    void serializeEntity(std::ostream &outStream,
                         const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vert);
        
        twoPTwoCLocalJacobian_.serializeEntity(outStream, vert);
    };

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     */
    void deserializeEntity(std::istream &inStream,
                           const Vertex &vert)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, vert);

        twoPTwoCLocalJacobian_.deserializeEntity(inStream, vert);
    };


private:
    // calculates the jacobian matrix at a given position
    LocalJacobian twoPTwoCLocalJacobian_;
};
}

#endif
