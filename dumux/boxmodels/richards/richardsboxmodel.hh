/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
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
#ifndef DUMUX_RICHARDS_BOX_MODEL_HH
#define DUMUX_RICHARDS_BOX_MODEL_HH

#include <dumux/boxmodels/boxscheme/boxscheme.hh>

#include "richardsboxjacobian.hh"

namespace Dune
{
/*!
 * \brief Adaption of the BOX scheme to the isothermal richards model.
 */
template<class TypeTag >
class RichardsBoxModel : public BoxScheme<TypeTag,  RichardsBoxModel<TypeTag> >
{
    typedef RichardsBoxModel<TypeTag>          ThisType;
    typedef BoxScheme<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))        Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))         Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian))  LocalJacobian;

public:
    RichardsBoxModel(Problem &prob)
        : ParentType(prob, onePLocalJacobian_),
          onePLocalJacobian_(prob)
    {
    }

    /*!
     * \brief All relevant primary and secondary of the current
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer) 
    {
        onePLocalJacobian_.addVtkFields(writer, this->curSolFunction());
    }

private:
    // calculates the jacobian matrix at a given position
    LocalJacobian  onePLocalJacobian_;
};
}

#endif
