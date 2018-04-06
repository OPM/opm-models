// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::BlackOilNewtonMethod
 */
#ifndef EWOMS_BLACK_OIL_NEWTON_METHOD_HH
#define EWOMS_BLACK_OIL_NEWTON_METHOD_HH

#include "blackoilproperties.hh"

#include <ewoms/common/signum.hh>

#include <opm/material/common/Unused.hpp>

namespace Ewoms {

/*!
 * \ingroup BlackOilModel
 *
 * \brief A newton solver which is specific to the black oil model.
 */
template <class TypeTag>
class BlackOilNewtonMethod : public GET_PROP_TYPE(TypeTag, DiscNewtonMethod)
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscNewtonMethod) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;

    static const unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    BlackOilNewtonMethod(Simulator& simulator) : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the immiscible model.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
    }

    /*!
     * \brief Returns the number of degrees of freedom for which the
     *        interpretation has changed for the most recent iteration.
     */
    unsigned numPriVarsSwitched() const
    { return numPriVarsSwitched_; }

protected:
    friend NewtonMethod<TypeTag>;
    friend ParentType;

    /*!
     * \copydoc FvBaseNewtonMethod::beginIteration_
     */
    void beginIteration_()
    {
        numPriVarsSwitched_ = 0;
        ParentType::beginIteration_();
    }

    /*!
     * \copydoc FvBaseNewtonMethod::endIteration_
     */
    void endIteration_(SolutionVector& uCurrentIter,
                       const SolutionVector& uLastIter)
    {
#if HAVE_MPI
        // in the MPI enabled case we need to add up the number of DOF
        // for which the interpretation changed over all processes.
        int localSwitched = numPriVarsSwitched_;
        MPI_Allreduce(&localSwitched,
                      &numPriVarsSwitched_,
                      /*num=*/1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);
#endif // HAVE_MPI

        this->simulator_.model().newtonMethod().endIterMsg()
            << ", num switched=" << numPriVarsSwitched_;

        ParentType::endIteration_(uCurrentIter, uLastIter);
    }

    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        const auto& comm = this->simulator_.gridView().comm();

        int succeeded;
        try {
            ParentType::update_(nextSolution,
                                currentSolution,
                                solutionUpdate,
                                currentResidual);
            succeeded = 1;
        }
        catch (...) {
            std::cout << "Newton update threw an exception on rank "
                      << comm.rank() << "\n";
            succeeded = 0;
        }
        succeeded = comm.min(succeeded);

        if (!succeeded)
            throw Opm::NumericalIssue("A process did not succeed in adapting the primary variables");

        numPriVarsSwitched_ = comm.sum(numPriVarsSwitched_);
    }

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     */
    void updatePrimaryVariables_(unsigned globalDofIdx,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual)
    {
        currentValue.checkDefined();
        Opm::Valgrind::CheckDefined(update);
        Opm::Valgrind::CheckDefined(currentResidual);

        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            // calculate the update of the current primary variable. For the
            // black-oil model we limit the pressure and saturation updates, but do
            // we not clamp anything after the specified number of iterations was
            // reached
            Scalar delta = update[eqIdx];

            // limit changes in water saturation to 20%
            if (eqIdx == Indices::waterSaturationIdx
                && std::abs(delta) > 0.2)
            {
                delta = Ewoms::signum(delta)*0.2;
            }
            else if (eqIdx == Indices::compositionSwitchIdx) {
                // the switching primary variable for composition is tricky because the
                // "reasonable" value ranges it exhibits vary widely depending on its
                // interpretation (it can represent Sg, Rs or Rv).  we limit changes in
                // gas saturation to 20%
                if (currentValue.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg
                    && std::abs(delta) > 0.2)
                {
                    delta = Ewoms::signum(delta)*0.2;
                }
            }

            // do the actual update
            nextValue[eqIdx] = currentValue[eqIdx] - delta;
        }

        // switch the new primary variables to something which is physically meaningful
        if (nextValue.adaptPrimaryVariables(this->problem(), globalDofIdx))
            ++ numPriVarsSwitched_;

        nextValue.checkDefined();
    }

private:
    int numPriVarsSwitched_;
};
} // namespace Ewoms

#endif
