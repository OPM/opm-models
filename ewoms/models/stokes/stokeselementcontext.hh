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
 * \copydoc Ewoms::StokesElementContext
 */
#ifndef EWOMS_STOKES_ELEMENT_CONTEXT_HH
#define EWOMS_STOKES_ELEMENT_CONTEXT_HH

#include <ewoms/disc/common/fvbaseelementcontext.hh>
#include <ewoms/common/alignedallocator.hh>

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief This class stores an array of IntensiveQuantities objects, one
 *        intensive quantities object for each of the element's vertices
 */
template<class TypeTag>
class StokesElementContext : public FvBaseElementContext<TypeTag>
{
    typedef FvBaseElementContext<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    typedef typename GridView::template Codim<0>::Entity Element;

    // we don't allow copies of element contexts!
    StokesElementContext(const StokesElementContext& ) = delete;

public:
    /*!
     * \brief The constructor.
     */
    explicit StokesElementContext(const Simulator& simulator)
        : ParentType(simulator)
    { }

    static void *operator new(size_t size)
    { return Ewoms::aligned_alloc(alignof(StokesElementContext), size); }

    static void operator delete(void *ptr)
    { Ewoms::aligned_free(ptr); }

    /*!
     * \brief Compute the finite volume geometry for an element.
     *
     * \param elem The grid element for which the finite volume geometry ought to be
     *             computed.
     */
    void updateStencil(const Element& elem)
    {
        ParentType::updateStencil(elem);

#if HAVE_DUNE_LOCALFUNCTIONS
        this->stencil_.updateCenterGradients();
#else
        // center gradients require dune-localfunctions
        assert(false);
#endif
    }

    /*!
     * \brief Compute the extensive quantities of all sub-control volume
     *        faces of the current element for a single time index.
     *
     * \param timeIdx The index of the solution vector used by the
     *                time discretization.
     */
    void updateExtensiveQuantities(unsigned timeIdx)
    {
        // the history size of the time discretization in number of steps
        enum { historySize = GET_PROP_VALUE(TypeTag, TimeDiscHistorySize) };

        if (timeIdx == 0) {
            for (unsigned tIdx = 0; tIdx < historySize; ++ tIdx) {
                // update the gradients inside a sub control volume. These gradients is
                // essentially extensive quantities, but they are attributed the sub-control
                // volumes, which makes them odd fellows.
                size_t nDof = this->numDof(timeIdx);
                for (unsigned gradDofIdx = 0; gradDofIdx < nDof; gradDofIdx++) {
                    auto& iq = this->dofVars_[gradDofIdx].intensiveQuantities[tIdx];
                    iq.updateScvGradients(/*context=*/*this, gradDofIdx, tIdx);
                }
            }
        }

        ParentType::updateExtensiveQuantities(timeIdx);
    }
};

} // namespace Ewoms

#endif
