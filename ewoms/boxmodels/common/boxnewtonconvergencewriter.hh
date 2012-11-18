// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::BoxNewtonConvergenceWriter
 */
#ifndef EWOMS_BOX_NEWTON_CONVERGENCE_WRITER_HH
#define EWOMS_BOX_NEWTON_CONVERGENCE_WRITER_HH

#include <ewoms/io/vtkmultiwriter.hh>
#include <ewoms/common/propertysystem.hh>

namespace Ewoms {
namespace Properties {
// forward declaration of the required property tags
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(NewtonMethod);
NEW_PROP_TAG(SolutionVector);
NEW_PROP_TAG(GlobalEqVector);
}

/*!
 * \ingroup BoxModel
 * \ingroup Newton
 *
 * \brief Writes the intermediate solutions during the Newton scheme
 *        for models using the box scheme
 */
template <class TypeTag>
class BoxNewtonConvergenceWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;

    typedef Ewoms::VtkMultiWriter<GridView>  VtkMultiWriter;

public:
    BoxNewtonConvergenceWriter(NewtonMethod &nm)
        : newtonMethod_(nm)
    {
        timeStepIdx_ = 0;
        iteration_ = 0;
        vtkMultiWriter_ = 0;
    }

    ~BoxNewtonConvergenceWriter()
    { delete vtkMultiWriter_; }

    /*!
     * \brief Called by the Newton method before the actual algorithm
     *        is started for any given timestep.
     */
    void beginTimeStep()
    {
        ++timeStepIdx_;
        iteration_ = 0;
    }

    /*!
     * \brief Called by the Newton method before an iteration of the
     *        Newton algorithm is started.
     */
    void beginIteration()
    {
        ++ iteration_;
        if (!vtkMultiWriter_)
            vtkMultiWriter_ = new VtkMultiWriter(newtonMethod_.problem().gridView(), "convergence");
        vtkMultiWriter_->beginWrite(timeStepIdx_ + iteration_ / 100.0);
    }

    /*!
     * \brief Write the Newton update to disk.
     *
     * Called after the linear solution is found for an iteration.
     *
     * \param uLastIter The solution vector of the previous iteration.
     * \param deltaU The negative difference between the solution
     *        vectors of the previous and the current iteration.
     */
    void writeFields(const SolutionVector &uLastIter,
                     const GlobalEqVector &deltaU)
    {
        try {
            newtonMethod_.problem().model().addConvergenceVtkFields(*vtkMultiWriter_, uLastIter, deltaU);
        }
        catch (...) {
            std::cout << "oops: exception thrown on rank " << newtonMethod_.problem().gridView().comm().rank() << " while writing the convergence\n";
        };

    }

    /*!
     * \brief Called by the Newton method after an iteration of the
     *        Newton algorithm has been completed.
     */
    void endIteration()
    { vtkMultiWriter_->endWrite(); }

    /*!
     * \brief Called by the Newton method after Newton algorithm
     *        has been completed for any given timestep.
     *
     * This method is called regardless of whether the Newton method
     * converged or not.
     */
    void endTimeStep()
    { iteration_ = 0; }

private:
    int timeStepIdx_;
    int iteration_;
    VtkMultiWriter *vtkMultiWriter_;
    NewtonMethod &newtonMethod_;
};

} // namespace Ewoms

#endif
