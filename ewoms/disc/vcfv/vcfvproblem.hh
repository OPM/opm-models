// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \copydoc Ewoms::VcfvProblem
 */
#ifndef EWOMS_VCFV_PROBLEM_HH
#define EWOMS_VCFV_PROBLEM_HH

#include "vcfvproperties.hh"

#include <ewoms/io/vtkmultiwriter.hh>
#include <ewoms/io/restart.hh>
#include <dune/common/fvector.hh>

#include <iostream>
#include <limits>
#include <string>

namespace Ewoms {

/*!
 * \ingroup VcfvModel
 *
 * \brief Base class for all problems which use the VCVF discretization.
 *
 * \note All quantities are specified assuming a threedimensional
 *       world. Problems discretized using 2D grids are assumed to be
 *       extruded by \f$1 m\f$ and 1D grids are assumed to have a
 *       cross section of \f$1m \times 1m\f$.
 */
template<class TypeTag>
class VcfvProblem
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;

    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GridView::Grid::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    // copying a problem is not a good idea
    VcfvProblem(const VcfvProblem &);

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     *
     * \param timeManager The time manager of the simulation
     * \param gridView The view on the DUNE grid which ought to be
     *                 used (normally the leaf grid view)
     */
    VcfvProblem(TimeManager &timeManager, const GridView &gridView)
        : gridView_(gridView)
        , bboxMin_(std::numeric_limits<double>::max())
        , bboxMax_(-std::numeric_limits<double>::max())
        , elementMapper_(gridView)
        , vertexMapper_(gridView)
        , timeManager_(&timeManager)
        , newtonMethod_(asImp_())
    {
        init_();
    }

    ~VcfvProblem()
    { delete resultWriter_; }

    /*!
     * \brief Registers all available parameters for the problem and
     *        the model.
     */
    static void registerParameters()
    {
        Model::registerParameters();
        REGISTER_PARAM(TypeTag, Scalar, MaxTimeStepSize, "The maximum size to which all time steps are limited to [s]");
        REGISTER_PARAM(TypeTag, Scalar, MinTimeStepSize, "The minimum size to which all time steps are limited to [s]");
        REGISTER_PARAM(TypeTag, unsigned, MaxTimeStepDivisions, "The maximum number of divisions by two of the timestep size before the simulation bails out");
    }

    /*!
     * \brief Called by the Ewoms::TimeManager in order to
     *        initialize the problem.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        // set the initial condition of the model
        model().init(asImp_());

        assembleTime_ = 0.0;
        solveTime_ = 0.0;
        updateTime_ = 0.0;
    }


    /*!
     * \brief Called after the simulation has been finished
     *        sucessfully.
     */
    void finalize()
    {
        if (gridView().comm().rank() == 0) {
            Scalar totalTime = std::max(1e-100, assembleTime_ + solveTime_ + updateTime_);
            int numCores = this->gridView().comm().size();
            std::cout << "Simulation of problem '" << asImp_().name() << "' finished.\n"
                      << "Timing receipt [s] (solve total/assemble/linear solve/update): "
                      << totalTime  << " (" << totalTime*numCores << " cummulative, " << numCores <<" processes) / "
                      << assembleTime_  << " (" << assembleTime_/totalTime*100 << "%) / "
                      << solveTime_ << " (" << solveTime_/totalTime*100 << "%) / "
                      << updateTime_ << " (" << updateTime_/totalTime*100 << "%)"
                      << "\n";
        }
    }

    /*!
     * \brief Returns the total wall time spend on solving the
     *        system [s].
     */
    Scalar solveTime() const
    { return solveTime_; }

    /*!
     * \brief Returns the total wall time spend on updating the
     *        iterative solutions [s].
     */
    Scalar updateTime() const
    { return updateTime_; }

    /*!
     * \brief Evaluate the boundary conditions for a boundary segment.
     *
     * \param values Stores the fluxes over the boundary segment.
     * \param context The object representing the execution context from which this method is called.
     * \param spaceIdx The local index of the spatial entity which represents the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void boundary(BoundaryRateVector &values,
                  const Context &context,
                  int spaceIdx, int timeIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Problem does not provide a boundary() method"); }

    /*!
     * \brief Evaluate the constraints for a control volume.
     *
     * \param constraints Stores the values of the primary variables at a given spatial and temporal location.
     * \param context The object representing the execution context from which this method is called.
     * \param spaceIdx The local index of the spatial entity which represents the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void constraints(Constraints &constraints,
                     const Context &context,
                     int spaceIdx, int timeIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Problem does not provide a constraints() method"); }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param rate Stores the values of the volumetric creation/anihilition rates of the conserved quantities.
     * \param context The object representing the execution context from which this method is called.
     * \param spaceIdx The local index of the spatial entity which represents the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void source(RateVector &rate,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Problem does not provide a source() method"); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values Stores the primary variables.
     * \param context The object representing the execution context from which this method is called.
     * \param spaceIdx The local index of the spatial entity which represents the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    void initial(PrimaryVariables &values,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    { DUNE_THROW(Dune::InvalidStateException, "Problem does not provide a initial() method"); }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     *
     * \param context The object representing the execution context from which this method is called.
     * \param spaceIdx The local index of the spatial entity which represents the boundary segment.
     * \param timeIdx The index used for the time discretization
     */
    template <class Context>
    Scalar extrusionFactor(const Context &context,
                           int spaceIdx, int timeIdx) const
    { return asImp_().extrusionFactor(); }

    Scalar extrusionFactor() const
    { return 1.0; }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Called by the time manager before the time integration.
     */
    void preTimeStep()
    {}

    /*!
     * \brief Called by Ewoms::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
        const int maxFails = GET_PARAM(TypeTag, unsigned, MaxTimeStepDivisions);
        const Scalar minTimeStepSize = GET_PARAM(TypeTag, Scalar, MinTimeStepSize);

        // if the time step size of the time manager is smaller than
        // the specified minimum size and we're not going to finish
        // the simulation or an episode, try with the minimum size.
        if (timeManager().timeStepSize() < minTimeStepSize &&
            !timeManager().episodeWillBeOver() &&
            !timeManager().willBeFinished())
        {
            timeManager().setTimeStepSize(minTimeStepSize);
        }

        for (int i = 0; i < maxFails; ++i) {
            if (model_.update(newtonMethod_)) {
                assembleTime_ += newtonMethod_.assembleTime();
                solveTime_ += newtonMethod_.solveTime();
                updateTime_ += newtonMethod_.updateTime();

                return;
            }

            assembleTime_ += newtonMethod_.assembleTime();
            solveTime_ += newtonMethod_.solveTime();
            updateTime_ += newtonMethod_.updateTime();

            Scalar dt = timeManager().timeStepSize();
            Scalar nextDt = dt / 2;
            if (nextDt < minTimeStepSize)
                break; // give up: we can't make the time step smaller anymore!
            timeManager().setTimeStepSize(nextDt);

            // update failed
            if (gridView().comm().rank() == 0)
                std::cout << "Newton solver did not converge with "
                          << "dt=" << dt << " seconds. Retrying with time step of "
                          << nextDt << " seconds\n";
        }

        DUNE_THROW(Dune::MathError,
                   "Newton solver didn't converge after "
                   << maxFails
                   << " time-step divisions. dt="
                   << timeManager().timeStepSize());
    }

    /*!
     * \brief Returns the newton method object
     */
    NewtonMethod &newtonMethod()
    { return newtonMethod_; }

    /*!
     * \copydoc newtonMethod()
     */
    const NewtonMethod &newtonMethod() const
    { return newtonMethod_; }

    /*!
     * \brief Called by Ewoms::TimeManager whenever a solution for a
     *        time step has been computed and the simulation time has
     *        been updated.
     */
    Scalar nextTimeStepSize()
    {
        return std::min(GET_PARAM(TypeTag, Scalar, MaxTimeStepSize),
                        newtonMethod_.suggestTimeStepSize(timeManager().timeStepSize()));
    }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behavior is to write one restart file every 10 time
     * steps. This file is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteRestartFile() const
    {
        return timeManager().timeStepIndex() > 0 &&
            (timeManager().timeStepIndex() % 10 == 0);
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behavior is to write out every the solution for
     * very time step. This file is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteOutput() const
    { return true; }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    { }

    /*!
     * \brief Called by the time manager after everything which can be
     *        done about the current time step is finished and the
     *        model should be prepared to do the next time integration.
     */
    void advanceTimeLevel()
    { model_.advanceTimeLevel(); }

    /*!
     * \brief Called when the end of an simulation episode is reached.
     *
     * Typically a new episode should be started in this method.
     */
    void episodeEnd()
    {
        std::cerr << "The end of an episode is reached, but the problem "
                  << "does not override the episodeEnd() method. "
                  << "Doing nothing!\n";
    }
    // \}

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     * It could be either overwritten by the problem files, or simply
     * declared over the setName() function in the application file.
     */
    const char *name() const
    { return simName_.c_str(); }

    /*!
     * \brief Set the problem name.
     *
     * This static method sets the simulation name, which should be
     * called before the application problem is declared! If not, the
     * default name "sim" will be used.
     *
     * \param newName The problem's name
     */
    void setName(const char *newName)
    { simName_ = newName; }

    /*!
     * \brief The GridView which used by the problem.
     */
    const GridView &gridView() const
    { return gridView_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the smallest values.
     */
    const GlobalPosition &bboxMin() const
    { return bboxMin_; }

    /*!
     * \brief The coordinate of the corner of the GridView's bounding
     *        box with the largest values.
     */
    const GlobalPosition &bboxMax() const
    { return bboxMax_; }

    /*!
     * \brief Returns the mapper for vertices to indices.
     */
    const VertexMapper &vertexMapper() const
    { return vertexMapper_; }

    /*!
     * \brief Returns the mapper for elements to indices.
     */
    const ElementMapper &elementMapper() const
    { return elementMapper_; }

    /*!
     * \brief Returns TimeManager object used by the simulation
     */
    TimeManager &timeManager()
    { return *timeManager_; }

    /*!
     * \copydoc timeManager()
     */
    const TimeManager &timeManager() const
    { return *timeManager_; }

    /*!
     * \brief Returns numerical model used for the problem.
     */
    Model &model()
    { return model_; }

    /*!
     * \copydoc model()
     */
    const Model &model() const
    { return model_; }
    // \}

    /*!
     * \name Restart mechanism
     */
    // \{

    /*!
     * \brief This method writes the complete state of the simulation
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.ers</tt>. (Ewoms ReStart
     * file.)  See Ewoms::Restart for details.
     */
    void serialize()
    {
        typedef Ewoms::Restart Restarter;
        Restarter res;
        res.serializeBegin(asImp_());
        if (gridView().comm().rank() == 0)
            std::cout << "Serialize to file '" << res.fileName() << "'"
                      << ", time step size: " << asImp_().timeManager().timeStepSize() << "\n";

        timeManager().serialize(res);
        asImp_().serialize(res);
        res.serializeEnd();
    }

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.ers</tt>. (Ewoms ReStart
     * file.)  See Ewoms::Restart for details.
     *
     * \tparam Restarter The serializer type
     *
     * \param res The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        createResultWriter_();
        resultWriter_->serialize(res);
        model().serialize(res);
    }

    /*!
     * \brief Load a previously saved state of the whole simulation
     *        from disk.
     *
     * \param tRestart The simulation time on which the program was
     *                 written to disk.
     */
    void restart(Scalar tRestart)
    {
        typedef Ewoms::Restart Restarter;

        Restarter res;

        res.deserializeBegin(asImp_(), tRestart);
        if (gridView().comm().rank() == 0)
            std::cout << "Deserialize from file '" << res.fileName() << "'\n";
        timeManager().deserialize(res);
        asImp_().deserialize(res);
        res.deserializeEnd();
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        createResultWriter_();
        resultWriter_->deserialize(res);
        model().deserialize(res);
    }

    // \}

    /*!
     * \brief Write the relevant secondary variables of the current
     *        solution into an VTK output file.
     *
     * \param verbose If true, then a message will be printed to stdout if a file is written
     */
    void writeOutput(bool verbose = true)
    {
        // write the current result to disk
        if (asImp_().shouldWriteOutput()) {
            if (verbose && gridView().comm().rank() == 0)
                std::cout << "Writing result file for \"" << asImp_().name() << "\"\n";

            // calculate the time _after_ the time was updated
            Scalar t = timeManager().time() + timeManager().timeStepSize();
            createResultWriter_();
            resultWriter_->beginWrite(t);
            model().addOutputVtkFields(*resultWriter_);
            resultWriter_->endWrite();
        }
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void init_()
    {
        // calculate the bounding box of the local partition of the grid view
        VertexIterator vIt = gridView_.template begin<dim>();
        const VertexIterator vEndIt = gridView_.template end<dim>();
        for (; vIt!=vEndIt; ++vIt) {
            for (int i=0; i<dim; i++) {
                bboxMin_[i] = std::min(bboxMin_[i], vIt->geometry().corner(0)[i]);
                bboxMax_[i] = std::max(bboxMax_[i], vIt->geometry().corner(0)[i]);
            }
        }

        // communicate to get the bounding box of the whole domain
        for (int i = 0; i < dim; ++i) {
            bboxMin_[i] = gridView_.comm().min(bboxMin_[i]);
            bboxMax_[i] = gridView_.comm().max(bboxMax_[i]);
        }

        // set a default name for the problem
        simName_ = "sim";

        resultWriter_ = NULL;
    }

    // makes sure that the result writer exists
    void createResultWriter_()
    { if (!resultWriter_) resultWriter_ = new VtkMultiWriter(gridView_, asImp_().name()); }

    //! Returns the applied VTK-writer for the output
    VtkMultiWriter& resultWriter()
    {
        createResultWriter_();
        return *resultWriter_;
    }

    // CPU time keeping
    Scalar assembleTime_;
    Scalar solveTime_;
    Scalar updateTime_;

    std::string simName_;
    const GridView gridView_;

    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;

    ElementMapper elementMapper_;
    VertexMapper vertexMapper_;

    TimeManager *timeManager_;

    Model model_;

    NewtonMethod newtonMethod_;

    VtkMultiWriter *resultWriter_;
};

}

#endif
