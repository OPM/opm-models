/*
  Copyright (C) 2009-2013 by Andreas Lauser
  Copyright (C) 2011-2012 by Markus Wolff
  Copyright (C) 2011 by Benjamin Faigle

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
*/
/*!
 * \file
 * \copydoc Ewoms::Simulator
 */
#ifndef EWOMS_SIMULATOR_HH
#define EWOMS_SIMULATOR_HH

#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/io/restart.hh>
#include <ewoms/common/parametersystem.hh>

#include <opm/core/utility/PropertySystem.hpp>

#include <ewoms/common/timer.hh>

#include <iostream>
#include <iomanip>
#include <memory>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GridManager);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Model);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(EndTime);
NEW_PROP_TAG(RestartTime);
NEW_PROP_TAG(InitialTimeStepSize);
}}

namespace Ewoms {
/*!
 * \ingroup Simulator
 *
 * \brief Manages the initializing and running of time dependent
 *        problems.
 *
 * This class instantiates the grid, the model and the problem to be
 * simlated and runs the simulation loop. The time axis is treated as
 * a sequence of "episodes" which are defined as time intervals for
 * which the problem exhibits boundary conditions and source terms
 * that do not depend on time.
 */
template <class TypeTag>
class Simulator
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridManager) GridManager;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;


public:
    // do not allow to copy simulators around
    Simulator(const Simulator &) = delete;

    Simulator(bool verbose = true)
    {
        setupTimer_.start();

        verbose_ = verbose && Ewoms::MpiHelper::getCollectiveCommunication().rank() == 0;

        timeStepIdx_ = 0;
        startTime_ = 0.0;
        time_ = 0.0;
        endTime_ = EWOMS_GET_PARAM(TypeTag, Scalar, EndTime);
        timeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);

        episodeIdx_ = 0;
        episodeStartTime_ = 0;
        episodeLength_ = 1e100;

        finished_ = false;

        if (verbose_)
            std::cout << "Allocating the grid\n" << std::flush;
        gridManager_.reset(new GridManager(*this));

        if (verbose_)
            std::cout << "Distributing the grid\n" << std::flush;
        gridManager_->loadBalance();

        if (verbose_)
            std::cout << "Allocating the model\n" << std::flush;
        model_.reset(new Model(*this));

        if (verbose_)
            std::cout << "Allocating the problem\n" << std::flush;
        problem_.reset(new Problem(*this));

        if (verbose_)
            std::cout << "Finish init of the model\n" << std::flush;
        model_->finishInit();

        if (verbose_)
            std::cout << "Finish init of the problem\n" << std::flush;
        problem_->finishInit();

        if (verbose_)
            std::cout << "Construction of simulation done\n" << std::flush;
    }

    /*!
     * \brief Registers all runtime parameters used by the simulation.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EndTime,
                             "The simulation time at which the simulation is finished [s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitialTimeStepSize,
                             "The size of the initial time step [s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, RestartTime,
                             "The simulation time at which a restart should be attempted [s]");

        GridManager::registerParameters();
        Model::registerParameters();
        Problem::registerParameters();
    }

    /*!
     * \brief Return a reference to the grid manager of simulation
     */
    GridManager &gridManager()
    { return *gridManager_; }

    /*!
     * \brief Return a reference to the grid manager of simulation
     */
    const GridManager &gridManager() const
    { return *gridManager_; }

    /*!
     * \brief Return the grid view for which the simulation is done
     */
    GridView &gridView()
    { return gridManager_->gridView(); }

    /*!
     * \brief Return the grid view for which the simulation is done
     */
    const GridView &gridView() const
    { return gridManager_->gridView(); }

    /*!
     * \brief Return the physical model used in the simulation
     */
    Model &model()
    { return *model_; }

    /*!
     * \brief Return the physical model used in the simulation
     */
    const Model &model() const
    { return *model_; }

    /*!
     * \brief Return the object which specifies the pysical setup of
     *        the simulation
     */
    Problem &problem()
    { return *problem_; }

    /*!
     * \brief Return the object which specifies the pysical setup of
     *        the simulation
     */
    const Problem &problem() const
    { return *problem_; }

    /*!
     * \brief Set the time of the start of the simulation.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     */
    void setStartTime(Scalar t)
    { startTime_ = t; }

    /*!
     * \brief Return the time of the start of the simulation.
     */
    Scalar startTime() const
    { return startTime_; }

    /*!
     * \brief Set the current simulated time, don't change the current
     *        time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     */
    void setTime(Scalar t)
    { time_ = t; }

    /*!
     * \brief Set the current simulated time and the time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     * \param stepIdx The new time step index
     */
    void setTime(Scalar t, int stepIdx)
    {
        time_ = t;
        timeStepIdx_ = stepIdx;
    }

    /*!
     * \brief Return the number of seconds of simulated time which have elapsed since the
     *        start time.
     *
     * To get the time after the time integration, you have to add
     * timeStepSize() to time().
     */
    Scalar time() const
    { return time_; }

    /*!
     * \brief Set the time of simulated seconds at which the simulation runs.
     *
     * \param t The time \f$\mathrm{[s]}\f$ at which the simulation is finished
     */
    void setEndTime(Scalar t)
    { endTime_ = t; }

    /*!
     * \brief Returns the number of (simulated) seconds which the simulation
     *        runs.
     */
    Scalar endTime() const
    { return endTime_; }

    /*!
     * \brief Returns the current wall time required by actually running the simulation.
     */
    const Ewoms::Timer& timer() const
    { return executionTimer_; }

    /*!
     * \brief Returns total wall clock time required to write the visualization files of
     *        all time steps to disk
     */
    Scalar totalWriteTime() const
    { return totalWriteTime_; }

    /*!
     * \brief Returns the wall time required by setting up and initializing the simulation.
     */
    double setupTime() const
    { return setupTimer_.realTimeElapsed(); }

    /*!
     * \brief Set the current time step size to a given value.
     *
     * If the step size would exceed the length of the current
     * episode, the timeStep() method will take care that the step
     * size won't exceed the episode or the end of the simulation,
     * though.
     *
     * \param timeStepSize The new value for the time step size \f$\mathrm{[s]}\f$
     */
    void setTimeStepSize(Scalar timeStepSize)
    { timeStepSize_ = timeStepSize; }

    /*!
     * \brief Returns the time step length \f$\mathrm{[s]}\f$ so that we
     *        don't miss the beginning of the next episode or cross
     *        the end of the simlation.
     */
    Scalar timeStepSize() const
    {
        Scalar maximumTimeStepSize =
            std::min(episodeMaxTimeStepSize(),
                     std::max<Scalar>(0.0, endTime() - this->time()));

        return std::min(timeStepSize_, maximumTimeStepSize);
    }

    /*!
     * \brief Returns number of time steps which have been
     *        executed since the beginning of the simulation.
     */
    int timeStepIndex() const
    { return timeStepIdx_; }

    /*!
     * \brief Specify whether the simulation is finished
     *
     * \param yesno If true the simulation is considered finished
     *              before the end time is reached, else it is only
     *              considered finished if the end time is reached.
     */
    void setFinished(bool yesno = true)
    { finished_ = yesno; }

    /*!
     * \brief Returns true if the simulation is finished.
     *
     * This is the case if either setFinished(true) has been called or
     * if the end time is reached.
     */
    bool finished() const
    {
        assert(timeStepSize_ >= 0.0);
        Scalar eps = std::max(std::abs(this->time()), timeStepSize())*1e-8;
        return finished_ || (this->time() + eps >= endTime());
    }

    /*!
     * \brief Returns true if the simulation is finished after the
     *        time level is incremented by the current time step size.
     */
    bool willBeFinished() const
    {
        return
            finished_ ||
            (this->time() + std::max(std::abs(this->time()), timeStepSize())*1e-8
             + timeStepSize_ >= endTime());
    }

    /*!
     * \brief Aligns the time step size to the episode boundary and to
     *        the end time of the simulation.
     */
    Scalar maxTimeStepSize() const
    {
        if (finished())
            return 0.0;

        return std::min(episodeMaxTimeStepSize(),
                        std::max<Scalar>(0.0, endTime() - this->time()));
    }

    /*!
     * \brief Change the current episode of the simulation.
     *
     * \param episodeStartTime Time when the episode began \f$\mathrm{[s]}\f$
     * \param episodeLength Length of the episode \f$\mathrm{[s]}\f$
     */
    void startNextEpisode(Scalar episodeStartTime, Scalar episodeLength)
    {
        ++episodeIdx_;
        episodeStartTime_ = episodeStartTime;
        episodeLength_ = episodeLength;
    }

    /*!
     * \brief Start the next episode, but don't change the episode
     *        identifier.
     *
     * \param len Length of the episode \f$\mathrm{[s]}\f$, infinite if not
     *            specified.
     */
    void startNextEpisode(Scalar len = 1e100)
    {
        ++episodeIdx_;
        episodeStartTime_ = startTime_ + time_;
        episodeLength_ = len;
    }

    /*!
     * \brief Sets the index of the current episode.
     *
     * Use this method with care!
     */
    void setEpisodeIndex(int episodeIdx)
    { episodeIdx_ = episodeIdx; }

    /*!
     * \brief Returns the index of the current episode.
     *
     * The first episode has the index 0.
     */
    int episodeIndex() const
    { return episodeIdx_; }

    /*!
     * \brief Returns the absolute time when the current episode
     *        started \f$\mathrm{[s]}\f$.
     */
    Scalar episodeStartTime() const
    { return episodeStartTime_; }

    /*!
     * \brief Returns the length of the current episode in
     *        simulated time \f$\mathrm{[s]}\f$.
     */
    Scalar episodeLength() const
    { return episodeLength_; }

    /*!
     * \brief Returns true if the current episode has just been started
     */
    bool episodeBegins() const
    { return startTime() + this->time() < episodeStartTime_ + episodeLength()*1e-8; }

    /*!
     * \brief Returns true if the current episode is finished at the
     *        current time.
     */
    bool episodeIsOver() const
    { return startTime() + this->time() > episodeStartTime_ + episodeLength()*(1 - 1e-8); }

    /*!
     * \brief Returns true if the current episode will be finished
     *        after the current time step.
     */
    bool episodeWillBeOver() const
    {
        return
            startTime()
            + this->time()
            + timeStepSize()
            >=
            episodeStartTime_
            + episodeLength()*(1 - 1e-8);
    }

    /*!
     * \brief Aligns the time step size to the episode boundary if the
     *        current time step crosses the boundary of the current episode.
     */
    Scalar episodeMaxTimeStepSize() const
    {
        // if the current episode is over and the simulation
        // wants to give it some extra time, we will return
        // the time step size it suggested instead of trying
        // to align it to the end of the episode.
        if (episodeIsOver())
            return 0.0;

        // make sure that we don't exceed the end of the
        // current episode.
        return std::max<Scalar>(0.0, episodeLength() - (this->time() - episodeStartTime()));
    }

    /*
     * \}
     */

    /*!
     * \brief Runs the simulation using a given problem class.
     *
     * This method makes sure that time steps sizes are aligned to
     * episode boundaries, amongst other stuff.
     */
    void run()
    {
        Scalar restartTime = EWOMS_GET_PARAM(TypeTag, Scalar, RestartTime);
        if (restartTime > -1e100) {
            // try to restart a previous simulation
            time_ = restartTime;

            Ewoms::Restart res;
            res.deserializeBegin(*this, time_);
            if (verbose_)
                std::cout << "Deserialize from file '" << res.fileName() << "'\n" << std::flush;
            this->deserialize(res);
            problem_->deserialize(res);
            model_->deserialize(res);
            res.deserializeEnd();
        }
        else {
            // if no restart is done, apply the initial solution
            if (verbose_)
                std::cout << "Applying the initial solution of the \"" << problem_->name()
                          << "\" problem\n" << std::flush;

            Scalar oldTimeStepSize = timeStepSize_;
            int oldTimeStepIdx = timeStepIdx_;
            timeStepSize_ = 0.0;
            timeStepIdx_ = -1;

            model_->applyInitialSolution();

            // write initial condition
            if (problem_->shouldWriteOutput())
                problem_->writeOutput();

            timeStepSize_ = oldTimeStepSize;
            timeStepIdx_ = oldTimeStepIdx;
        }

        setupTimer_.stop();

        totalWriteTime_ = 0.0;
        executionTimer_.start();

        // do the time steps
        while (!finished()) {
            if (episodeBegins())
                // notify the problem that a new episode has just been
                // started.
                problem_->beginEpisode();

            // pre-process the current solution
            problem_->beginTimeStep();

            // execute the time integration scheme
            problem_->timeIntegration();

            // post-process the current solution
            problem_->endTimeStep();

            // write the result to disk
            writeTimer_.start();
            if (problem_->shouldWriteOutput())
                problem_->writeOutput();
            writeTimer_.stop();
            totalWriteTime_ += writeTimer_.realTimeElapsed();

            // do the next time integration
            Scalar oldDt = timeStepSize();
            problem_->advanceTimeLevel();

            // advance the simulated time by the current time step size
            time_ += oldDt;
            ++timeStepIdx_;

            // notify the problem if an episode is finished
            if (episodeIsOver()) {
                // make the linearization not usable for the first
                // iteration of the next time step at the end of each
                // episode. Strictly speaking, this is a layering
                // violation as the simulator is not supposed to know
                // any specifics of the non-linear solver, but this
                // call is too easy to forget within the problem's
                // endEpisode(), so we do it here anyway!
                model_->jacobianAssembler().setLinearizationReusable(false);

                // Notify the problem about the end of the current
                // episode. This needs to define what to the next
                // episode or an exception is thrown...
                problem_->endEpisode();
            }
            else
                // ask the problem to provide the next time step size
                setTimeStepSize(problem_->nextTimeStepSize());

            // write restart file if mandated by the problem
            if (problem_->shouldWriteRestartFile())
                serialize();

            if (verbose_) {
                std::cout << "Time step " << timeStepIndex() << " done. "
                          << "Execution time: " << executionTimer_.realTimeElapsed() << " seconds" << humanReadableTime(executionTimer_.realTimeElapsed())
                          << ", simulation time: " << this->time() << " seconds" << humanReadableTime(this->time())
                          << ", time step size: " << oldDt << " seconds" << humanReadableTime(oldDt)
                          << ", next time step size: " << timeStepSize() << " seconds" << humanReadableTime(timeStepSize())
                          << "\n" << std::flush;
            }
        }

        executionTimer_.stop();

        problem_->finalize();
    }

    /*!
     * \brief Given a time step size in seconds, return it in a format which is more
     *        easily parsable by humans.
     *
     * e.g. 874000.0 will become "10.12 days"
     */
    static std::string humanReadableTime(Scalar timeInSeconds, bool isAmendment=true)
    {
        std::ostringstream oss;
        oss << std::setprecision(2);
        oss << std::fixed;
        if (isAmendment)
            oss << " (";
        if (timeInSeconds >= 365.25*24*60*60)
            oss << timeInSeconds/(365.25*24*60*60) << " years";
        else if (timeInSeconds >= 24.0*60*60)
            oss << timeInSeconds/(24.0*60*60) << " days";
        else if (timeInSeconds >= 60.0*60)
            oss << timeInSeconds/(60.0*60) << " hours";
        else if (timeInSeconds >= 60.0)
            oss << timeInSeconds/(60.0) << " minutes";
        else if (timeInSeconds < 1e-3)
            oss << timeInSeconds/(1e-6) << " microseconds";
        else if (timeInSeconds <= 1.0)
            oss << timeInSeconds/(1e-3) << " milliseconds";
        else if (!isAmendment)
            oss << timeInSeconds << " seconds";
        else
            return "";
        if (isAmendment)
            oss << ")";

        return oss.str();
    }

    /*!
     * \name Saving/restoring the simulation state
     * \{
     */

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
        res.serializeBegin(*this);
        if (gridView().comm().rank() == 0)
            std::cout << "Serialize to file '" << res.fileName() << "'"
                      << ", next time step size: " << timeStepSize()
                      << "\n" << std::flush;

        this->serialize(res);
        problem_->serialize(res);
        model_->serialize(res);
        res.serializeEnd();
    }

    /*!
     * \brief Write the time manager's state to a restart file.
     *
     * \tparam Restarter The type of the object which takes care to serialize
     *                   data
     * \param restarter The serializer object
     */
    template <class Restarter>
    void serialize(Restarter &restarter)
    {
        restarter.serializeSectionBegin("Simulator");
        restarter.serializeStream()
            << episodeIdx_ << " "
            << episodeStartTime_ << " "
            << episodeLength_ << " "
            << startTime_ << " "
            << time_ << " "
            << timeStepIdx_ << " ";
        restarter.serializeSectionEnd();
    }

    /*!
     * \brief Read the time manager's state from a restart file.
     *
     * \tparam Restarter The type of the object which takes care to deserialize
     *                   data
     * \param restarter The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter &restarter)
    {
        restarter.deserializeSectionBegin("Simulator");
        restarter.deserializeStream()
            >> episodeIdx_
            >> episodeStartTime_
            >> episodeLength_
            >> startTime_
            >> time_
            >> timeStepIdx_;
        restarter.deserializeSectionEnd();
    }

private:
    std::unique_ptr<GridManager> gridManager_;
    std::unique_ptr<Model> model_;
    std::unique_ptr<Problem> problem_;

    int episodeIdx_;
    Scalar episodeStartTime_;
    Scalar episodeLength_;

    Ewoms::Timer setupTimer_;
    Ewoms::Timer executionTimer_;
    Ewoms::Timer writeTimer_;
    Scalar totalWriteTime_;
    Scalar startTime_;
    Scalar time_;
    Scalar endTime_;

    Scalar timeStepSize_;
    int timeStepIdx_;
    bool finished_;
    bool verbose_;
};
} // namespace Ewoms

#endif
