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

#include <opm/core/utility/PropertySystem.hpp>

#include <dune/common/timer.hh>

#include <iostream>
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
        verbose_ = verbose && Dune::MPIHelper::getCollectiveCommunication().rank() == 0;

        // deal with the restart stuff
        Scalar restartTime = EWOMS_GET_PARAM(TypeTag, Scalar, RestartTime);
        if (restartTime > -1e100) {
            time_ = restartTime;
            doRestart_ = true;
        }
        else {
            time_ = 0.0;
            doRestart_ = false;
        }

        timeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
        endTime_ = EWOMS_GET_PARAM(TypeTag, Scalar, EndTime);

        episodeIdx_ = 0;
        episodeLength_ = 1e100;
        episodeStartTime_ = 0;

        finished_ = false;

        if (verbose_)
            std::cout << "Allocating the grid\n"
                      << std::flush;

        gridManager_.reset(new GridManager(*this));

        if (verbose_)
            std::cout << "Distributing the grid\n"
                      << std::flush;
        gridManager_->loadBalance();

        if (verbose_)
            std::cout << "Allocating the problem and the model\n"
                      << std::flush;
        model_.reset(new Model(*this));
        problem_.reset(new Problem(*this));
    }

    /*!
     * \brief Registers all runtime parameters used by the simulation.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EndTime,
                             "The time at which the simulation is finished [s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitialTimeStepSize,
                             "The size of the initial time step [s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, RestartTime,
                             "The time time at which a simulation restart should "
                             "be attempted [s]");

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
     * \brief Return the time \f$\mathrm{[s]}\f$ before the time integration.
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
     * \brief Returns the current wall time (cpu time).
     */
    double wallTime() const
    { return timer_.elapsed(); }

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
        return finished_
            || (this->time() + std::max(std::abs(this->time()), timeStepSize())*1e-8 >= endTime());
    }

    /*!
     * \brief Returns true if the simulation is finished after the
     *        time level is incremented by the current time step size.
     */
    bool willBeFinished() const
    {
        return
            finished_ ||
            (this->time() + std::max(std::abs(this->time()), timeStepSize())*1e-8 + timeStepSize_
             >= endTime());
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
        episodeStartTime_ = time_;
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
    { return this->time() < episodeStartTime_ + episodeLength()*1e-8; }

    /*!
     * \brief Returns true if the current episode is finished at the
     *        current time.
     */
    bool episodeIsOver() const
    { return this->time() > episodeStartTime_ + episodeLength()*(1 - 1e-8); }

    /*!
     * \brief Returns true if the current episode will be finished
     *        after the current time step.
     */
    bool episodeWillBeOver() const
    {
        return
            this->time() + timeStepSize()
            >= episodeStartTime_ + episodeLength() * (1 - 1e-8);
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
        return std::max<Scalar>(0.0,
                                episodeLength()
                                - (this->time() - episodeStartTime()));
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
        init_();

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
            Scalar dt = timeStepSize();

            // post-process the current solution
            problem_->endTimeStep();

            // write the result to disk
            if (problem_->shouldWriteOutput())
                problem_->writeOutput();

            // prepare the model for the next time integration
            problem_->advanceTimeLevel();

            // advance the simulated time by the current time step size
            time_ += dt;
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
            else {
                // notify the problem that the timestep is done and ask it
                // for a suggestion for the next timestep size
                // set the time step size for the next step
                setTimeStepSize(problem_->nextTimeStepSize());
            }

            // write restart file if mandated by the problem
            if (problem_->shouldWriteRestartFile())
                serialize();

            if (verbose_) {
                std::cout << "Time step " << timeStepIndex() << " done. "
                          << "Wall time:" << timer_.elapsed()
                          << ", time:" << this->time() << ", time step size:" << timeStepSize()
                          << "\n" << std::flush;
            }
        }

        problem_->finalize();
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
     * \brief Load a previously saved state of the whole simulation
     *        from disk.
     *
     * \param restartTime The simulation time on which the program was
     *                    written to disk.
     */
    void restart(Scalar restartTime)
    {
        typedef Ewoms::Restart Restarter;

        Restarter res;

        res.deserializeBegin(*this, restartTime);
        if (gridView().comm().rank() == 0)
            std::cout << "Deserialize from file '" << res.fileName() << "'\n" << std::flush;
        this->deserialize(res);
        problem_->deserialize(res);
        model_->deserialize(res);
        res.deserializeEnd();
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
    /*!
     * \brief Apply the initial condition and write it to disk.
     *
     * This method also deserializes a previously saved state of the
     * simulation from disk if instructed to do so.
     */
    void init_()
    {
        timer_.reset();

        // apply the initial solution
        if (verbose_)
            std::cout << "Applying the initial solution of the \"" << problem_->name() << "\" problem\n"
                      << std::flush;
        problem_->init();

        // try restart the simulation if requested
        if (doRestart_) {
            if (verbose_)
                std::cout << "Loading the previously saved state for t=" << time_ << "\n"
                          << std::flush;

            restart(time_);
        }
        else {
            // write initial condition (if problem is not restarted)
            time_ -= timeStepSize_;
            if (problem_->shouldWriteOutput())
                problem_->writeOutput();
            time_ += timeStepSize_;
        }
    }

    std::unique_ptr<GridManager> gridManager_;
    std::unique_ptr<Model> model_;
    std::unique_ptr<Problem> problem_;

    int episodeIdx_;
    Scalar episodeStartTime_;
    Scalar episodeLength_;

    Dune::Timer timer_;
    Scalar startTime_;
    Scalar time_;
    Scalar endTime_;

    Scalar timeStepSize_;
    int timeStepIdx_;
    bool finished_;
    bool verbose_;
    bool doRestart_;
};
} // namespace Ewoms

#endif
