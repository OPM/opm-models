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
 * \brief Provides a mechanism to dispatch work to separate threads
 */
#ifndef EWOMS_TASKLETS_HH
#define EWOMS_TASKLETS_HH

#include <stdexcept>
#include <cassert>
#include <thread>
#include <mutex>
#include <queue>

namespace Ewoms {

/*!
 * \brief The base class for tasklets.
 *
 * Tasklets are a generic mechanism for potentially running work in a separate thread.
 */
class TaskletInterface
{
public:
    TaskletInterface() {}
    virtual ~TaskletInterface() {}
    virtual void run() = 0;
    virtual bool isEndMarker () const { return false; }
};

/*!
 * \brief Handles where a given tasklet is run.
 *
 * Depending on the "runAsync" constructor parameter, a tasklet can either be a separate
 * worker thread or it can be the main thread.
 */
class TaskletRunner
{
    /// \brief Implements a barrier. This class can only be used in the asynchronous case.
    class BarrierTasklet : public TaskletInterface
    {
    public:
        BarrierTasklet()
        { barrierMutex.lock(); }

        void run()
        { barrierMutex.unlock(); }

        std::mutex barrierMutex;
    };

    /// \brief TerminateThreadTasklet class
    /// Empty tasklet marking thread termination.
    class TerminateThreadTasklet : public TaskletInterface
    {
    public:
        void run()
        { }

        bool isEndMarker() const
        { return true; }
    };

public:
    // prohibit copying of tasklet runners
    TaskletRunner(const TaskletRunner&) = delete;

    TaskletRunner(const bool runAsync)
    {
        if (runAsync) {
            // make sure that the runner thread blocks when the tasklet queue is empty
            runnerMutex_.lock();

            // create a worker thread
            thread_.reset(new std::thread(startThread_, this));
        }
    }

    /*!
     * \brief Destructor
     *
     * If a worker thread was created to run the tasklets, this method waits until the
     * worker thread has been terminated, i.e. all scheduled tasklets are guaranteed to
     * be completed.
     */
    ~TaskletRunner()
    {
        if (thread_) {
            // dispatch a tasklet which will terminate the worker thread
            dispatch(std::make_shared<TerminateThreadTasklet>());

            // wait until the thread has been terminated
            thread_->join();
        }
    }

    /*!
     * \brief Add a new tasklet.
     *
     * The tasklet is either run immediately or deferred to a separate thread.
     */
    void dispatch(std::shared_ptr<TaskletInterface> tasklet)
    {
        if (!thread_)
            // run the tasklet immediately in synchronous mode.
            tasklet->run();
        else {
            // lock mutex for the tasklet queue to make sure that nobody messes with the
            // task queue
            taskletQueueMutex_.lock();

            // add the tasklet to the queue
            taskletQueue_.push(tasklet);
            // fire up the worker thread
            runnerMutex_.unlock();

            taskletQueueMutex_.unlock();
        }
    }

    /*!
     * \brief Make sure that all tasklets have been completed after this method has been called
     */
    void barrier()
    {
        if (!thread_)
            // nothing needs to be done to implement a barrier in synchronous mode
            return;

        // dispatch a barrier tasklet and wait until it has been run by the worker thread
        auto barrierTasklet = std::make_shared<BarrierTasklet>();
        dispatch(barrierTasklet);
        barrierTasklet->barrierMutex.lock();
    }

protected:
    // main function of the worker thread
    static void startThread_(TaskletRunner* taskletRunner)
    { taskletRunner->run_(); }

    //! do the work until the queue received an end tasklet
    void run_()
    {
        while (true) {
            // wait until tasklets have been pushed to the queue.
            //
            // The unlocking is done in the thread that adds a task
            runnerMutex_.lock();

            // lock mutex for access to taskletQueue_
            taskletQueueMutex_.lock();

            // remove tasklet from queue
            std::shared_ptr<TaskletInterface> tasklet = taskletQueue_.front();
            taskletQueue_.pop();

            // if the queue is not yet empty, make sure that we are will process the next
            // tasklet immediately after we're finished with the current one.
            if (!taskletQueue_.empty())
                runnerMutex_.unlock();

            // unlock mutex for access to taskletQueue_
            taskletQueueMutex_.unlock();

            // if tasklet is an end marker, terminate the thread
            if (tasklet->isEndMarker()) {
                if(!taskletQueue_.empty())
                    throw std::logic_error("TaskletRunner: Not all queued tasklets were executed");
                return;
            }

            // execute tasklet action
            tasklet->run();
        }
    }

    std::unique_ptr<std::thread> thread_;
    std::queue<std::shared_ptr<TaskletInterface> > taskletQueue_;
    std::mutex taskletQueueMutex_;
    std::mutex runnerMutex_;
};

} // end namespace Opm
#endif
