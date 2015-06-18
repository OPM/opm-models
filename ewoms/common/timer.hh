// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2014 by Andreas Lauser

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
 */
#ifndef EWOMS_TIMER_HH
#define EWOMS_TIMER_HH

#include <time.h>

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Ewoms {
/*
 * \brief Provides an encapsulation to measure the system time
 *
 * This means the wall clock time used by the simulation, the CPU time
 * used by all threads of a single process and the CPU time used by
 * the overall simulation. (i.e., the time used by all threads of all
 * involved processes.)
 */
class Timer
{
    struct TimeData
    {
        struct timespec realtimeData;
        struct timespec cputimeData;
    };
public:
    Timer()
    { halt(); }

    /*!
     * \brief Start counting the time resources used by the simulation.
     */
    void start()
    {
        isStopped_ = false;
        measure_(startTime_);
    }

    /*!
     * \brief Stop counting the time resources used by the simulation.
     */
    void stop()
    {
        isStopped_ = true;
        measure_(stopTime_);
    }

    /*!
     * \brief Stop the measurement and always return 0 for all timing values
     */
    void halt()
    {
        isStopped_ = true;

        measure_(startTime_);
        stopTime_ = startTime_;
    }

    /*!
     * \brief Return the real time [s] elapsed
     *
     * If stop() was not yet called, this returns the time elapsed
     * since the last call to start().
     */
    double realTimeElapsed() const
    {
        TimeData stopTime(stopTime_);

        if (!isStopped_)
            measure_(stopTime);

        const auto &t1 = startTime_.realtimeData;
        const auto &t2 = stopTime.realtimeData;

        return
            static_cast<double>(t2.tv_sec - t1.tv_sec)
            + static_cast<double>(t2.tv_nsec - t1.tv_nsec)/1e9;
    }

    /*!
     * \brief Return the CPU time [s] used by all threads of the local process
     *
     * If stop() was not yet called, this returns the time elapsed
     * since the last call to start().
     */
    double cpuTimeElapsed() const
    {
        TimeData stopTime(stopTime_);

        if (!isStopped_)
            measure_(stopTime);

        const auto &t1 = startTime_.cputimeData;
        const auto &t2 = stopTime.cputimeData;

        return
            static_cast<double>(t2.tv_sec - t1.tv_sec)
            + static_cast<double>(t2.tv_nsec - t1.tv_nsec)/1e9;
    }

    /*!
     * \brief Return the CPU time [s] used by all threads of the all processes of the simulation
     *
     * If stop() was not yet called, this returns the time elapsed
     * since the last call to start(). Note that this method must be
     * called synchronously by all processes of the simulation...
     */
    double globalCpuTimeElapsed() const
    {
        TimeData stopTime(stopTime_);

        if (!isStopped_)
            measure_(stopTime);

        double val = cpuTimeElapsed();
        double globalVal = val;

#if HAVE_MPI
        MPI_Reduce(&val,
            &globalVal,
            /*count=*/1,
            MPI_DOUBLE,
            MPI_SUM,
            /*rootRank=*/0,
            MPI_COMM_WORLD);
#endif
        return globalVal;
    }

private:
    // measure the current time and put it into the object passed via
    // the argument.
    static void measure_(TimeData& timeData)
    {
        // This method is Linux (or at least POSIX) specific. TODO:
        // add implementations for different operating systems.

        // measure the real time
        clock_gettime(CLOCK_MONOTONIC, &timeData.realtimeData);

        // measure the CPU time
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timeData.cputimeData);
    }

    bool isStopped_;
    TimeData startTime_;
    TimeData stopTime_;
};
} // namespace Ewoms

#endif
