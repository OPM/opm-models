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
 * \copydoc Ewoms::ThreadManager
 */
#ifndef EWOMS_THREAD_MANAGER_HH
#define EWOMS_THREAD_MANAGER_HH

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ewoms/parallel/locks.hh>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/common/propertysystem.hh>

#include <dune/common/version.hh>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(ThreadsPerProcess);
}

/*!
 * \brief Simplifies multi-threaded capabilities.
 */
template <class TypeTag>
class ThreadManager
{
public:
    enum {
#if defined(_OPENMP) || DOXYGEN
        //! Specify whether OpenMP is really available or not
        isFake = false
#else
        isFake = true
#endif
    };

    /*!
     * \brief Register all run-time parameters of the thread manager.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, ThreadsPerProcess,
                             "The maximum number of threads to be instantiated per process "
                             "('-1' means 'automatic')");
    }

    static void init()
    {
        numThreads_ = EWOMS_GET_PARAM(TypeTag, int, ThreadsPerProcess);

        // some safety checks. This is pretty ugly macro-magic, but so what?
#if !defined(_OPENMP)
        if (numThreads_ != 1 && numThreads_ != -1)
            OPM_THROW(std::invalid_argument,
                      "OpenMP is not available. The only valid values for "
                      "threads-per-process is 1 and -1 but it is " << numThreads_ << "!");
        numThreads_ = 1;
#elif !DUNE_VERSION_NEWER(DUNE_COMMON, 2,4) && !defined NDEBUG
        if (numThreads_ != 1)
            OPM_THROW(std::invalid_argument,
                      "You seem to be using Dune "
                      <<DUNE_COMMON_VERSION_MAJOR<<"."<<DUNE_COMMON_VERSION_MINOR
                      << " in debug mode and with the number of OpenMP threads larger than 1. "
                      "This Dune version does not support thread parallelism in debug mode!");
        numThreads_ = 1;

#elif DUNE_VERSION_NEWER(DUNE_COMMON, 2,4) && !defined NDEBUG && defined DUNE_INTERFACECHECK
        if (numThreads_ != 1)
            OPM_THROW(std::invalid_argument,
                      "You explicitly enabled Barton-Nackman interface checking in Dune. "
                      "The Dune implementation of this is currently incompatible with "
                      "thread parallelism!");
        numThreads_ = 1;
#else
        if (numThreads_ == 0)
            OPM_THROW(std::invalid_argument,
                      "Zero threads per process are not possible: It must be at least 1, "
                      "(or -1 for 'automatic')!");
#endif

#ifdef _OPENMP
        // actually limit the number of threads and get the number of threads which are
        // used in the end.
        if (numThreads_ > 0)
            omp_set_num_threads(numThreads_);

        numThreads_ = omp_get_max_threads();
#endif
    }

    /*!
     * \brief Return the maximum number of threads of the current process.
     */
    static int maxThreads()
    { return numThreads_; }

    /*!
     * \brief Return the index of the current OpenMP thread
     */
    static int threadId()
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

private:
    static int numThreads_;
};

template <class TypeTag>
int ThreadManager<TypeTag>::numThreads_ = 1;
} // namespace Ewoms

#endif
