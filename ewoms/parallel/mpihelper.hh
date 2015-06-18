// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
 * \copydoc Ewoms::MpiHelper
 */
#ifndef EWOMS_MPI_HELPER_HH
#define EWOMS_MPI_HELPER_HH

#include <dune/common/deprecated.hh>

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/ErrorMacros.hpp>

// this is a bit hacky, but required for DUNE compatibility
#ifdef DUNE_MPIHELPER
#error "You must include this file _before_ DUNE's mpihelper.hh!"
#endif
#define DUNE_MPIHELPER

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#else
#include <dune/common/collectivecommunication.hh>
#include <dune/common/mpicollectivecommunication.hh>
#endif

#if HAVE_MPI
#include <mpi.h>
#endif

namespace Ewoms {

/*!
 * \brief Initializes MPI at the beginning of the simulation and
 *        cleans up afterwards.
 *
 * This is basically the same idea as Dune's MPIHelper class, but
 * without the weird singleton API which might lead to race conditions
 * for some MPI libraries.  This class is intended to be used like
 * this:
 * \code
 * int main(int argc, char** argv)
 * {
 *   Ewoms::MpiHelper mpiHelper(argc, argv);
 *   // program code
 *   ...
 * }
 * \endcode
 *
 * Except that it requires to be instanced in the main function, it
 * is designed to be 100% compatible with Dune's MPIHelper.
 */
class MpiHelper
{
public:
    enum {
#if HAVE_MPI
        //! Specify whether MPI is really available or not
        isFake = false
#else
        //! Specify whether MPI is really available or not
        isFake = true
#endif
    };

/// The type of the MPI communicator object
#if HAVE_MPI
    typedef MPI_Comm MPICommunicator;
#else
    typedef Dune::No_Comm MPICommunicator;
#endif

    MpiHelper(int &argc, char **&argv)
    {
        // make it only possible to instanciate this class once during
        // the lifetime of the program
        static bool wasInstanciated = false;
        if (wasInstanciated)
            OPM_THROW(std::logic_error,
                      "Ewoms::MpiHelper may only be instanciated once!");
        wasInstanciated = true;

#if HAVE_MPI
        if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
            OPM_THROW(std::logic_error, "Initialization of MPI failed!");

        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
#endif
    }

    ~MpiHelper()
    {
#if HAVE_MPI
        MPI_Finalize();
#endif
    }

    /*!
     * \brief The deprecated singleton DUNE interface
     */
    DUNE_DEPRECATED_MSG(
        "Construct an instance of this class in your main function instead!")
    static MpiHelper &instance(int &argc, char **&argv)
    {
        static MpiHelper singleton(argc, argv);
        return singleton;
    }

    /*!
     * \brief Return the MPI rank of the running process.
     */
    static int rank()
    { return rank_; }

    /*!
     * \brief Return the number of MPI processes of the simulation.
     */
    static int size()
    { return 1; }

    /*!
     * \brief Returns the communicator which can be used to talk to
     *        peer processes.
     */
    static MPICommunicator getCommunicator()
    {
#if HAVE_MPI
        return MPI_COMM_WORLD;
#else
        static MPICommunicator comm;
        return comm;
#endif
    }

    /*!
     * \brief Returns the communicator which can be used for soliloquizing.
     */
    static MPICommunicator getLocalCommunicator()
    {
#if HAVE_MPI
        return MPI_COMM_SELF;
#else
        return getCommunicator();
#endif
    }

    /*!
     * \brief Returns an object which can be used for collective
     * communications.
     */
    static Dune::CollectiveCommunication<MPICommunicator>
    getCollectiveCommunication()
    {
        return Dune::CollectiveCommunication<MPICommunicator>(getCommunicator());
    }

private:
    MpiHelper &operator=(const MpiHelper);
    static int size_;
    static int rank_;
};

int MpiHelper::size_ = 1;
int MpiHelper::rank_ = 0;

} // namespace Ewoms

// this is a bit hacky, but required for DUNE compatibility
namespace Dune {
typedef Ewoms::MpiHelper MPIHelper;
}

#endif
