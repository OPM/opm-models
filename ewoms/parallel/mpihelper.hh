// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2013 by Andreas Lauser                                    *
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
 * \copydoc Ewoms::MpiHelper
 */
#ifndef EWOMS_MPI_HELPER_HH
#define EWOMS_MPI_HELPER_HH

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
 * should be compatible with Dune's MPIHelper class.
 */
class MpiHelper
{
public:
    enum{
#if HAVE_MPI
        isFake = false //!< Returns whether MPI is really available or not
#else
        isFake = true //!< Returns whether MPI is really available or not
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
        static bool wasInstanced = false;
        if (wasInstanced)
            DUNE_THROW(Dune::InvalidStateException,
                       "Ewoms::MpiHelper may only be instanciated once!");
        wasInstanced = true;

#if HAVE_MPI
        if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
            DUNE_THROW(Dune::InvalidStateException,
                       "Initialization of MPI failed!");
#endif
    }

    ~MpiHelper()
    {
#if HAVE_MPI
        MPI_Finalize();
#endif
    }


    /*!
     * \brief Return the MPI rank of the running process.
     */
    static int rank()
    {
#if HAVE_MPI
        int tmp;
        MPI_Comm_rank(MPI_COMM_WORLD, &tmp);
        return tmp;
#else
        return 0;
#endif
    }

    /*!
     * \brief Return the number of MPI processes of the simulation.
     */
    static int size()
    {
#if HAVE_MPI
        int tmp;
        MPI_Comm_size(MPI_COMM_WORLD, &tmp);
        return tmp;
#else
        return 1;
#endif
    }

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
    static Dune::CollectiveCommunication<MPICommunicator> getCollectiveCommunication()
    { return Dune::CollectiveCommunication<MPICommunicator>(getCommunicator()); }

private:
    MpiHelper& operator=(const MpiHelper);
};

}

// this is a bit hacky, but required for DUNE compatibility
#define DUNE_MPIHELPER
namespace Dune {
    typedef Ewoms::MpiHelper MPIHelper;
}

#endif
