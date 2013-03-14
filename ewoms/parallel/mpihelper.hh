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
        MPI_Init(&argc, &argv);
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
};

}

#endif
