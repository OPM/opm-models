/*
  Copyright (C) 2013-2014 by Andreas Lauser
  Copyright (c) 2013 by SINTEF ICT, Applied Mathematics.
  Copyright (c) 2013 by Uni Research AS

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
 *
 * \brief This file implements several wrapper classes around the
 *        opaque ERT data types.
 *
 * These class are shamelessly ripped-off from opm-core and are
 * required to make writing Eclipse files exception safe...
 */
#ifndef EWOMS_ERT_WRAPPERS_HH
#define EWOMS_ERT_WRAPPERS_HH

#if HAVE_ERT && HAVE_DUNE_CORNERPOINT

#include <ert/ecl/fortio.h>
#include <ert/ecl/ecl_endian_flip.h>
#include <ert/ecl/ecl_grid.h>
#include <ert/ecl/ecl_kw_magic.h>
#include <ert/ecl/ecl_kw.h>
#include <ert/ecl/ecl_sum.h>
#include <ert/ecl/ecl_util.h>
#include <ert/ecl/ecl_init_file.h>
#include <ert/ecl/ecl_file.h>
#include <ert/ecl/ecl_rst_file.h>

#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

namespace Ewoms {
/*!
 * \brief This is a smart pointer class for ERT's ecl_kw_type structure.
 */
template <typename T>
class ErtKeyword
{
public:
#if HAVE_ERT
    typedef ecl_kw_type ErtHandleType;
#else
    typedef int ErtHandleType;
#endif

    // don't allow copies for objects of this class
    ErtKeyword(const ErtKeyword &) = delete;

    // Default constructor
    ErtKeyword()
        : ertHandle_(0)
    {}

    /// Initialization from single-precision array.
    ErtKeyword(const std::string& name,
               const std::vector<float>& data)
        : ertHandle_(0)
    { set(name, data); }

    /// Initialization from double-precision array.
    ErtKeyword(const std::string& name,
               const std::vector<double>& data)
        : ertHandle_(0)
    { set(name, data); }

    /// Initialization from double-precision array.
    ErtKeyword(const std::string& name,
               const std::vector<int>& data)
        : ertHandle_(0)
    { set(name, data); }

    ~ErtKeyword()
    {
#if HAVE_ERT
        if (ertHandle_)
            ecl_kw_free(ertHandle_);
#endif
    }

    template <class DataElementType>
    void set(const std::string name, const std::vector<DataElementType>& data)
    {
#if HAVE_ERT
        if(ertHandle_) {
            ecl_kw_free(ertHandle_);
        }

        ertHandle_ = ecl_kw_alloc(name.c_str(),
                                  data.size(),
                                  ertType_());

        // number of elements to take
        const int numEntries = data.size();

        // fill it with values
        T* target = static_cast<T*>(ecl_kw_get_ptr(ertHandle()));
        for (int i = 0; i < numEntries; ++i) {
            target[i] = static_cast<T>(data[i]);
        }
#endif
    }

    ErtHandleType *ertHandle() const
    { return ertHandle_; }

private:
#if HAVE_ERT
    static ecl_type_enum ertType_()
    {
        if (std::is_same<T, float>::value)
        { return ECL_FLOAT_TYPE; }
        if (std::is_same<T, double>::value)
        { return ECL_DOUBLE_TYPE; }
        if (std::is_same<T, int>::value)
        { return ECL_INT_TYPE; }

        OPM_THROW(std::logic_error,
                  "Unhandled type for data elements in ErtKeyword");
    }
#endif

    ErtHandleType *ertHandle_;
};

/*!
 * \brief This is a smart pointer class for ERT's ecl_grid_type structure.
 *
 * The class is shamelessly ripped of from opm-core and is required to
 * make writing Eclipse files exception safe...
 */
class ErtGrid
{
    enum { dim = 3 };

public:
#if HAVE_ERT
    typedef ecl_grid_type ErtHandleType;
#else
    typedef int ErtHandleType;
#endif

    ErtGrid(const ErtGrid& ) = delete;

    /*!
     * \brief Create an ERT grid based an Opm::EclipseGrid.
     */
    ErtGrid(Opm::EclipseGridConstPtr eclGrid)
    {
#if HAVE_ERT
        std::vector<double> mapaxesData;
        std::vector<double> coordData;
        std::vector<double> zcornData;
        std::vector<int> actnumData;

        eclGrid->exportMAPAXES(mapaxesData);
        eclGrid->exportCOORD(coordData);
        eclGrid->exportZCORN(zcornData);
        eclGrid->exportACTNUM(actnumData);

        ErtKeyword<float> mapaxesKeyword("MAPAXES", mapaxesData);
        ErtKeyword<float> coordKeyword("COORD", coordData);
        ErtKeyword<float> zcornKeyword("ZCORN", zcornData);
        ErtKeyword<int> actnumKeyword("ACTNUM", actnumData);

        ertHandle_ = ecl_grid_alloc_GRDECL_kw(eclGrid->getNX(),
                                              eclGrid->getNY(),
                                              eclGrid->getNZ(),
                                              zcornKeyword.ertHandle(),
                                              coordKeyword.ertHandle(),
                                              actnumKeyword.ertHandle(),
                                              mapaxesKeyword.ertHandle());
#endif // HAVE_ERT && HAVE_DUNE_CORNERPOINT
    }

    ~ErtGrid()
    {
#if HAVE_ERT
        ecl_grid_free(ertHandle_);
#endif
    }


    /*!
     * \brief Save the grid to an .EGRID file.
     */
    void write(const std::string& fileName, int reportStepIdx)
    {
#if HAVE_ERT
        ecl_grid_fwrite_EGRID(ertHandle(), fileName.c_str());
#endif
    }

    ErtHandleType *ertHandle() const
    { return ertHandle_; }

private:
    ErtHandleType *ertHandle_;
};

} // namespace Ewoms

#endif // HAVE_ERT && HAVE_DUNE_CORNERPOINT

#endif // EWOMS_ERT_WRAPPERS_HH
