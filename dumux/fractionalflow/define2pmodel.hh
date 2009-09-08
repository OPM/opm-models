// $Id:$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_DEFINE2PMODEL_HH
#define DUNE_DEFINE2PMODEL_HH

namespace Dune
{
struct DefineModel
{
    enum
    {
        //saturation flags
        saturationW = 0,
        saturationNW = 1,
        //pressure flags
        pressureW = 0,
        pressureNW =1,
        pressureGlobal = 2,
        //velocity flags
        velocityW = 0,
        velocityNW = 1,
        velocityTotal = 2
    };
    int saturationType;
    int pressureType;
    int velocityType;
    bool compressibility;
    DefineModel():saturationType(saturationW),pressureType(pressureW),velocityType(velocityTotal),compressibility(false)
    {}
};
}
#endif
