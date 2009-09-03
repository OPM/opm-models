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
#ifndef DUNE_PRIMARYVARIABLES_HH
#define DUNE_PRIMARYVARIABLES_HH

namespace Dune
{
struct PrimaryVariablesDefinition
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
        noPressure = 3,
        //velocity flags
        velocityW = 0,
        velocityNW = 1,
        velocityTotal = 2
    };
    int saturationType;
    int pressureType;
    int velocityType;
} defaultImplementation = {0, 0, 2};
}
#endif
