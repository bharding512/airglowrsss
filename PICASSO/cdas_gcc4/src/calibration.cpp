/***************************************************************************
                          calibration.cpp  -  description
                             -------------------
    begin                : Sat Jan 10 2004
    copyright            : (C) 2004 by Jonathan Makela
    email                : jmakela@nrl.navy.mil
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
 
#include "calibration.h"

Calibration::Calibration(){
  pi = 4. * atan(1.);
  myCenterEl = 89.9 * pi / 180.;
  myCenterAz = 0;
  myIsAllSky = true;
}
Calibration::~Calibration(){
}
/** No descriptions */
void Calibration::SetCoefficients(){
}
