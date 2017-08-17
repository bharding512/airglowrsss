/***************************************************************************
                          calibration.h  -  description
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

#ifndef CALIBRATION_H
#define CALIBRATION_H


/**
  *@author Jonathan Makela
  */

class Calibration {
public: 
	Calibration();
	~Calibration();
  /** No descriptions */
  void SetCoefficients();
public: // Public attributes
  /**  */
  double a0;
  double a1;
  double a2;
  double b0;
  double b1;
  double b2;
  bool myIsAllSky;
  double myCenterAz;
  double myCenterEl;
  double pi;
};

#endif
