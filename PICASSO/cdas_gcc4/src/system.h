/***************************************************************************
                          system.h  -  description
                             -------------------
    begin                : Wed Aug 6 2003
    copyright            : (C) 2003 by Jonathan Makela
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

#ifndef SYSTEM_H
#define SYSTEM_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "configuration.h"

using namespace std;

/**The System Class is used to contain information about the status of the entire camera system, such as the location, the filters, the camera, etc...
  *@author Jonathan Makela
  */

class System {
public: 
	System();
	~System();
  
public: // Public attributes
  /** The configuration being used. */
  Configuration myConfiguration;
  /**  */
  string myConfFileName;
  /**  */
  bool myIsRunning;
  /**  */
  bool myServerRunning;
  /**  */
  bool myControlRunning;
  /**  */
  bool myIsExposing;
  /**  */
  string myXMLFilename;
  
};

#endif
