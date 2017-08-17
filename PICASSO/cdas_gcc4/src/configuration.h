/***************************************************************************
                          configuration.h  -  description
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

#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "user.h"
#include "site.h"
#include "camera.h"
//#include "camerasbig.h"
#include "camerasbiguniv.h"
#include "cameraandor.h"
#include "camerafli.h"
//#include "cameraapogee.h"
//#include "camerapvcam.h"
#include "filterwheel.h"
#include "filterwheelcfwl.h"
#include "filterwheelfli.h"
#include "filterwheelkeo.h"
#include "filterwheelserial.h"
#include "logger.h"
#include "schedule.h"
#include "directories.h"
#include "image.h"
#include "sequencer.h"

using namespace std;

/**The Configuration class contains information on the system configuration.
  *@author Jonathan Makela
  */

class Configuration {
public: 
	Configuration();
	~Configuration();
  /** Copy constructor */
   Configuration(Configuration& c);
  /** Loads a configuration in from filename */
  bool LoadConfiguration(string filename);
  /** Saves the configuration information */
  bool SaveConfiguration(string filename);
  /** Saves the configuration out to an XML file. */
  bool SaveXML(string filename);
  /** Parses a line of XML.  Basically returns the string between <XXX>YYY</XXX>, returning YYY.  Does absolutely no error checking on if it is a valid XML string. */
  string ParseXML(string input);
public: // Public attributes
  /**  */
  /**  */
  User myUser;
  Site mySite;
  /**  */
  Directories myDirectories;
  /** pointer to the Camera */
  Camera *myCamera;
  /** the filter wheel. */
  FilterWheel *myFilterWheel;
  /** The logger being used */
  Logger myLogger;
  /** The schedule being used */
  Schedule mySchedule;
  /** A pointer to the image currently being used */
  Image *myImage;
  /**  */
  int myCameraType;
  int myFilterWheelType;
  /**  */
  Sequencer mySequencer;
};

#endif
