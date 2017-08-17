/***************************************************************************
                          filterwheel.h  -  description
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

#ifndef FILTERWHEEL_H
#define FILTERWHEEL_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

/**The FilterWheel class contains the routines to control the filterwheel.
  *@author Jonathan Makela
  */

#define SBIGCFWL 0
#define FLIFWL 1
#define KEOFW 2
#define KEOSER 3

class FilterWheel {
public: 
	FilterWheel();
	virtual ~FilterWheel();
  /** Read filterwheel class from file */
  virtual void Read(ifstream &is);
  /** Write the filter class to config file */
  virtual void Write(ofstream &os);
  /** No descriptions */
  virtual void WriteXML(ofstream &os);
  /** Returns whether the shutter is open or closed. */
  virtual bool IsShutterOpen();
  /** Initializes the filterwheel. */
  virtual bool SetUp_FilterWheel();
  /** Closes the filterwheel shutter. */
  virtual void CloseShutter();
  /** Moves the filterwheel to the position given by pos. */
  virtual int MoveToPosition(int pos);
  /** Returns the current position of the filterwheel. */
  virtual int GetCurrentPosition();
  /** No descriptions */
  virtual void SetShutterPosition(int num);
  /** Returns the shutter position. */
  virtual int GetShutterPosition();
  /** Opens the filterwheel shutter */
  virtual void OpenShutter();
  virtual void SetFWDevice(char *s);
  virtual void SetTempControlDevice(char *s);
  virtual void SetNumFilters(int num);
  virtual int GetNumFilters();

public: // Public attributes
  /**  */
  bool myFilterOK;
  /** The number of filters in the filterwheel */
  int myNumFilters;
  /** If using the SBIG CFW-L, we have a filter position dedicated to the shutter. */
  int myShutterPos;
  int myTempControlDevice;
  int myFWDevice;
};

#endif
