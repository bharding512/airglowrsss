/***************************************************************************
                          filterWheelCFWL.h  -  description
                             -------------------
    begin                : Wed Jun 16 2004
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

#ifndef FILTERWheelCFWL_H
#define FILTERWheelCFWL_H

#ifndef _PARDRV_
  #include "sbigudrv.h"
#endif

#include <filterwheel.h>

/**Implements the SBIG CFW-L 5-position filterwheel in the SBIG Universal Driver.
  *@author Jonathan Makela
  */

class FilterWheelCFWL : public FilterWheel  {
public:
	FilterWheelCFWL();
	~FilterWheelCFWL();
  /** No descriptions */
  void WriteXML(ofstream & os);
  /** Write the filter class to config file */
  void Write(ofstream & os);
  /** Read filterwheel class from file */
  void Read(ifstream & is);
  /** Initializes the filterwheel. */
  bool SetUp_FilterWheel();
  /** Moves the filterwheel to the position given by pos. */
  int MoveToPosition(int pos);
  /** Returns whether the shutter is open or closed. */
  bool IsShutterOpen();
  /** Returns the current position of the filterwheel. */
  int GetCurrentPosition();
  /** Closes the filterwheel shutter. */
  void CloseShutter();
  /** Return a string describing the last error encountered */
  string GetErrorString();
  /** Sets the filter position containing the shutter. */
  void SetShutterPosition(int num);
  /** Returns the shutter position. */
  int GetShutterPosition();
    void OpenShutter();
  void SetNumFilters(int num);
  int GetNumFilters();
public: // Public attributes
  /**  */
  bool myFilterOK;
  /** The number of filters in the current filterwheel. */
  int myNumFilters;
  /** The location of the shutter position. */
  int myShutterPos;
    int myTempControlDevice;
    int myFWDevice;
private: // Private methods
  /** Bottleneck function for all calls to the driver that logs the
    command and error.  First it activates our handle and then it
    calls the driver.  Activating the handle first allows having
    multpile instances of this class dealing with multple cameras
    on different communications port.

    Also allows direct access to the SBIG Universal Driver after the driver has been opened.  */
  PAR_ERROR SBIGUnivDrvCommand(short command, void *Params, void *Results);
  /** Holds a reference to the last error encountered */
  PAR_ERROR myLastError;
  /** Holds a reference to the last command given */
  PAR_COMMAND myLastCommand;
};

#endif
