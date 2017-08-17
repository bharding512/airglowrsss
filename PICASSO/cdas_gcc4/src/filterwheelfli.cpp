/***************************************************************************
                          filterWheelFLI.cpp  -  description
                             -------------------
    begin                : Fri Oct 21 2005
    copyright            : (C) 2004,5 by Jonathan Makela and Ethan Miller
    email                : jmakela@uiuc.edu, esmiller@uiuc.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "filterwheelfli.h"
#include "system.h"
#include "libfli.h"

#ifndef INVALID_HANDLE_VALUE
  #define INVALID_HANDLE_VALUE -1
#endif

extern System mySystem;

FilterWheelFLI::FilterWheelFLI(){
  myFilterOK = true;
  myNumFilters = 5;
  myShutterPos = 1;

  int i, j;
  char **list;
  long err;
  /** Poll FLI driver for attached USB Filterwheel devices */
  err = FLIList(FLIDOMAIN_USB | FLIDEVICE_FILTERWHEEL, &list);
  /** Reformat returned list into C strings with \0 terminators */
  for (i = 0; list[i] != NULL; i++) {
    for (j = 0; list[i][j] != '\0'; j++)
      if (list[i][j] == ';') {
        list[i][j] = '\0';
        break;
      }
  }
  /** Select first operational USB filterwheel for our use */
  for (i = 0; list[i] != NULL; i++) {
    err = FLIOpen(&filterwheelDevice, list[i], FLIDOMAIN_USB | FLIDEVICE_FILTERWHEEL);
    if (err == 0) {
      break;
    }
  }
	FLIFreeList(list);
  CloseShutter();
}

FilterWheelFLI::~FilterWheelFLI(){
  myFilterOK = false;
  CloseShutter();
	FLIClose(filterwheelDevice);
}

/** Read filterwheel class from file */
void FilterWheelFLI::Read(ifstream & is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#FilterWheel" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    mySystem.myConfiguration.myLogger.WriteErrorString("End of configuration file reached before #FilterWheel found");
    return;
  }

  // Read the data
  is >> mySystem.myConfiguration.myFilterWheelType;
  is >> myNumFilters;
  is >> myShutterPos;
  is >> myFWDevice;
  is >> myTempControlDevice;

  // Close the shutter now that we've read in the position
  CloseShutter();
}

/** Write the filter class to config file */
void FilterWheelFLI::Write(ofstream & os){
  // Put the header label
  os << "#FilterWheel" << endl;

  // Write the CCD settemp
  os << mySystem.myConfiguration.myFilterWheelType << endl;
  os << myNumFilters << endl;
  os << myShutterPos << endl;
  os << myFWDevice << endl;
  os << myTempControlDevice << endl;
}

/** No descriptions */
void FilterWheelFLI::WriteXML(ofstream & os){
}

/** Closes the filterwheel shutter. */
void FilterWheelFLI::CloseShutter(){
  // Send the command to move to the shutter closed position
  MoveToPosition(myShutterPos);

  while(IsShutterOpen())
    sleep(1);
}

/** Returns the current position of the filterwheel. */
int FilterWheelFLI::GetCurrentPosition() {
    int err;
    long currentposn;
    err = FLIGetFilterPos(filterwheelDevice, &currentposn);
    if (err != 0) {
      string str("Cannot read filter position: ");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myFilterOK = false;
      return -999;
    }
    return (int)currentposn+1;
}

/** Returns whether the shutter is open or closed. */
bool FilterWheelFLI::IsShutterOpen(){
  if(GetCurrentPosition() == myShutterPos)
    return false;
  else
    return true;
}

/** Moves the filterwheel to the position given by pos. */
int FilterWheelFLI::MoveToPosition(int pos){
    flidev_t dev;
    int err;
    long currentposn;
    err = FLISetFilterPos(filterwheelDevice, (long)pos-1);
    if (err != 0) {
      string str("Cannot change filter position: ");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myFilterOK = false;
      return -999;
    }
    err = FLIGetFilterPos(filterwheelDevice, &currentposn);
    if (currentposn != ((long)pos-1)) {
    	string str("New position not equal to requested position");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myFilterOK = false;
      return -999;
    }
    return (int)currentposn+1;
}

/** Initializes the filterwheel. */
bool FilterWheelFLI::SetUp_FilterWheel(){
}

/** Return a string describing the last error encountered */
string FilterWheelFLI::GetErrorString(){
  string s("General filterwheel error");
  return s;
}

/** Sets the filter position containing the shutter. */
void FilterWheelFLI::SetShutterPosition(int num){
  myShutterPos = num;
}
/** Returns the shutter position. */
int FilterWheelFLI::GetShutterPosition(){
  return myShutterPos;
}


/*!
    \fn FilterWheelFLI::OpenShutter()
 */
void FilterWheelFLI::OpenShutter()
{
  // The FLI filterwheel uses a filter blank as a shutter, so there is nothing to do.
  return;
}

void FilterWheelFLI::SetNumFilters(int num){
  myNumFilters = num;
}

int FilterWheelFLI::GetNumFilters(){
  return myNumFilters;
}
