/***************************************************************************
                          filterwheel.cpp  -  description
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

#include "filterwheel.h"
#include "system.h"

extern System mySystem;

FilterWheel::FilterWheel(){
myFilterOK = false;
//strcpy(myFWDevice,"");
myNumFilters = 0;
myShutterPos = 0;
//strcpy(myTempControlDevice,"");
myFWDevice = 0;
myTempControlDevice=0;
}
FilterWheel::~FilterWheel(){
}

/** Write the filter class to config file */
void FilterWheel::Write(ofstream &os){
  // Put the header label
  os << "#FilterWheel" << endl;

  // Write the CCD settemp
  os << mySystem.myConfiguration.myFilterWheelType << endl;
  os << myNumFilters << endl;
  os << myShutterPos << endl;
  os << myFWDevice << endl;
  os << myTempControlDevice << endl;
}

/** Read filterwheel class from file */
void FilterWheel::Read(ifstream &is){
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
}
/** No descriptions */
void FilterWheel::WriteXML(ofstream &os){
}
/** Returns the current position of the filterwheel. */
int FilterWheel::GetCurrentPosition(){
}
/** Moves the filterwheel to the position given by pos. */
int FilterWheel::MoveToPosition(int pos){
}
/** Closes the filterwheel shutter. */
void FilterWheel::CloseShutter(){
}
/** Initializes the filterwheel. */
bool FilterWheel::SetUp_FilterWheel(){
}
/** Returns whether the shutter is open or closed. */
bool FilterWheel::IsShutterOpen(){
}
/** No descriptions */
void FilterWheel::SetShutterPosition(int num){
  myShutterPos = num;
}
/** Returns the shutter position. */
int FilterWheel::GetShutterPosition(){
  return myShutterPos;
}
/** Opens the filterwheel shutter. */
void FilterWheel::OpenShutter(){
}

void FilterWheel::SetFWDevice(char *s){
//   strcpy(myFWDevice,s);;
}

void FilterWheel::SetTempControlDevice(char *s){
//   strcpy(myTempControlDevice,s);
}

void FilterWheel::SetNumFilters(int num){
  myNumFilters = num;
}

int FilterWheel::GetNumFilters(){
  return myNumFilters;
}
