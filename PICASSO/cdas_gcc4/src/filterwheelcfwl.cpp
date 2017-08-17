/***************************************************************************
                          filterWheelCFWL.cpp  -  description
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

#include "filterwheelcfwl.h"
#include "system.h"
#include "sbigudrv.h"

#ifndef INVALID_HANDLE_VALUE
  #define INVALID_HANDLE_VALUE -1
#endif

extern System mySystem;

FilterWheelCFWL::FilterWheelCFWL(){
  myFilterOK = true;
  myNumFilters = 5;
  myShutterPos = 1;

  if(mySystem.myConfiguration.myCameraType != SBIGUNIV)
  {
    // The camera isn't set up, so fail
    myFilterOK = false;
    mySystem.myConfiguration.myFilterWheelType = ERROR;
  }
  else
  {
    CloseShutter();
  }
}

FilterWheelCFWL::~FilterWheelCFWL(){
  myFilterOK = false;
  CloseShutter();
}

/** Read filterwheel class from file */
void FilterWheelCFWL::Read(ifstream & is){
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
void FilterWheelCFWL::Write(ofstream & os){
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
void FilterWheelCFWL::WriteXML(ofstream & os){
}

/** Closes the filterwheel shutter. */
void FilterWheelCFWL::CloseShutter(){
  // Send the command to move to the shutter closed position
  MoveToPosition(myShutterPos);

  while(IsShutterOpen())
    sleep(1);
}

/** Returns the current position of the filterwheel. */
int FilterWheelCFWL::GetCurrentPosition(){
  CFWParams cfwp;
  CFWResults cfwr;

  // Setup the parameters for the request
  cfwp.cfwModel = 4;    // 4 is a CFW-L filterwheel
  cfwp.cfwCommand = 0;  // 0 is the query command

  short res;

  // Execuate the command
  res = SBIGUnivDrvCommand(CC_CFW, &cfwp, &cfwr);

  if(res != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot read filter position: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myFilterOK = false;
    return res;
  }

  if(cfwr.cfwStatus != CFWS_IDLE)
  {
    // The filterwheel is still moving
    return -999;
  }

  return cfwr.cfwPosition;
}

/** Returns whether the shutter is open or closed. */
bool FilterWheelCFWL::IsShutterOpen(){
  if(GetCurrentPosition() == myShutterPos || GetCurrentPosition() == 15)
    return false;
  else
    return true;
}

/** Moves the filterwheel to the position given by pos. */
int FilterWheelCFWL::MoveToPosition(int pos){
  CFWParams cfwp;

  // Check to make sure the filter position is within bounds
  if(pos < 1 || pos > myNumFilters)
    return ERROR;

  // Setup the parameters for the request
  cfwp.cfwModel = 4;    // 4 is a CFW-L filterwheel
  cfwp.cfwCommand = 1;  // 1 is the Goto command
  cfwp.cwfParam1 = pos; // Send it to the pos

  // The results structure
  CFWResults cfwr;

  short res;

  // Execute the command
  res = SBIGUnivDrvCommand(CC_CFW, &cfwp, &cfwr);

  if(res != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot change filter position: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myFilterOK = false;
    return myLastError;
  }

    // Setup the parameters for the request
  cfwp.cfwModel = 4;    // 4 is a CFW-L filterwheel
  cfwp.cfwCommand = 0;  // 0 is the query command

  // Execuate the command
  res = SBIGUnivDrvCommand(CC_CFW, &cfwp, &cfwr);

  if(res != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot read filter position: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myFilterOK = false;
    return res;
  }

  if(cfwr.cfwStatus != CFWS_IDLE)
  {
    // The filterwheel is still moving
    return -999;
  }

  // It is done moving, so return the position
  return cfwr.cfwPosition;    
}

/** Initializes the filterwheel. */
bool FilterWheelCFWL::SetUp_FilterWheel(){
}

/** Bottleneck function for all calls to the driver that logs the
    command and error.  First it activates our handle and then it
    calls the driver.  Activating the handle first allows having
    multpile instances of this class dealing with multple cameras
    on different communications port.

    Also allows direct access to the SBIG Universal Driver after the driver has been opened.  */
PAR_ERROR FilterWheelCFWL::SBIGUnivDrvCommand(short command, void *Params, void *Results){
  SetDriverHandleParams sdhp;

  // make sure we have a valid handle to the driver
  myLastCommand = (PAR_COMMAND) command;
  if(mySystem.myConfiguration.myCamera->GetCamHandle() == INVALID_HANDLE_VALUE)
    myLastError = CE_DRIVER_NOT_OPEN;
  else
  {
    // hanlde is valid so install it in the driver
    sdhp.handle = mySystem.myConfiguration.myCamera->GetCamHandle();
    if((myLastError = (PAR_ERROR)::SBIGUnivDrvCommand(CC_SET_DRIVER_HANDLE, &sdhp, NULL)) == CE_NO_ERROR)
      // call the desired command
      myLastError = (PAR_ERROR)::SBIGUnivDrvCommand(command, Params, Results);
  }

  return myLastError;
}

/** Return a string describing the last error encountered */
string FilterWheelCFWL::GetErrorString(){
  GetErrorStringParams gesp;
  GetErrorStringResults gesr;
  string s;

  gesp.errorNo = myLastError;
  SBIGUnivDrvCommand(CC_GET_ERROR_STRING, &gesp, &gesr);
  s = gesr.errorString;

  return s;  
}

/** Sets the filter position containing the shutter. */
void FilterWheelCFWL::SetShutterPosition(int num){
  myShutterPos = num;
}
/** Returns the shutter position. */
int FilterWheelCFWL::GetShutterPosition(){
  return myShutterPos;
}


/*!
    \fn FilterWheelCFWL::OpenShutter()
 */
void FilterWheelCFWL::OpenShutter()
{
    // The CFWL uses a filter blank as a shutter, so there is nothing to do.
    return;
}

void FilterWheelCFWL::SetNumFilters(int num){
  myNumFilters = num;
}

int FilterWheelCFWL::GetNumFilters(){
  return myNumFilters;
}
