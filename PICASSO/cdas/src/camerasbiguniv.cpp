/***************************************************************************
                          camerasbiguniv.cpp  -  description
                             -------------------
    begin                : Mon Apr 26 2004
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

#include "camerasbiguniv.h"
#include "system.h"
#include <sstream>
#include "sbigudrv.h"

#ifndef INVALID_HANDLE_VALUE
  #define INVALID_HANDLE_VALUE -1
#endif

#define MAX_DIG_WIDTH 4096

// Temperature conversion constants
// Defined in the SBIG Universal driver documentation
#define T0      25.0
#define R0      3.0
#define DT_CCD  25.0
#define DT_AMB  45.0
#define RR_CCD  2.57
#define RR_AMB  7.791
#define RB_CCD  10.0
#define RB_AMB  3.0
#define MAX_AD  4096

extern System mySystem;

CameraSBIGUniv::CameraSBIGUniv(){
  // Initialize variables
  myCameraOK = false;
  myCCDSetTemp = 0.;
  myCameraHandle = 0;
  myCCDTemp = -999.;
  myAmbientTemp = -999.;
  myFanPower = -999;
  myBinX = 1;
  myBinY = 1;
  myTop = 0;
  myLeft = 0;
  myWidth = 1;
  myHeight = 1;
  myExposureTime = 0.;
  myLastError = CE_NO_ERROR;
  myLastCommand = CC_NULL;
  myCameraType = NO_CAMERA;
  myCCD = CCD_IMAGING;
  myABGState = ABG_CLK_MED7;
  myOpenShutter = 2;
  myReadoutMode = 0;
  myDarkData = NULL;
  myUSBgain = 1.2;
  myAutoTempRegulation = false;

  if(!InitCamera())
  {
    // There was an error
    mySystem.myConfiguration.myCameraType = ERROR;
  }
  else
    mySystem.myConfiguration.myCameraType = SBIGUNIV;
}

CameraSBIGUniv::~CameraSBIGUniv(){
  if(myDarkData != NULL)
    free(myDarkData);
  
  // Destroy the CameraSBIGUniv object
  if(mySystem.myIsExposing)
    KillExposure();

  ShutDownCamera();
  mySystem.myConfiguration.myCameraType = ERROR;
}

/** Initializes the camera */
bool CameraSBIGUniv::InitCamera(){
  // First open the device and get the camera handle
  OpenDeviceParams odp;

//  InitCamera();
  odp.deviceType = DEV_USB;
  if(OpenDriver() == CE_NO_ERROR)
    myLastError = OpenDevice(odp);
  else {
      // There was an error
      string str("Cannot open SBIG USB camera: ");
      str = str + GetErrorString();
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
  }

  // Establish a link
  if(EstablishLink() != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot establish link to SBIG USB camera: ");
    str = str + GetErrorString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Set safe defaults for the CCD imaging region
  myCCD = 0;        // We only have an imaging CCD
  myABGState = ABG_CLK_MED7;   // ABG shut off
  myOpenShutter = 2;  // Default to keep the shutter closed (for saftety reasons)
  myExposureTime = MIN_ST7_EXPOSURE;
  myReadoutMode = 0;
  myTop = 0;
  myLeft = 0;
  myWidth = 0;
  myHeight = 0;
  mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (InitCamera)");
  mySubtract = 0;

  // Get the driver information using a standard request
  if(GetDriverInfo((DRIVER_REQUEST) 0) != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot get driver info for SBIG USB camera: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Set the miscellaneous control parameters
  myMCP.fanEnable = TRUE;   // Default to having the fan run
  myMCP.ledState = 0;       // Keep the LED off
  myMCP.shutterCommand = 2; // Default have the shutter closed
  if(SBIGUnivDrvCommand(CC_MISCELLANEOUS_CONTROL, &myMCP, NULL) != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot setting miscellaneous control info for SBIG USB camera: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Get CCD Info
  GetCCDInfoParams gcip;
  gcip.request = (myCCD == CCD_IMAGING) ? 0:1;
  if(SBIGUnivDrvCommand(CC_GET_CCD_INFO, &gcip, &myGCCDIR0) != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot get CCD info for SBIG USB camera: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Set the camera gain
  // This should be able to be set in the FLTKCLIENT package, but not yet
  USBADControlParams ucp;
  ucp.command = USB_AD_IMAGING_GAIN;
//  ucp.data = myUSBgain;  // gain = 6.0 / (1.0 + 5.0*((63-data)/63))
  ucp.data = (ushort) 63 - (6.0/myUSBgain -1) * 63/5.0;
  if(SBIGUnivDrvCommand(CC_USB_AD_CONTROL, &ucp, NULL) != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot set the gain for SBIG USB camera: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
  }

  // We have sucessfully initialized the camera
  myCameraOK = true;

  // Set the CCD temperature
  SetCCDTemp(myCCDSetTemp);

  // Find the initial state of the temperature status
  ReadTemperatureStatus();
  
  return myCameraOK;
}

/** Abort the current exposure */
void CameraSBIGUniv::KillExposure(){
  if(EndExposure() != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot kill exposure: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
  } 
}

/** Bottleneck function for all calls to the driver that logs the
    command and error.  First it activates our handle and then it
    calls the driver.  Activating the handle first allows having
    multpile instances of this class dealing with multple cameras
    on different communications port.

    Also allows direct access to the SBIG Universal Driver after the driver has been opened. */
PAR_ERROR CameraSBIGUniv::SBIGUnivDrvCommand(short command, void *Params, void *Results){
  SetDriverHandleParams sdhp;

  // make sure we have a valid handle to the driver
  myLastCommand = (PAR_COMMAND) command;
  if(myCameraHandle == INVALID_HANDLE_VALUE)
    myLastError = CE_DRIVER_NOT_OPEN;
  else
  {
    // hanlde is valid so install it in the driver
    sdhp.handle = myCameraHandle;
    if((myLastError = (PAR_ERROR)::SBIGUnivDrvCommand(CC_SET_DRIVER_HANDLE, &sdhp, NULL)) == CE_NO_ERROR)
      // call the desired command
      myLastError = (PAR_ERROR)::SBIGUnivDrvCommand(command, Params, Results);
  }

  return myLastError;
}

/** Return a string describing the last error encountered */
string CameraSBIGUniv::GetErrorString(){
  GetErrorStringParams gesp;
  GetErrorStringResults gesr;
  string s;

  gesp.errorNo = myLastError;
  SBIGUnivDrvCommand(CC_GET_ERROR_STRING, &gesp, &gesr);
  s = gesr.errorString;

  return s;  
}

/** Open the drive.  Must be made before any other calls and
    should be valled only once per instance of the camera class.
    Based on the results of the open call to the driver, this can
    open a new handle to the driver. */
PAR_ERROR CameraSBIGUniv::OpenDriver(){
  short res;
  GetDriverHandleResults gdhr;
  SetDriverHandleParams sdhp;

  // call the driver directly so doesn't install our handle
  res = ::SBIGUnivDrvCommand(myLastCommand = CC_OPEN_DRIVER, NULL, NULL);
  if(res == CE_DRIVER_NOT_CLOSED)
  {
    // The driver is open already which we interpret as having been
    // opened by anoter instance of the class, so get the driver to
    // allocate a new handle and then record it
    sdhp.handle = INVALID_HANDLE_VALUE;
    res = ::SBIGUnivDrvCommand(CC_SET_DRIVER_HANDLE, &sdhp, NULL);
    if(res == CE_NO_ERROR)
    {
      res = ::SBIGUnivDrvCommand(CC_OPEN_DRIVER, NULL, NULL);
      if(res == CE_NO_ERROR)
      {
        res = ::SBIGUnivDrvCommand(CC_GET_DRIVER_HANDLE, NULL, &gdhr);
        if(res == CE_NO_ERROR)
          myCameraHandle = gdhr.handle;
      }
    }
  }
  else if(res == CE_NO_ERROR)
  {
    // The driver was not open, so record the driver handle so we
    // can support multiple instances of this class talking to
    // multiple cameras
    res = ::SBIGUnivDrvCommand(CC_GET_DRIVER_HANDLE, NULL, &gdhr);
    if(res == CE_NO_ERROR)
      myCameraHandle = gdhr.handle;
  }

  return myLastError = (PAR_ERROR) res;
}

/** Call once to open a particular port (USB, LPT, ETH, etc).  Must be
    balanced with a call to CloseDevice. */
PAR_ERROR CameraSBIGUniv::OpenDevice(OpenDeviceParams odp){
  return SBIGUnivDrvCommand(CC_OPEN_DEVICE, &odp, NULL);
}

/** Shuts down the camera by closing the device and driver */
int CameraSBIGUniv::ShutDownCamera(){
  CloseDevice();
  CloseDriver();
  myCameraOK = false;

  return 0;  
}

/** Should be called for every call to OpenDriverr.  Standard destructor
    does this for you as well.  Closing the driver multiple times will
    not hurt but will return an error. */
PAR_ERROR CameraSBIGUniv::CloseDriver(){
  PAR_ERROR res;

  res = SBIGUnivDrvCommand(CC_CLOSE_DRIVER, NULL, NULL);
  if(res == CE_NO_ERROR)
    myCameraHandle = INVALID_HANDLE_VALUE;

  return res;
}

/** Closes whichever device was opened by OpenDriver */
PAR_ERROR CameraSBIGUniv::CloseDevice(){
  return SBIGUnivDrvCommand(CC_CLOSE_DEVICE, NULL, NULL);
}

/** Once the driver and device are open, call this to establish a
    communications link with the camera. */
PAR_ERROR CameraSBIGUniv::EstablishLink(){
  PAR_ERROR res;
  EstablishLinkResults elr;
  EstablishLinkParams elp;

  res = SBIGUnivDrvCommand(CC_ESTABLISH_LINK, &elp, &elr);
  if(res == CE_NO_ERROR)
    myCameraType = (CAMERA_TYPE) elr.cameraType;

  return res;
}

/** Get the requested driver info for the passed request.
    This call only works with the DRIVER_STANDARD and DRIVER_EXTENDED
    requests as you pass it a results reference tat only works with those 2 requests. */
PAR_ERROR CameraSBIGUniv::GetDriverInfo(DRIVER_REQUEST request){
  GetDriverInfoParams gdip;

  gdip.request = request;
  myLastCommand = CC_GET_DRIVER_INFO;
  if(request > DRIVER_EXTENDED)
    return myLastError = CE_BAD_PARAMETER;
  else
    return SBIGUnivDrvCommand(CC_GET_DRIVER_INFO, &gdip, &myGDIR0);
}

/** Converts from AD temperatures to degrees C */
double CameraSBIGUniv::ADToDegreesC(unsigned short ad, MY_LOGICAL ccd){
  double r, degC;

  if(ad < 1)
    ad = 1;
  else if (ad >= MAX_AD -1)
    ad = MAX_AD - 1;
  if(ccd)
  {
    r = RB_CCD/(((double)MAX_AD/ad) - 1.0);
    degC = T0 - DT_CCD*(log(r/R0)/log(RR_CCD));
  }
  else
  {
    r = RB_AMB/(((double)MAX_AD/ad) - 1.0);
    degC = T0 - DT_AMB*(log(r/R0)/log(RR_AMB));
  }

  return degC;
}

/** Reads the temperature status and fills the variables */
void CameraSBIGUniv::ReadTemperatureStatus(){
  if(SBIGUnivDrvCommand(CC_QUERY_TEMPERATURE_STATUS, NULL, &myQTSR) != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot query temperature for SBIG USB camera: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
  }
  else
  {
    // Read the parameters
    myCCDTemp = ADToDegreesC(myQTSR.ccdThermistor, TRUE);
    myFanPower = (int) (myQTSR.power / 255.0);
    myAmbientTemp = ADToDegreesC(myQTSR.ambientThermistor, FALSE);
  }
}

/** Sets the setpoint for the CCD */
void CameraSBIGUniv::SetCCDTemp(double CCDTemp){
  
  // Place some error bounds on it
  if(CCDTemp > 20.0)
    CCDTemp = 20.0;

  if(myCameraOK)
  {
    SetTemperatureRegulationParams strp;
    if(CheckLink())
    {
      // We are linked and ready to go.  Update the temperature setpoint
      strp.ccdSetpoint = DegreesCToAD(CCDTemp, TRUE);
      strp.regulation = REGULATION_ON;
      if(SBIGUnivDrvCommand(CC_SET_TEMPERATURE_REGULATION, &strp, NULL) != CE_NO_ERROR)
      {
        // There was an error
        string str("Error setting CCD temperature setpoint for SBIG USB camera: ");
        str = str + string(GetErrorString());
        mySystem.myConfiguration.myLogger.WriteErrorString(str);
        return;
      }
      myCCDSetTemp = CCDTemp;
      isTempRegulating = true;
    }
  }
}

/** Convert temperature in degrees C to camera AD setpoint */
unsigned short CameraSBIGUniv::DegreesCToAD(double degC, MY_LOGICAL ccd){
  double r;
  unsigned short setpoint;

  if(degC < -50.0)
    degC = -50.0;
  else if(degC > 35.0)
    degC = 35.0;

  if(ccd)
  {
    r = R0 * exp(log(RR_CCD)*(T0-degC)/DT_CCD);
    setpoint = (unsigned short)(MAX_AD/((RB_CCD/r) + 1.0) + 0.5);
  }
  else
  {
    r = R0 * exp(log(RR_AMB)*(T0-degC)/DT_AMB);
    setpoint = (unsigned short)(MAX_AD/((RB_AMB/r) + 1.0) + 0.5);
  }

  return setpoint;
}

/** If a link has been established to a camera, return TRUE.
    Otherwise try to establish a link and if successful, return TRUE.
    If fails, return FALSE */
MY_LOGICAL CameraSBIGUniv::CheckLink(){
  if(myCameraType != NO_CAMERA || EstablishLink() == CE_NO_ERROR)
    return TRUE;
  else
    return FALSE;
}

/** Gets the flag on temperature regulation */
bool CameraSBIGUniv::GetTemperatureRegulation(){
  return isTempRegulating;
}

/** Returns the current CCD temperature */
double CameraSBIGUniv::GetCCDTemp(){
  return myCCDTemp;
}

/** Returns the estimated amount of time before the current exposure is done */
int CameraSBIGUniv::GetExposureTimeLeft(){
  return myExposureTimeLeft;
}

/** Returns the ambient temperature on the CCD chip */
double CameraSBIGUniv::GetAmbientTemp(){
  return myAmbientTemp;
}

/** Gets the fan power on the CCD */
int CameraSBIGUniv::GetFanPower(){
  return myFanPower;
}

/** Reads the camera information from the configuration file */
void CameraSBIGUniv::Read(ifstream &is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#Camera" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    mySystem.myConfiguration.myLogger.WriteErrorString("End of configuration file reached before #Camera found");
    return;
  }

  // Read the data
  is >> myCCDSetTemp;
  is >> mySystem.myConfiguration.myCameraType;
  is >> myBinX;
  is >> myBinY;
  is >> myHeight;
  mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (Read)");
  is >> myWidth;
  is >> myTop;
  is >> myLeft;
  is >> myExposureTime;
  is >> myDark;
  is >> myUSBgain;
  is >> myAutoTempRegulation;
}

/** Writes the camera information to the configuration file */
void CameraSBIGUniv::Write(ofstream &os){
  // Put the header label
  os << "#Camera" << endl;

  // Write the data
  os << myCCDSetTemp << endl;
  os << mySystem.myConfiguration.myCameraType << endl;
  os << myBinX << endl;
  os << myBinY << endl;
  os << myHeight << endl;
  os << myWidth << endl;
  os << myTop << endl;
  os << myLeft << endl;
  os << myExposureTime << endl;
  os << myDark << endl;
  os << myUSBgain << endl;
  os << myAutoTempRegulation << endl;
}

/** Writes the XML information out to be read by the webserver */
void CameraSBIGUniv::WriteXML(ofstream &os){
  // Write the header
  os << "<camera>" << endl;

  // Write the data
  os << "<ccdsettemp>" << myCCDSetTemp << "</ccdsettemp>" << endl;
  os << "<cameratype>" << mySystem.myConfiguration.myCameraType << "</cameratype>" << endl;
  os << "<binx>" << myBinX << "</binx>" << endl;
  os << "<biny>" << myBinY << "</biny>" << endl;
  os << "<imageheight>" << myHeight << "</imageheight>" << endl;
  os << "<imagewidth>" << myWidth << "</imagewidth>" << endl;
  os << "<imagetop>" << myTop << "</imagetop>" << endl;
  os << "<imageleft>" << myLeft << "</imageleft>" << endl;
  os << "<exposuretime>" << myExposureTime << "</exposuretime>" << endl;
  os << "<dark>" << myDark << "</dark>" << endl;
  os << "<usbgain>" << myUSBgain << "</usbgain>" << endl;
  os << "<regulation>" << myAutoTempRegulation << "</regulation>" << endl;

  os << "</camera>" << endl;
}

/** Set the size of the image */
void CameraSBIGUniv::SetImageSize(int width, int height){
  myWidth = width;
  myHeight = height;
mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (SetImageSize)");
}

/** Set the binning in the x and y directions */
void CameraSBIGUniv::SetBinning(int X, int Y){
  myBinX = X;
  myBinY = Y;
}

/** Set the offsets from the top-left corner of the CCD */
void CameraSBIGUniv::SetOffsets(int top, int left){
  myTop = top;
  myLeft = left;
}

/** Sets the exposure time */
void CameraSBIGUniv::SetExpTime(float exp){
  myExposureTime = exp;
}

/** Sets the binning in the x-direction only */
void CameraSBIGUniv::SetBinX(int X){
  myBinX = X;
}

/** Sets binning in the y-direction only */
void CameraSBIGUniv::SetBinY(int Y){
  myBinY = Y;
}

/** Sets the width of the image */
void CameraSBIGUniv::SetWidth(int width){
  myWidth = width;
}

/** Sets the height of the image */
void CameraSBIGUniv::SetHeight(int height){
  myHeight = height;
//  mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (SetHeight)");
}

/** Sets the left offset of the CCD */
void CameraSBIGUniv::SetLeft(int left){
  myLeft = left;
}

/** Sets the top offset of the CCD */
void CameraSBIGUniv::SetTop(int top){
  myTop = top;
}

/** Returns the width of the image */
int CameraSBIGUniv::GetWidth(){
  return myWidth;
}

/** Returns the height of the image */
int CameraSBIGUniv::GetHeight(){
  return myHeight;
}

/** Returns the left offset of the image */
int CameraSBIGUniv::GetLeft(){
  return myLeft;
}

/** Returns the top offset of the image */
int CameraSBIGUniv::GetTop(){
  return myTop;
}

/** Returns the binning in the x direction */
int CameraSBIGUniv::GetBinX(){
  return myBinX;
}

/** Returns the binning in the y direction */
int CameraSBIGUniv::GetBinY(){
  return myBinY;
}

/** Returns the CCD temperature setpoint */
double CameraSBIGUniv::GetCCDSetTemp(){
  return myCCDSetTemp;
}

/** Sets the dark image flag */
void CameraSBIGUniv::SetDark(bool dark){
  myDark = dark;
}

/** Returns the current exposure time */
float CameraSBIGUniv::GetExposureTime(){
  return myExposureTime;
}


/** Captures an image given the input parameters:
	*myImage - apointer to the CImage structure in which the image data will be saved
	filename - the filename for the image to be saved */
bool CameraSBIGUniv::CaptureImage(Image *myImage, string filename, bool singleimage)
{ 
  // First check to see if the camera is initialized
  if(!myCameraOK)
  {
    // Camera isn't initialized, so try to initialize it
    if(!InitCamera())
    {
      // There was an error
      string str("Error initialized camera: ");
      str = str + GetErrorString();
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

  // Update the temperature
  ReadTemperatureStatus();

  // Get the number of columns and rows in our image
  myRows = myHeight / myBinY;
  myCols = myWidth / myBinX;

  // Reflect this change in the image structure
  myImage->SetRows(myRows);
  myImage->SetCols(myCols);

  myCameraData = myImage->GetDataPointer();

  if(myDark)
  {
    // We want to perform a dark image
    myOpenShutter = 2;
  }
  else
  {
    // We want to perform a normal image
    myOpenShutter = 1;
  }

  // Figure out the readout mode
  // As a default use medium resolution mode
  myReadoutMode = 1;

  if(myBinX == 1 && myBinY == 1)
  {
    // High resolution mode
    myReadoutMode = 0;
  }
  if(myBinX == 2 && myBinY == 2)
  {
    // Medium resolution mode
    myReadoutMode = 1;
  }
  if(myBinX == 3 && myBinY == 3)
  {
    // Low resolution mode
    myReadoutMode = 2;
  }

  // Begin exposing
  mySystem.myIsExposing = true;

  // Set the times in the image header
//  mySystem.myConfiguration.myImage->SetTimesNow();
  myImage->SetTimesNow();

  // Write out to the logger
  ostringstream sout;
  sout << "Capturing " << filename << " w/ exposure time " << myExposureTime;
  mySystem.myConfiguration.myLogger.WriteSystemLogString(sout.str());
  
  if(!DoExposureCollect())
  {
    // Error in the exposure.  The error log should have been set in the function
    mySystem.myIsExposing = false;
    return false;
  }

  mySystem.myIsExposing = false;

  // See if we need to save the dark data
  if(myDark)
  {
    SetDarkData(myCameraData, myRows, myCols);
  }

  // Save the header information
  myImage->SetGain((int) (GetGain() * 10));
  myImage->SetCCDTemp((float)myCCDTemp);
  myImage->SetBinning(myBinX, myBinY);
  myImage->SetExposureTime(myExposureTime);
  myImage->SetSite(&mySystem.myConfiguration.mySite);
  myImage->SetFilterName(mySystem.myConfiguration.mySequencer.GetFilterName(mySystem.myConfiguration.mySequencer.GetCurrentFilter()));
  myImage->AutoContrast();
  if(singleimage)
    myImage->Normalize();
  else
    myImage->Normalize(!myDark);

  // Save the file to disk
//  ofstream myFile;
  string tifname;
  tifname = filename + ".tif";
//  myFile.open(tifname.c_str(), ios::out | ios::binary);
//  myImage->SaveTIF(myFile);
//  myFile.close();
  myImage->SaveTIF(tifname);

  return true;
}

/** Actually performs the exposure */
bool CameraSBIGUniv::DoExposureCollect(){
  // Make sure that the camera is still there
  if(!CheckLink())
  {
    // There was an error
    string str("Error linking to the SBIG camera: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check to make sure we're using a valid CCD
  if(myCCD > CCD_TRACKING)
  {
    // There was an error
    string str("Invalid CCD parameter: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check the anti-blooming parameter
  if(myABGState > ABG_CLK_HI7)
  {
    // There was an error
    string str("Invalid anti-blooming parameter: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check the shutter parameter
  if(myOpenShutter > SC_INITIALIZE_SHUTTER)
  {
    // There was an error
    string str("Invalid shutter parameter: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Get the CCD info
  GetCCDInfoParams gcip;
  gcip.request = (myCCD == CCD_IMAGING) ? 0:1;

  if(SBIGUnivDrvCommand(CC_GET_CCD_INFO, &gcip, &myGCCDIR0) != CE_NO_ERROR)
  {
    // There was an error
    string str("Error getting CCD info: ");
    str = str + GetErrorString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check the readout mode
  short myYBin = myReadoutMode >> 8;

  if(myCCD == CCD_TRACKING && myReadoutMode > 1)
  {
    // There was an error
    string str("Invalid readout mode for tracking CCD: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  if(myReadoutMode > 9)
  {
    // There is an error
    string str("Invalid readout mode: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Get the width & height of the imaged portion of the CCD
  unsigned short myCCDHeight, myCCDWidth;

  if(myReadoutMode >= 3 && myReadoutMode <= 5)
  {
    // Nx1,2,3 binning modes
    if(!myYBin)
    {
      // There was an error
      string str("Invalid Y binning parameter: ");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }

    myCCDWidth = myGCCDIR0.readoutInfo[myReadoutMode - 3].width;
    myCCDHeight = myGCCDIR0.readoutInfo[0].height >> myYBin;
  }
  else
  {
    myCCDWidth = myGCCDIR0.readoutInfo[myReadoutMode].width;
    myCCDHeight = myGCCDIR0.readoutInfo[myReadoutMode].height;
  }

  // Check the parameters
  if(myTop > (myCCDHeight - 1))
  {
    // There was an error
    string str("Image top off of CCD: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  if(myLeft > (myCCDWidth - 1))
  {
    // There was an error
    string str("Iamge left off of CCD: ");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  if(myWidth == 0)
  {
    // Use the entire CCD width
    myWidth = myCCDWidth;
    myLeft = 0;
  }
  else
  {
    if((myLeft + myCols) > myCCDWidth)
    {
      // There was an error
      string str("Image right off of CCD: ");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

  if(myHeight == 0)
  {
    // Use the entire CCD height
    myHeight = myCCDHeight;
    mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (DoExposureCollect)");
    myTop = 0;
  }
  else
  {
    if((myTop + myRows) > myCCDHeight)
    {
      // There was an error
      string str("Image bottom off of CCD: ");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

  // End any exposure that might have been in progress
  EndExposure();
  if(myLastError != CE_NO_ERROR && myLastError != CE_NO_EXPOSURE_IN_PROGRESS)
  {
    // There was an error
    string str("Error ending exposure: ");
    str = str + GetErrorString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    return false;
  }

  // Start the exposure
  if(StartExposure() != CE_NO_ERROR)
  {
     // There was an error
     string str("Error starting exposure: ");
     str = str + GetErrorString();
     mySystem.myConfiguration.myLogger.WriteErrorString(str);

     // Send an end exposure
     EndExposure();
     myCameraOK = false;
     return false;
  }

  // Calculate the start and stop times
  time_t start, stop, position;
  time(&start);
  stop = (time_t) (start + myExposureTime);

  // Monitor the exposure
  MY_LOGICAL expComp;
  PAR_ERROR err;
  do
  {
    time(&position);
    myExposureTimeLeft = (int) (stop - position);
  } while ((err = IsExposureComplete(expComp)) == CE_NO_ERROR && !expComp);

  EndExposure();
  if(err != CE_NO_ERROR)
  {
    // There was an error
    string str("Error ending the exposure: ");
    str = str + GetErrorString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Dump lines off of the CCD if necesary
  if(myTop > 0)
  {
    if(DumpLines(myTop) != CE_NO_ERROR)
    {
      // There was an error
      string str("Error dumping lines: ");
      str = str + GetErrorString();
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

  // Variables for calculation mean and std
  double sum1, sum2;
  sum1 = sum2 = 0.;
  unsigned short myLocalBuffer[MAX_DIG_WIDTH];

  // Start reading out lines
  StartReadoutParams srp;
  srp.ccd = myCCD;
  srp.left = (unsigned short) myLeft;
  srp.top = (unsigned short) myTop;
  srp.width = (unsigned short) myWidth;
  srp.height = (unsigned short) myHeight;
  srp.readoutMode = myReadoutMode;
  if((err = StartReadout(srp)) == CE_NO_ERROR)
  {
    int row, col;

    // Freeze the TE Cooler
    if(FreezeCooler() != CE_NO_ERROR)
    {
      // There was an error
      string str("Error freezing controller: ");
      str = str + GetErrorString();
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
 
    // Now readout each row
    ReadoutLineParams rlp;
    rlp.ccd = myCCD;
    rlp.pixelLength = myWidth;
    rlp.pixelStart = myLeft;
    rlp.readoutMode = myReadoutMode;

    // Begin the readout
    for(row = 0; row < myRows; row++)
    {
      if(ReadoutLine(rlp, FALSE, myLocalBuffer) != CE_NO_ERROR)
      {
        // There was an error
        string str("Error reading lines off of CCD: ");
        str = str + GetErrorString();
        mySystem.myConfiguration.myLogger.WriteErrorString(str);
        myCameraOK = false;

        // Unfreeze the CCD
        UnfreezeCooler();
        return false;
      }

      // Copy the data over into the image
      for(col = 0; col < myCols; col++)
      {
        myCameraData[row*myCols + col] = (IMAGETYPE) myLocalBuffer[col];
        sum1 += (double) myLocalBuffer[col];
        sum2 += (double) (((double) myLocalBuffer[col]) * myLocalBuffer[col]);
      }
    }
  }
  else
  {
    // There was an error
    string str("Error starting the readout: ");
    str = str + GetErrorString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // End the readout
  if(EndReadout() != CE_NO_ERROR)
  {
    // There was an error
    string str("Error ending the readout: ");
    str = str + GetErrorString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // calculate the mean and stddev
  double mean;
  mean = sum1 / (myRows * myCols);
  double temp;
  temp = sum2 / (myRows * myCols);
  double stddev;
  if(temp < mean * mean)
    stddev = 0.;
  else
    stddev = sqrt(temp - mean*mean);

  mySystem.myConfiguration.myImage->SetMean(mean);
  mySystem.myConfiguration.myImage->SetStdDev(stddev);

  // Finally, unfreeze the cooler
  if(UnfreezeCooler() != CE_NO_ERROR)
  {
    // There was an error
    string str("Error unfreezing the cooler: ");
    str = str + GetErrorString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  return true;    
}

/** Start an exposure in the camera.  Should be matched with an EndExposure call */
PAR_ERROR CameraSBIGUniv::StartExposure(){

  StartExposureParams sep;

  // Check the exposure time
  sep.exposureTime = (unsigned long) (myExposureTime * 100.0 + 0.5);
  if(sep.exposureTime < 1)
  {
    sep.exposureTime = 1;
    myExposureTime = ((float) 1 - 0.5)/100.;
  }
  sep.abgState = myABGState;
  sep.openShutter = myOpenShutter;
  sep.ccd = myCCD;

  if(CheckLink())
    return SBIGUnivDrvCommand(CC_START_EXPOSURE, &sep, NULL);
  else
    return myLastError;
}

/** End or abort an exposure in the camera.  Should be matched with a StartExposure call */
PAR_ERROR CameraSBIGUniv::EndExposure(){

  EndExposureParams eep;

  eep.ccd = myCCD;
  if(CheckLink() && myCameraOK)
    return SBIGUnivDrvCommand(CC_END_EXPOSURE, &eep, NULL);
  else
    return myLastError;
}

/** Query the camera to see if the exposure in progress is complete */
PAR_ERROR CameraSBIGUniv::IsExposureComplete(MY_LOGICAL &complete){
  QueryCommandStatusParams qcsp;
  QueryCommandStatusResults qcsr;

  complete = FALSE;
  if(CheckLink())
  {
    qcsp.command = CC_START_EXPOSURE;
    if(SBIGUnivDrvCommand(CC_QUERY_COMMAND_STATUS, &qcsp, &qcsr) == CE_NO_ERROR)
    {
      if(myCCD == CCD_IMAGING)
        complete = (qcsr.status & 0x03) != 0x02;
      else
        complete = (qcsr.status & 0x0C) != 0x08;
    }
  }

  return myLastError;
}
/** Discards lines from the top of the CCD */
PAR_ERROR CameraSBIGUniv::DumpLines(int lines){
    DumpLinesParams dlp;
    dlp.ccd = myCCD;
    dlp.readoutMode = myReadoutMode;
    dlp.lineLength = lines;

    if(CheckLink())
      return SBIGUnivDrvCommand(CC_DUMP_LINES, &dlp, NULL);
    else
      return myLastError;
}

/** Start the readout process.  This should be called after EndExposure
    and should be matched with an EndReadout call. */
PAR_ERROR CameraSBIGUniv::StartReadout(StartReadoutParams srp){
  if(CheckLink())
    return SBIGUnivDrvCommand(CC_START_READOUT, &srp, NULL);
  else
    return myLastError;
}
/** Freezes the TE cooler for quieter readouts */
PAR_ERROR CameraSBIGUniv::FreezeCooler(){
  SetTemperatureRegulationParams strp;
  strp.ccdSetpoint = DegreesCToAD(myCCDSetTemp, TRUE);
  strp.regulation = 3;

    if(CheckLink())
      return SBIGUnivDrvCommand(CC_SET_TEMPERATURE_REGULATION, &strp, NULL);
    else
      return myLastError;  
}

/** Unfreezes the cooler. */
PAR_ERROR CameraSBIGUniv::UnfreezeCooler(){
  SetTemperatureRegulationParams strp;
  strp.ccdSetpoint = DegreesCToAD(myCCDSetTemp, TRUE);
  strp.regulation = 4;

    if(CheckLink())
      return SBIGUnivDrvCommand(CC_SET_TEMPERATURE_REGULATION, &strp, NULL);
    else
      return myLastError;  
}

/** Readout a line of data from the camera, optionally performing a dark subtraction, placing the data at dest. */
PAR_ERROR CameraSBIGUniv::ReadoutLine(ReadoutLineParams rlp, MY_LOGICAL darkSubtract, unsigned short *dest){
  if(CheckLink())
  {
    if(darkSubtract)
      return SBIGUnivDrvCommand(CC_READ_SUBTRACT_LINE, &rlp, dest);
    else
      return SBIGUnivDrvCommand(CC_READOUT_LINE, &rlp, dest);
  }
  else
    return myLastError;
}

/** End a readout started with StartReadout.  Don't forget
    to make this call to prepare the CCD for idling. */
PAR_ERROR CameraSBIGUniv::EndReadout(){
  EndReadoutParams erp;

  erp.ccd = myCCD;

  if(CheckLink())
    return SBIGUnivDrvCommand(CC_END_READOUT, &erp, NULL);
  else
    return myLastError;
}
/** Returns the camera handle. */
short CameraSBIGUniv::GetCamHandle(){
  return myCameraHandle;
}
/** Sets the myDarkData variable. */
void CameraSBIGUniv::SetDarkData(IMAGETYPE * data, int width, int height){
  if(myDarkData != NULL)
    free(myDarkData);

  // allocate the array data
  myDarkData = (IMAGETYPE *) new BYTE[sizeof(IMAGETYPE)*width*height];

  // Copy over the data
  for(int i = 0; i < width; i++)
  {
    for(int j = 0; j < height; j++)
    {
        myDarkData[i*height + j] = data[i*height + j];
    }
  }
}
/** Returns the pointer to myDarkData. */
IMAGETYPE * CameraSBIGUniv::GetDarkDataPointer(){
  return myDarkData;
}
/** Sets the gain for the SBIG USB device. */
void CameraSBIGUniv::SetGain(float gain){
  if(gain < 1)
    gain = 1;
  if(gain > 6)
    gain = 6;
  
  myUSBgain = gain;

  // Set the camera gain
  // This should be able to be set in the FLTKCLIENT package, but not yet
  USBADControlParams ucp;
  ucp.command = USB_AD_IMAGING_GAIN;
//  ucp.data = myUSBgain;  // gain = 6.0 / (1.0 + 5.0*((63-data)/63))
  ucp.data = (ushort) 63 - (6.0/myUSBgain -1) * 63/5.0;
  if(SBIGUnivDrvCommand(CC_USB_AD_CONTROL, &ucp, NULL) != CE_NO_ERROR)
  {
    // There was an error
    string str("Cannot set the gain for SBIG USB camera: ");
    str = str + string(GetErrorString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
  }

}
/** Gets the USB gain for the SBIG camera. */
float CameraSBIGUniv::GetGain(){
  return myUSBgain;
}
/** Turns the temperature regulation of the CCD on. */
void CameraSBIGUniv::SetTemperatureRegulationOn(){
if(myCameraOK)
  {
    SetTemperatureRegulationParams strp;
    if(CheckLink())
    {
      // We are linked and ready to go.  Update the temperature setpoint
      strp.ccdSetpoint = DegreesCToAD(myCCDSetTemp, TRUE);
      strp.regulation = REGULATION_ON;
      if(SBIGUnivDrvCommand(CC_SET_TEMPERATURE_REGULATION, &strp, NULL) != CE_NO_ERROR)
      {
        // There was an error
        string str("Error setting CCD temperature setpoint for SBIG USB camera: ");
        str = str + string(GetErrorString());
        mySystem.myConfiguration.myLogger.WriteErrorString(str);
        return;
      }
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation ON");
      isTempRegulating = true;
    }
  }
}
/** Turns the CCD temperature regulation off. */
void CameraSBIGUniv::SetTemperatureRegulationOff(){
if(myCameraOK)
  {
    SetTemperatureRegulationParams strp;
    if(CheckLink())
    {
      // We are linked and ready to go.  Update the temperature setpoint
      strp.ccdSetpoint = DegreesCToAD(myCCDSetTemp, TRUE);
      strp.regulation = REGULATION_OFF;
      if(SBIGUnivDrvCommand(CC_SET_TEMPERATURE_REGULATION, &strp, NULL) != CE_NO_ERROR)
      {
        // There was an error
        string str("Error setting CCD temperature setpoint for SBIG USB camera: ");
        str = str + string(GetErrorString());
        mySystem.myConfiguration.myLogger.WriteErrorString(str);
        return;
      }
      	mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation OFF");
	isTempRegulating = false;
    }
  }
}
/** No descriptions */
void CameraSBIGUniv::SetAutoTemp(bool autoReg){
	myAutoTempRegulation = autoReg;
}
/** No descriptions */
bool CameraSBIGUniv::GetAutoTemp(){
	return myAutoTempRegulation;
}
