/***************************************************************************
                          camerasbiguniv.h  -  description
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

#ifndef CAMERASBIGUNIV_H
#define CAMERASBIGUNIV_H

#ifndef _PARDRV_
  #include "sbigudrv.h"
#endif

#include <camera.h>

/**Implements the SBIG Universal Driver/Library as described in the version 4.35 documentation.
  *@author Jonathan Makela
  */

class CameraSBIGUniv : public Camera  {
public: 
	CameraSBIGUniv();
	~CameraSBIGUniv();
  /** Return a string describing the last error encountered */
  string GetErrorString();
  /** Reads the temperature status and fills the variables */
  void ReadTemperatureStatus();
  /** Sets the setpoint for the CCD */
  void SetCCDTemp(double CCDTemp);
  /** Gets the flag on temperature regulation */
  bool GetTemperatureRegulation();
  /** Returns the current CCD temperature */
  double GetCCDTemp();
  /** Reads the camera information from the configuration file */
  void Read(ifstream &is);
  /** Gets the fan power on the CCD */
  int GetFanPower();
  /** Returns the ambient temperature on the CCD chip */
  double GetAmbientTemp();
  /** Returns the estimated amount of time before the current exposure is done */
  int GetExposureTimeLeft();
  /** Writes the camera information to the configuration file */
  void Write(ofstream &os);
  /** Writes the XML information out to be read by the webserver */
  void WriteXML(ofstream &os);
  /** Sets the width of the image */
  void SetWidth(int width);
  /** Sets binning in the y-direction only */
  void SetBinY(int Y);
  /** Sets the binning in the x-direction only */
  void SetBinX(int X);
  /** Sets the exposure time */
  void SetExpTime(float exp);
  /** Set the offsets from the top-left corner of the CCD */
  void SetOffsets(int top, int left);
  /** Set the binning in the x and y directions */
  void SetBinning(int X, int Y);
  /** Set the size of the image */
  void SetImageSize(int width, int height);
  /** Returns the binning in the y direction */
  int GetBinY();
  /** Returns the binning in the x direction */
  int GetBinX();
  /** Returns the top offset of the image */
  int GetTop();
  /** Returns the left offset of the image */
  int GetLeft();
  /** Returns the height of the image */
  int GetHeight();
  /** Returns the width of the image */
  int GetWidth();
  /** Sets the top offset of the CCD */
  void SetTop(int top);
  /** Sets the left offset of the CCD */
  void SetLeft(int left);
  /** Sets the height of the image */
  void SetHeight(int height);
  /** Returns the current exposure time */
  float GetExposureTime();
  /** Sets the dark image flag */
  void SetDark(bool dark);
  /** Returns the CCD temperature setpoint */
  double GetCCDSetTemp();
  /** Captures an image given the input parameters:
	*myImage - apointer to the CImage structure in which the image data will be saved
	filename - the filename for the image to be saved */
  bool CaptureImage(Image *myImage, string filename, bool singleimage);
  /** End or abort an exposure in the camera.  Should be matched with a StartExposure call */
  PAR_ERROR EndExposure();
  /** Discards lines from the top of the CCD */
  PAR_ERROR DumpLines(int lines);
  /** Readout a line of data from the camera, optionally performing a dark subtraction, placing the data at dest. */
  PAR_ERROR ReadoutLine(ReadoutLineParams rlp, MY_LOGICAL darkSubtract, unsigned short *dest);
  /** Returns the camera handle. */
  short GetCamHandle();
  /** Sets the myDarkData variable. */
  void SetDarkData(IMAGETYPE * data, int width, int height);
  /** Returns the pointer to myDarkData. */
  IMAGETYPE * GetDarkDataPointer();
  /** Gets the USB gain for the SBIG camera. */
  float GetGain();
  /** Sets the gain for the SBIG USB device. */
  void SetGain(float gain);
  /** Turns the CCD temperature regulation off. */
  void SetTemperatureRegulationOff();
  /** Turns the temperature regulation of the CCD on. */
  void SetTemperatureRegulationOn();
  /** No descriptions */
  void SetAutoTemp(bool autoReg);
  /** No descriptions */
  bool GetAutoTemp();
  /** Sets the desired run temperature for the CCD system. */

  /** Returns the desired idle temperature of the CCD. */
//Del by KDevelop: 
  /** Returns the desired running temperature of the system. */
//Del by KDevelop: 
  /** Sets the idle temperature for the CCD.  This is the temperature the CCD will be held at when the system is not running during AUTO sequencing. */
//Del by KDevelop: 
public: // Public attributes
  /** Flag for whether tha camera is initialized (TRUE) or not (FALSE) */
  bool myCameraOK;
//  /** Number of rows and columns in the image */
  int myRows;
  int myCols;
  /** The handle to the camera currently being used.  For now, we only support one camera at a time. */
  short myCameraHandle;
  /**  */
  float myUSBgain;
  /**  */
  bool isTempRegulating;
  /**  */
  bool myAutoTempRegulation;
private: // Private attributes
  /** Number of rows and columns in the image */
//  int myRows;
//  int myCols;
  /** Saves the data from the previous dark image. */
  IMAGETYPE* myDarkData;
  /** Holds the temperature that the CCD will try to maintain */
  double myCCDSetTemp;
  /** The current temperature of the CCD */
  double myCCDTemp;
  /** The current ambient temperature (not implemented in most versions of the SBIG camera) */
  double myAmbientTemp;
  /** The current power to the CCD cooling fan */
  int myFanPower;
  /** The X and Y binning to be used */
  int myBinX;
  int myBinY;
  /** The corners of the imaged area of the CCD */
  int myTop;
  int myLeft;
  int myWidth;
  int myHeight;
  /** The exposure time to be used for the image */
  float myExposureTime;
  /** The exposure time remaining */
  int myExposureTimeLeft;
  /** Holds a reference to the last error encountered */
  PAR_ERROR myLastError;
  /** Holds a reference to the last command given */
  PAR_COMMAND myLastCommand;
  /** The type of camera we are connected to */
  CAMERA_TYPE myCameraType;
  /** Determines which CCD to use (Imaging or Tracking) */
  int myCCD;
  /** Determines the anti-blooming gate state */
  int myABGState;
  /** Determines the shutter mode on the CCD */
  int myOpenShutter;
  /** Determines the readout mode for the CCD */
  int myReadoutMode;
  /** Determines whether the CCD should perform subtraction */
  int mySubtract;
  /** Contains the results of calls to GetDriverInfo */
  GetDriverInfoResults0 myGDIR0;
  /** The parameters controlling the fan, shutter, and LED */
  MiscellaneousControlParams myMCP;
  /** Contains information from querying the temperature status */
  QueryTemperatureStatusResults myQTSR;
  /** Contains results from GetCCDInfo calls */
  GetCCDInfoResults0 myGCCDIR0;
  /** Flag for if a dark image is to be taken */
  bool myDark;
  /** A pointer to hold the data read from the camera */
  IMAGETYPE *myCameraData;
private: // Private methods
  /** Initializes the camera */
  bool InitCamera();
  /** Abort the current exposure */
  void KillExposure();
  /** Bottleneck function for all calls to the driver that logs the command and error.  First it activates our handle and then it calls the driver.  Activating the handle first allows having multpile instances of this class dealing with multple cameras on different communications port.

Also allows direct access to the SBIG Universal Driver after the driver has been opened. */
  PAR_ERROR SBIGUnivDrvCommand(short command, void *Params, void *Results);
  /** Open the drive.  Must be made before any other calls and should be valled only once per instance of the camera class.  Based on the results of the open call to the driver, this can open a new handle to the driver. */
  PAR_ERROR OpenDriver();
  /** Call once to open a particular port (USB, LPT, ETH, etc).  Must be balanced with a call to CloseDevice. */
  PAR_ERROR OpenDevice(OpenDeviceParams odp);
  /** Closes whichever device was opened by OpenDriver */
  PAR_ERROR CloseDevice();
  /** Should be called for every call to OpenDriverr.  Standard destructor does this for you as well.  Closing the driver multiple times will not hurt but will return an error. */
  PAR_ERROR CloseDriver();
  /** Shuts down the camera by closing the device and driver */
  int ShutDownCamera();
  /** Once the driver and device are open, call this to establish a communications link with the camera. */
  PAR_ERROR EstablishLink();
  /** Get the requested driver info for the passed request.  This call only works with the DRIVER_STANDARD and DRIVER_EXTENDED requests as you pass it a results reference tat only works with those 2 requests. */
  PAR_ERROR GetDriverInfo(DRIVER_REQUEST request);
  /** Converts from AD temperatures to degrees C */
  double ADToDegreesC(unsigned short ad, MY_LOGICAL ccd);
  /** Convert temperature in degrees C to camera AD setpoint */
  unsigned short DegreesCToAD(double degC, MY_LOGICAL ccd);
  /** If a link has been established to a camera, return TRUE.  Otherwise try to establish a link and if successful, return TRUE.  If fails, return FALSE */
  MY_LOGICAL CheckLink();
  /** Actually performs the exposure */
  bool DoExposureCollect();
  /** Start an exposure in the camera.  Should be matched with an EndExposure call */
  PAR_ERROR StartExposure();
  /** Query the camera to see if the exposure in progress is complete */
  PAR_ERROR IsExposureComplete(MY_LOGICAL &complete);
  /** Start the readout process.  This should be called after EndExposure and should be matched with an EndReadout call. */
  PAR_ERROR StartReadout(StartReadoutParams srp);
  /** Unfreezes the cooler. */
  PAR_ERROR UnfreezeCooler();
  /** Freezes the TE cooler for quieter readouts */
  PAR_ERROR FreezeCooler();
  /** End a readout started with StartReadout.  Don't forget to make this call to prepare the CCD for idling. */
  PAR_ERROR EndReadout();
};

#endif
