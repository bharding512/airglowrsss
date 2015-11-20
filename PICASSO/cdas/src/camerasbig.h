/***************************************************************************
                          camerasbig.h  -  description
                             -------------------
    begin                : Fri Jan 9 2004
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

#ifndef CAMERASBIG_H
#define CAMERASBIG_H

#include <camera.h>

#include "lsbigusbcam.h"
#include "lsbiglptcam.h"
#include "image.h"

/**Control routines for the SBIG camera
  *@author Jonathan Makela
  */

class CameraSBIG : public Camera  {
public: 
	CameraSBIG();
	virtual ~CameraSBIG();
  /** Gets the fan power on the CCD */
  int GetFanPower();
  /** Returns the ambient temperature on the CCD chip */
  double GetAmbientTemp();
  /** Reads camera information from the configuration file */
  void Read(ifstream & is);
  /** No descriptions */
  int GetWidth();
  /** No descriptions */
  void SetTop(int top);
  /** No descriptions */
  void SetLeft(int left);
  /** No descriptions */
  void SetHeight(int height);
  /** No descriptions */
  void SetWidth(int width);
  /** No descriptions */
  void SetBinY(int Y);
  /** No descriptions */
  void SetBinX(int X);
  /** No descriptions */
  void SetExpTime(float exp);
  /** No descriptions */
  void SetOffsets(int top, int left);
  /** No descriptions */
  void SetBinning(int X, int Y);
  /** No descriptions */
  void SetImagingSize(int width, int height);
  /** Writes the camera information to the configuration file */
  void Write(ofstream & os);
  /** No descriptions */
  int GetBinY();
  /** No descriptions */
  int GetBinX();
  /** No descriptions */
  int GetTop();
  /** No descriptions */
  int GetLeft();
  /** No descriptions */
  int GetHeight();
  /** No descriptions */
  double GetCCDSetTemp();
  /** No descriptions */
  void SetDark(bool dark);
  /** No descriptions */
  float GetExposureTime();
  /** No descriptions */
  void WriteXML(ofstream & os);
private: // Private methods
  /** Shuts down the camera */
  int ShutDownCamera();
  /** Sets the setpoint for the CCD */
  void SetCCDTemp(double CCDTemp);
  /** Aborts the current exposure */
  void KillExposure();
private: // Private methods
  /** Initializes the camera */
  bool InitCamera();
  /** Returns the estimated amount of time before the current exposure is done */
  int GetExposureTimeLeft();
private: // Private methods
  /** Returns the current CCD temperature */
  double GetCCDTemp();
private: // Private methods
  /** Actually performs the exposure */
  bool DoExposureCollect();
private: // Private methods
  /** Captures an image given the input parameters:
	*myImage - a pointer to the CImage structure in which the image data will be saved
	ExpTime - number of seconds for the current exposure */
  bool CaptureImage(Image * myImage, string filename);
  /** No descriptions */
  void ReadTemperatureStatus();
    /** Structure that holds input parameters for the SBIG USB camera */
  GrabImageParams myGIP;
  /** Structure to hold the Driver Info */
  GetDriverInfoParams myGDIP;
  /** Structure to hold the Driver Info Results */
  GetDriverInfoResults0 myGDIR0;
  /** Structure containing various control parameters */
  MiscellaneousControlParams myMCP;
  /** Structure containing CCD Info */
  GetCCDInfoParams myGCCDIP;
  /** Structure containing results of CCD Info */
  GetCCDInfoResults0 myGCCDIR0;
  /** Structre info for temperature regulation */
  SetTemperatureRegulationParams mySTRP;
  /** Structure for querying the temperature status of the CCD */
  QueryTemperatureStatusResults myQTSR;
private: // Private attributes
  /** A pointer to the SBIG camera */
  SbigUsbCam *myCameraHandle;
  /**  */
  double myCCDSetTemp;

  /**  */
  double myCCDTemp;
  /**  */
  double myAmbientTemp;
  /**  */
  int myFanPower;
  /**  */
  int myRows;
  /**  */
  float myExposureTime;
  /**  */
  int myCols;
  /**  */
  int myBinX;
  int myBinY;
  int myTop;
  int myLeft;
  int myWidth;
  int myHeight;
  /**  */
  IMAGETYPE *myCameraData;
  /**  */
  int myExposureTimeLeft;
  /**  */
  bool myDark;
public: // Public attributes
  /**  */
  bool myCameraOK;
};

#endif
