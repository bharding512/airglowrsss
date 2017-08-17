/***************************************************************************
                          camera.h  -  description
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

#ifndef CAMERA_H
#define CAMERA_H

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "image.h"
#include <math.h>

#define ERROR -999
//#define SBIGUSB 0
#define PVCAM 1
#define SBIGUNIV 0
#define ANDOR 2
#define FLIPROLINE 3
#define APOGEE 4

using namespace std;

/**The Camera class contains the control routines for the camera, such as how to take a picture, the size of the image, etc...
  *@author Jonathan Makela
  */

class Camera {
public: 
	Camera();
	virtual ~Camera();
  /** Reads camera information from the configuration file */
  virtual void Read(ifstream& is);
  /** Writes the camera information to the configuration file */
  virtual void Write(ofstream& os);
  /** Sets the setpoint for the CCD */
  virtual void SetCCDTemp(double CCDTemp);
  /** Gets the flag on temperature regulation */
  virtual bool GetTemperatureRegulation();
 public: // Public attributes
  /** boolean if the temperature is regulating or not */
  bool isTempRegulating;
  /**  */
  bool myCameraOK;
  /**  */
  double myCCDTemp;
  /** The handle to the camera currently being used.  For now, we only support one camera at a time. */
  short myCameraHandle;
  /** Saves the data from the previous dark image. */
  IMAGETYPE* myDarkData;
  /**  */
  float myUSBgain;
  /**  */
  bool myAutoTempRegulation;
  /**  */
  virtual void SetImagingSize(int width, int height);
 /** Returns the estimated amount of time before the current exposure is done */
  virtual int GetExposureTimeLeft();
  /** Initializes the camera */
  virtual bool InitCamera();
  /** Gets the fan power on the CCD */
  virtual int GetFanPower();
  /** Returns the ambient temperature on the CCD chip */
  virtual double GetAmbientTemp();
  /** No descriptions */
  virtual void SetBinning(int X, int Y);
  /** No descriptions */
  virtual void SetOffsets(int top, int left);
  /** No descriptions */
  virtual void SetExpTime(float exp);
  /** No descriptions */
  virtual void SetBinY(int Y);
  /** No descriptions */
  virtual void SetBinX(int X);
  /** No descriptions */
  virtual void SetHeight(int height);
  /** No descriptions */
  virtual void SetWidth(int width);
  /** No descriptions */
  virtual void SetLeft(int left);
  /** No descriptions */
  virtual void SetTop(int top);
  /** No descriptions */
  virtual int GetWidth();
  /** No descriptions */
  virtual int GetHeight();
  /** No descriptions */
  virtual int GetLeft();
  /** No descriptions */
  virtual int GetTop();
  /** No descriptions */
  virtual int GetBinY();
  /** No descriptions */
  virtual int GetBinX();
  /** Aborts the current exposure */
  virtual void KillExposure();
  /** Shuts down the camera */
  virtual int ShutDownCamera();
  /** Captures an image given the input parameters:
	*myImage - a pointer to the CImage structure in which the image data will be saved
	ExpTime - number of seconds for the current exposure
	gain - the value of the gain to be used in the image
	binx, biny - the x- and y-binning factors
	offsetx, offsety - the x- and y-offset factors for the imaging region
	sizex, sizey - the number of pixels in x,y to be exposed
	dark - flag to decide if a dark image is to be taken or not */
  virtual bool CaptureImage(Image *myImage, string filename, bool singleimage);
  /** Returns the current CCD temperature */
  virtual double GetCCDTemp();
  /** No descriptions */
  virtual double GetCCDSetTemp();
  /** No descriptions */
  virtual void SetDark(bool dark);
  /** No descriptions */
  virtual float GetExposureTime();
  /** No descriptions */
  virtual void ReadTemperatureStatus();
  /** No descriptions */
  virtual void WriteXML(ofstream &os);
  /** Returns the camera handle. */
  virtual short GetCamHandle();
  /** Sets the myDarkData variable. */
  virtual void SetDarkData(IMAGETYPE *data, int width, int height);
  /** Returns the pointer to myDarkData. */
  virtual IMAGETYPE* GetDarkDataPointer();
  /** Sets the gain for the SBIG USB device. */
  virtual void SetGain(float gain);
  /** Gets the USB gain for the SBIG camera. */
  virtual float GetGain();
  /** Turns the CCD temperature regulation off. */
  virtual void SetTemperatureRegulationOff();
  /** Turns the temperature regulation of the CCD on. */
  virtual void SetTemperatureRegulationOn();
  /** No descriptions */
  virtual void SetAutoTemp(bool autoReg);
  /** No descriptions */
  virtual bool GetAutoTemp();
private: // Private methods
  /** Actually performs the exposure */
  virtual bool DoExposureCollect();
private: // Private attributes
  /**  */
  double myCCDSetTemp;
  /**  */
  int myBinX;
  int myBinY;
  int myTop;
  int myLeft;
  int myWidth;
  int myHeight;
  /**  */
  bool myDark;
  /**  */
  float myExposureTime;
};

#endif
