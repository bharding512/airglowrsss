/***************************************************************************
                          camerapvcam  -  description
                             -------------------
    begin                : Fri Nov 21, 2014
    copyright            : (C) 2014 by Jonathan Makela
    email                : jmakela@illinois.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#ifndef CAMERAPVCAM_H
#define CAMERAPVCAM_H

#include <camera.h>

#include "atmcdLXd.h"

/**
	@author Jonathan Makela <jmakela@illinois.edu>
*/
class CameraPVCAM : public Camera
{
public:
    CameraPVCAM();

    ~CameraPVCAM();
    bool GetTemperatureRegulation();
    double GetCCDTemp();
    int GetExposureTimeLeft();
    void Read(ifstream &is);
    void Write(ofstream &os);
    void WriteXML(ofstream &os);

private:
    bool InitCamera();
    void KillExposure();
    int ShutDownCamera();
    void ReadTemperatureStatus();
    void SetCCDTemp(double CCDTemp);
    bool DoExposureCollect();
    unsigned long myLastError;
    int myBinX;
    int myBinY;
    int myTop;
    int myCCDSetTemp;
    int myLeft;
    int myHeight;
    int myWidth;
    float myExposureTime;
    double myCCDTemp;
    int myExposureTimeLeft;
    bool myDark;
    float myUSBgain;
    bool myAutoTempRegulation;
    IMAGETYPE *myCameraData;
    short myCameraHandle;
protected:
    bool myLibInitOK;
public:
    bool myCameraOK;
    bool isTempRegulating;
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
    bool CaptureImage(Image *myImage, string filename, bool singleimage);
    void SetDarkData(IMAGETYPE * data, int width, int height);
    IMAGETYPE * GetDarkDataPointer();
    void SetGain(float gain);
    float GetGain();
    void SetTemperatureRegulationOn();
    void SetTemperatureRegulationOff();
    void SetAutoTemp(bool autoReg);
    bool GetAutoTemp();
    int myRows;
    int myCols;
};

#endif
