/***************************************************************************
                          cameraapogee  -  description
                             -------------------
    begin                : Wed 11 May 2011
    copyright            : (C) 2011 by Jonathan Makela
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
#ifndef CAMERAAPOGEE_H
#define CAMERAAPOGEE_H

#include <camera.h>
#include "tcl.h"
#include "CameraIO_Linux.h"
#include "ccd.h"

/**
	@author Jonathan Makela <jmakela@uiuc.edu>
*/
class CameraAPOGEE : public Camera
{
public:
    CameraAPOGEE();

    ~CameraAPOGEE();
    bool GetTemperatureRegulation();
    double GetCCDTemp();
    int GetExposureTimeLeft();
    void Read(ifstream &is);
    void Write(ofstream &os);
    void WriteXML(ofstream &os);
    void SetImageSize(int width, int height);
    void SetBinnint(int X, int Y);
    void SetOffsets(int top, int left);
    void SetExpTime(float exp);
    void SetBinX(int X);
    void SetBinY(int Y);
    void SetHeight(int height);
    void SetWidth(int width);
    void SetTop(int top);
    void SetLeft(int left);
    int GetWidth();
    int GetHeight();
    int GetTop();
    int GetLeft();
    int GetBinX();
    int GetBinY();
    double GetCCDSetTemp();
    float GetExposureTime();
    void SetDark(bool dark);
    bool CaptureImage(Image *myImage, string filename, bool singleimage);
    IMAGETYPE * GetDarkDataPointer();
    void SetGain(float gain);
    float GetGain();
    void SetTemperatureRegulationOn();
    void SetTemperatureRegulationOff();
    bool GetAutoTemp();
    void SetAutoTemp(bool autReg);

private:
    bool InitCamera();
    void ReadTemperatureStatus();
    void KillExposure();
    int ShutDownCamera();
    void SetCCDTemp(double CCDTemp);
    bool DoExposureCollect();
    void SetDarkData(IMAGETYPE *data, int width, int height);
private:
    CCameraIO *cameraDevice;
    long myBottom;
    long myLeft;
    long myTop;
    long myRight;
    double myCCDSetTemp;
    int myBinX;
    int myBinY;
    long myWidth;
    long myHeight;
    bool myDark;
    int myExposureTimeLeft;
    int myCols;
    int myRows;
    IMAGETYPE *myCameraData;
    long myLeftOffset;
    long myTopOffset;
public:
    float myExposureTime;
    bool isTempRegulating;
    
protected:
    int attribute_2;
    int myCount;
	long config_load(char *cfgname, short BaseAddress, short RegOffset);
};

#endif
