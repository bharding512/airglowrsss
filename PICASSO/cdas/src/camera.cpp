/***************************************************************************
                          camera.cpp  -  description
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

#include "camera.h"
#include "system.h"

extern System mySystem;

Camera::Camera(){
  // Initialze variables
/*  myCameraOK = false;
  myCCDSetTemp = 0;
  myBinX = 1;
  myBinY = 1;
  myTop = 0;
  myLeft = 0;
  myWidth = 1;
  myHeight = 1;
  myExposureTime = 0.;
  myDarkData = NULL;
  myUSBgain = 1.2;
  isTempRegulating = false;*/
}

Camera::~Camera(){
/*  if(myDarkData != NULL)
    free(myDarkData);
    
  // Check if the system thinks it is exposing
  if(mySystem.myIsExposing)
    KillExposure();

  // Shut down the camera
  ShutDownCamera();*/
}

/** Writes the camera information to the configuration file */
void Camera::Write(ofstream& os){
  // Put the header label
  os << "#Camera" << endl;

  // Write the CCD settemp
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

/** Reads camera information from the configuration file */
void Camera::Read(ifstream& is){
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
  is >> myWidth;
  is >> myTop;
  is >> myLeft;
  is >> myExposureTime;
  is >> myDark;
  is >> myUSBgain;
  is >> myAutoTempRegulation;
}

/** Sets the setpoint for the CCD */
void Camera::SetCCDTemp(double CCDTemp){
  myCCDSetTemp = CCDTemp;
}

/** Gets the flag on temperature regulation */
bool Camera::GetTemperatureRegulation(){
  return isTempRegulating;
}

/** Returns the estimated amount of time before the current exposure is done */
int Camera::GetExposureTimeLeft(){
  return (int) 0;
}

/** Aborts the current exposure */
void Camera::KillExposure(){
}

/** Returns the current CCD temperature */
double Camera::GetCCDTemp(){
  return -1;
}

/** Captures an image given the input parameters:
	*myImage - a pointer to the CImage structure in which the image data will be saved
	ExpTime - number of seconds for the current exposure
	gain - the value of the gain to be used in the image
	binx, biny - the x- and y-binning factors
	offsetx, offsety - the x- and y-offset factors for the imaging region
	sizex, sizey - the number of pixels in x,y to be exposed
	dark - flag to decide if a dark image is to be taken or not */
bool Camera::CaptureImage(Image *myImage, string filename, bool singleimage){
  return false;
}

/** Shuts down the camera */
int Camera::ShutDownCamera(){
  return -1;
}

/** Initializes the camera */
bool Camera::InitCamera(){
  return false;
}

/** Actually performs the exposure */
bool Camera::DoExposureCollect(){
  return false;
}

/** Returns the ambient temperature on the CCD chip */
double Camera::GetAmbientTemp(){
  return -999.;
}

/** Gets the fan power on the CCD */
int Camera::GetFanPower(){
  return -999;
}
/** No descriptions */
void Camera::SetBinning(int X, int Y){
  myBinX = X;
  myBinY = Y;
}

/** No descriptions */
void Camera::SetImagingSize(int width, int height)
{
  myWidth = width;
  myHeight = height;
}

/** No descriptions */
void Camera::SetOffsets(int top, int left){
  myTop = top;
  myLeft = left;
}
/** No descriptions */
void Camera::SetExpTime(float exp){
  myExposureTime = exp;
}
/** No descriptions */
void Camera::SetBinX(int X){
  myBinX = X;
  if(myBinX < 1) myBinX = 1;
}
/** No descriptions */
void Camera::SetBinY(int Y){
  myBinY = Y;
  if(myBinY < 1) myBinY = 1;
}
/** No descriptions */
void Camera::SetTop(int top){
  myTop = top;
}
/** No descriptions */
void Camera::SetLeft(int left){
  myLeft = left;
}
/** No descriptions */
void Camera::SetWidth(int width){
  myWidth = width;
}
/** No descriptions */
void Camera::SetHeight(int height){
  myHeight = height;
}
/** No descriptions */
int Camera::GetBinX(){
  return myBinX;
}
/** No descriptions */
int Camera::GetBinY(){
  return myBinY;
}
/** No descriptions */
int Camera::GetTop(){
  return myTop;
}
/** No descriptions */
int Camera::GetLeft(){
  return myLeft;
}
/** No descriptions */
int Camera::GetHeight(){
  return myHeight;
}
/** No descriptions */
int Camera::GetWidth(){
  return myWidth;
}
/** No descriptions */
double Camera::GetCCDSetTemp(){
  return myCCDSetTemp;
}
/** No descriptions */
void Camera::SetDark(bool dark){
  myDark = dark;
}
/** No descriptions */
float Camera::GetExposureTime(){
  return myExposureTime;
}
/** No descriptions */
void Camera::ReadTemperatureStatus(){
}
/** Write XML information for the camera */
void Camera::WriteXML(ofstream &os){
    // Put the header label
  os << "<camera>" << endl;

  // Write the CCD settemp
  os << "<ccdsettemp>" << myCCDSetTemp << "</ccdsettemp>" << endl;
  os << "<cameratype>" << mySystem.myConfiguration.myCameraType << "</cameratype>" << endl;
  os << "<binx>" << myBinX << "</binx>" << endl;
  os << "<biny>" << myBinY << "</biny>" << endl;
  os << "<imageheight>" << myHeight << "</imageheight>" << endl;
  os << "<imagewidth>" << myWidth << "</imagewidth>" << endl;
  os << "<imagetop>" << myTop << "</imagetop>" << endl;
  os << "<imageleft>" << myLeft << "</imageleft>" << endl;
  os << "<exposuretime> " << myExposureTime << "</exposuretime>" << endl;
  os << "<dark>" << myDark << "</dark>" << endl;
  os << "<usbgain>" << myUSBgain << "</usbgain>" << endl;
  os << "<regulation>" << myAutoTempRegulation << "</regulation>" << endl;

  os << "</camera>" << endl;
}
/** Returns the camera handle. */
short Camera::GetCamHandle(){
  return myCameraHandle;
}
/** Sets the myDarkData variable. */
void Camera::SetDarkData(IMAGETYPE *data, int width, int height){
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
IMAGETYPE* Camera::GetDarkDataPointer(){
  return myDarkData;
}
/** Sets the gain for the SBIG USB device. */
void Camera::SetGain(float gain){
  if(gain < 1)
    gain = 1;
  if(gain > 6)
    gain = 6;
    
  myUSBgain = gain;
}
/** Gets the USB gain for the SBIG camera. */
float Camera::GetGain(){
  return myUSBgain;
}
/** Turns the temperature regulation of the CCD on. */
void Camera::SetTemperatureRegulationOn(){
}
/** Turns the CCD temperature regulation off. */
void Camera::SetTemperatureRegulationOff(){
}
/** No descriptions */
void Camera::SetAutoTemp(bool autoReg){
	myAutoTempRegulation = autoReg;
}
/** No descriptions */
bool Camera::GetAutoTemp(){
	return myAutoTempRegulation;
}
