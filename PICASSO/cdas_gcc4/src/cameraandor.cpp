/***************************************************************************
                          cameraandor  -  description
                             -------------------
    begin                : Thu Aug 2 2007
    copyright            : (C) 2007 by Jonathan Makela
    email                : jmakela@uiuc.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "cameraandor.h"
#include "system.h"
#include <sstream>

extern System mySystem;

CameraAndor::CameraAndor()
{
	// Initialize variables
	myCameraOK = false;
	myCCDSetTemp = 0;
	myCCDTemp = -999;
	myBinX = 1;
	myBinY = 1;
	myTop = 1;
	myLeft = 1;
	myWidth = 1;
	myHeight = 1;
	myExposureTime = 0.1;
	myDarkData = NULL;
	myAutoTempRegulation = false;
	myUSBgain = 1.0;
	myLastError = DRV_SUCCESS;
	myDark = false;
	
	if(!InitCamera())
	{
		// There was an error
		mySystem.myConfiguration.myCameraType = ERROR;
	}
	else
		mySystem.myConfiguration.myCameraType = ANDOR;
}

CameraAndor::~CameraAndor()
{
	if(myDarkData != NULL)
		free(myDarkData);
	
	// Destroy the CameraAndor object
	if(mySystem.myIsExposing)
		KillExposure();
	
	ShutDownCamera();
	mySystem.myConfiguration.myCameraType = ERROR;
}

/*!
    \fn CameraAndor::InitCamera()
 */
bool CameraAndor::InitCamera()
{
// 	at_32 lNumCameras;
// 	GetAvailableCameras(&lNumCameras);
// 	int iSelectedCamera = 0;
// 	
// 	if(iSelectedCamera < lNumCameras && iSelectedCamera >= 0)
// 	{
// 		GetCameraHandle(iSelectedCamera, &myCameraHandle);
// 		SetCurrentCamera(myCameraHandle);
// 	}
	
	myLastError = Initialize("/usr/local/etc/andor");
	if(myLastError != DRV_SUCCESS)
	{
		// There was an error
		string str("Cannot open ANDOR camera.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
		
	// Sleep to allow initialization to complete
	sleep(2);
	
	// Set Read Mode to IMAGE
	SetReadMode(4);
	
	// Set Acquisition Mode to SINGLE SCAN
	SetAcquisitionMode(1);
	
	// Set Initial Exposure Time
	SetExposureTime(myExposureTime);
	
	// Get Detector Dimensions
	GetDetector(&myWidth, &myHeight);
	
	// Initialize Shutter
	// Output TTL high (1)
	// Shutter mode Auto (0)
	// Closing time 10 ms
	// Opening time 10 ms
	SetShutter(1, 0, 10, 10);
	
	// Setup Image dimensions to something safe
	SetImage(myBinX, myBinY, myLeft, myWidth, myTop, myHeight);
	
	// We have successfully setup the camera
	myCameraOK = true;
	
	// Set CCD temperature
	SetCCDTemp(myCCDSetTemp);
	
	// Read the initial temperature
	ReadTemperatureStatus();

//	SetCoolerMode(0);
	
	// Set readout speeds
	SetHSSpeed(0,2);
	SetVSSpeed(2);
	
	// Set gain
	SetPreAmpGain(2);

	return myCameraOK;
}


/*!
    \fn CameraAndor::KillExposure()
 */
void CameraAndor::KillExposure()
{
	myLastError = AbortAcquisition();
	if(myLastError != DRV_SUCCESS && myLastError != DRV_IDLE)
	{
		string str("Cannot kill exposure.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
	}
}


/*!
    \fn CameraAndor::ShutDownCamera()
 */
int CameraAndor::ShutDownCamera()
{
	// Turn the cooler off
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Warming up sensor to 0 C");
	SetTemperatureRegulationOff();
	while(myCCDTemp < 0)
	{
		ReadTemperatureStatus();
		cout << myCCDTemp << endl;
		sleep(10);
	}
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Sensor warmed up successfully");
	
	SetCoolerMode(0);
	
	// Monitor the temperature until it is > 0 C
//	ReadTemperatureStatus();
	// Shute down the camera
	ShutDown();
	
	myCameraOK = false;
	
	return 0;
}


/*!
    \fn CameraAndor::ReadTemperatureStatus()
 */
void CameraAndor::ReadTemperatureStatus()
{
	int temp;
	
	myLastError = GetStatus(&temp);
	if(temp==DRV_ACQUIRING)
		// We are currently exposing.  Do not disturb the CCD!
		return;
		
	GetTemperature(&temp);
		
	myCCDTemp = (double) temp;
	
	return;
}


/*!
    \fn CameraAndor::SetCCDTemp(int CCDTemp)
 */
void CameraAndor::SetCCDTemp(double CCDTemp)
{
	// Get the range of valid values
	int minT, maxT;
	myLastError = GetTemperatureRange(&minT, &maxT);
	if(myLastError != DRV_SUCCESS)
	{
	        string str("Error getting temperature range for ANDOR camera");
                mySystem.myConfiguration.myLogger.WriteErrorString(str);
                myCameraOK = false;
		return;
	}

    // Place some error bounds on it
	if(CCDTemp > maxT)
		CCDTemp = maxT;
	
	if(myCameraOK)
	{
		// Set the temperature setting
		myLastError = SetTemperature((int) CCDTemp);
		if(myLastError != DRV_SUCCESS)
		{
			string str("Error setting CCD temperature setpoint for ANDOR camera.");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			myCameraOK = false;
			return;
		}
	
		// Turn the cooler on
		myLastError = CoolerON();
		if(myLastError != DRV_SUCCESS)
		{
			string str("Error turning ANDOR cooler on.");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			myCameraOK = false;
			return;
		}
		
		myCCDSetTemp = (int) CCDTemp;
		isTempRegulating = true;
	}
}


/*!
    \fn CameraAndor::GetTemperatureRegulation()
 */
bool CameraAndor::GetTemperatureRegulation()
{
	return isTempRegulating;
}


/*!
    \fn CameraAndor::GetCCDTemp()
 */
double CameraAndor::GetCCDTemp()
{
	return (double) myCCDTemp;
}


/*!
    \fn CameraAndor::GetExposureTimeLeft()
 */
int CameraAndor::GetExposureTimeLeft()
{
	return myExposureTimeLeft;
}


/*!
    \fn CameraAndor::Read(ifstream &is)
 */
void CameraAndor::Read(ifstream &is)
{
	// Reads the camera information from the configuration file
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


/*!
    \fn CameraAndor::Write(ofstream &os)
 */
void CameraAndor::Write(ofstream &os)
{
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
	os << "0" << endl;
	os << "0" << endl;
	os << myUSBgain << endl;
	os << myAutoTempRegulation << endl;
}


/*!
    \fn CameraAndor::WriteXML(ofstream &os)
 */
void CameraAndor::WriteXML(ofstream &os)
{
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
void CameraAndor::SetImageSize(int width, int height){
	myWidth = width;
	myHeight = height;
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (SetImageSize)");
}

/** Set the binning in the x and y directions */
void CameraAndor::SetBinning(int X, int Y){
	myBinX = X;
	myBinY = Y;
}

/** Set the offsets from the top-left corner of the CCD */
void CameraAndor::SetOffsets(int top, int left){
	myTop = top;
	myLeft = left;
}

/** Sets the exposure time */
void CameraAndor::SetExpTime(float exp){
	myExposureTime = exp;
}

/** Sets the binning in the x-direction only */
void CameraAndor::SetBinX(int X){
	myBinX = X;
}

/** Sets binning in the y-direction only */
void CameraAndor::SetBinY(int Y){
	myBinY = Y;
}

/** Sets the width of the image */
void CameraAndor::SetWidth(int width){
	myWidth = width;
}

/** Sets the height of the image */
void CameraAndor::SetHeight(int height){
	myHeight = height;
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (SetHeight)");
}

/** Sets the left offset of the CCD */
void CameraAndor::SetLeft(int left){
	myLeft = left;
}

/** Sets the top offset of the CCD */
void CameraAndor::SetTop(int top){
	myTop = top;
}

/** Returns the width of the image */
int CameraAndor::GetWidth(){
	return myWidth;
}

/** Returns the height of the image */
int CameraAndor::GetHeight(){
	return myHeight;
}

/** Returns the left offset of the image */
int CameraAndor::GetLeft(){
	return myLeft;
}

/** Returns the top offset of the image */
int CameraAndor::GetTop(){
	return myTop;
}

/** Returns the binning in the x direction */
int CameraAndor::GetBinX(){
	return myBinX;
}

/** Returns the binning in the y direction */
int CameraAndor::GetBinY(){
	return myBinY;
}

/** Returns the CCD temperature setpoint */
double CameraAndor::GetCCDSetTemp(){
	return myCCDSetTemp;
}

/** Sets the dark image flag */
void CameraAndor::SetDark(bool dark){
	myDark = dark;
}

/** Returns the current exposure time */
float CameraAndor::GetExposureTime(){
	return myExposureTime;
}


/*!
    \fn CameraAndor::CaptureImage(Image *myImage, string filename, bool singleimage)
 */
bool CameraAndor::CaptureImage(Image *myImage, string filename, bool singleimage)
{
	// First check to make sure that the camera is initialized
	if(!myCameraOK)
	{
		// Camera isn't initialized, so try to initialize it
		if(!InitCamera())
		{
			// There was an error
			string str("Error initializing ANDOR camera.");
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
	
	// Get the pointer to the data
	myCameraData = myImage->GetDataPointer();
	
	if(myDark)
	{
		// We want to perform a dark image
		SetShutter(1,2,10,10);
	}
	else
	{
		// We want to perform a normal image
		SetShutter(1,0,10,10);
	}
	
	// Begin exposing
	mySystem.myIsExposing = true;
	
	// Set the times in the image header
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
		SetDarkData(myCameraData, myRows, myCols);
	
	// Save the heaader information
	myImage->SetGain((int) GetGain());
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
	string tifname;
	tifname = filename + ".tif";
	myImage->SaveTIF(tifname);
	
	return true;
}

/*!
    \fn CameraAndor::DoExposureCollect()
 */
bool CameraAndor::DoExposureCollect()
{
	// Setup Image dimensions
	myLastError = SetImage(myBinX, myBinY, myLeft+1, myWidth, myTop+1, myHeight);
	if(myLastError != DRV_SUCCESS)
	{
		// There was an error
		string str("Error in image parameters.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// End and exposure that might have been in progress
	myLastError = AbortAcquisition();
	if(myLastError != DRV_SUCCESS && myLastError != DRV_IDLE)
	{
		// There was an error
		string str("Error ending exposure.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// Set the exposure time
	myLastError = SetExposureTime(myExposureTime);
	if(myLastError != DRV_SUCCESS)
	{	// There was an error
		string str("Error setting exposure time.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
	
	// Start the exposure
	myLastError = StartAcquisition();
	if(myLastError != DRV_SUCCESS)
	{	// There was an error
		string str("Error starting data acquisition.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
	
	// Calculate the start and stop times
	time_t start, stop, position;
	time(&start);
	stop = (time_t) (start+myExposureTime);
	
	// Monitor the exposure
	int status;
	do
	{
		time(&position);
		myExposureTimeLeft = (int) (stop-position);
		myLastError = GetStatus(&status);
	}while(status==DRV_ACQUIRING);
	
	if(myLastError != DRV_SUCCESS)
	{
		// There was an error
		string str("Error in the image acquisition.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
	
	// Read the CCD
	long datasize = (long) myRows*myCols;
	at_32* imageData = new at_32[datasize];
	myLastError = GetAcquiredData(imageData, datasize);
	if(myLastError != DRV_SUCCESS)
	{
		// There was an error
		string str("Error reading data off of camera.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		delete[] imageData;
		return false;
	}
	
	// Copy the data over into the image data pointer
	double sum1, sum2;
	sum1 = sum2 = 0.;
	for(int row = 0; row < myRows; row++)
	{
		for(int col = 0; col < myCols; col++)
		{
			myCameraData[row*myCols + col] = (IMAGETYPE) imageData[row*myCols + col];
			sum1 += (double) imageData[row*myCols + col];
			sum2 += (double) (((double)imageData[row*myCols + col]) * imageData[row*myCols + col]);
		}
	}
	
	delete[] imageData;
	
	// Calculate the mean and stddev
	double mean;
	mean = sum1 / (myRows * myCols);
	double temp;
	temp = sum2 / (myRows * myCols);
	double stddev;
	if(temp < mean*mean)
		stddev = 0.;
	else
		stddev = sqrt(temp - mean*mean);
	
	mySystem.myConfiguration.myImage->SetMean(mean);
	mySystem.myConfiguration.myImage->SetStdDev(stddev);
	
	return true;
}


/*!
    \fn CameraAndor::SetDarkData(IMAGETYPE * data, int width, int height)
 */
void CameraAndor::SetDarkData(IMAGETYPE * data, int width, int height)
{
	if(myDarkData != NULL)
		free(myDarkData);
	
	// Allocate the array data
	myDarkData = (IMAGETYPE *) new BYTE[sizeof(IMAGETYPE)*width*height];
	
	// Copy over the data
	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
			myDarkData[i*height+j] = data[i*height + j];
}


/*!
    \fn CameraAndor::GetDarkDataPointer()
 */
IMAGETYPE * CameraAndor::GetDarkDataPointer()
{
	return myDarkData;
}

/*!
    \fn CameraAndor::SetGain(float gain)
 */
void CameraAndor::SetGain(float gain)
{
	// Note that this is not used in the ANDOR camera but is here to keep everything working
	if(gain < 1)
		gain = 1;
	if(gain > 6)
		gain = 6;
	
	myUSBgain = gain;
}


/*!
    \fn CameraAndor::GetGain()
 */
float CameraAndor::GetGain()
{
	return myUSBgain;
}


/*!
    \fn CameraAndor::SetTemperatureRegulationOn()
 */
void CameraAndor::SetTemperatureRegulationOn()
{
	if(myCameraOK)
	{
		myLastError = CoolerON();
		if(myLastError != DRV_SUCCESS)
		{
			// There was an error
			string str("Error turning CCD cooler on.");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			return;
		}
	
		mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation ON");
		isTempRegulating = true;
	}
}


/*!
    \fn CameraAndor::SetTemperatureRegulationOff()
 */
void CameraAndor::SetTemperatureRegulationOff()
{
	if(myCameraOK)
	{
		myLastError = CoolerOFF();
		if(myLastError != DRV_SUCCESS)
		{
			// There was an error
			string str("Error turning CCD cooler off.");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			return;
		}
	
		mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation OFF");
		isTempRegulating = false;
	}
}


/*!
    \fn CameraAndor::SetAutoTemp(bool autoReg)
 */
void CameraAndor::SetAutoTemp(bool autoReg)
{
	myAutoTempRegulation = autoReg;
}


/*!
    \fn CameraAndor::GetAutoTemp()
 */
bool CameraAndor::GetAutoTemp()
{
	return myAutoTempRegulation;
}
