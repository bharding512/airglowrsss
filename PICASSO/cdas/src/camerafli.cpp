/***************************************************************************
                          camerafli  -  description
                             -------------------
    begin                : Fri Aug 17 2007
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
#include "camerafli.h"
#include "system.h"
#include <sstream>

extern System mySystem;

		
CameraFLI::CameraFLI()
{
	// Initialize variables
	myCameraOK = false;
	myCCDSetTemp = 0;
	myCCDTemp = 25;
	myBinX = 1;
	myBinY = 1;
	myTop = 1;
	myLeft = 1;
	myLeftOffset = 0;
	myTopOffset = 0;
	myWidth = 1;
	myHeight = 1;
	myExposureTime = 0.1;
	myDarkData = NULL;
	myAutoTempRegulation = false;
	myUSBgain = 1.0;
	myDark = false;
	myLastError = 0;
	
	if(!InitCamera())
		// There was an error
		mySystem.myConfiguration.myCameraType = ERROR;
	else
		// Alls well
		mySystem.myConfiguration.myCameraType = FLIPROLINE;
}


CameraFLI::~CameraFLI()
{
	if(myDarkData != NULL)
		free(myDarkData);
	
	// Stop the exposure, if it is occuring
	if(mySystem.myIsExposing)
		KillExposure();
	
	ShutDownCamera();
	
	mySystem.myConfiguration.myCameraType = ERROR;
}

/*!
    \fn CameraFLI::InitCamera
 */
bool CameraFLI::InitCamera()
{
	int i, j;
	char **list;
	
	/** Poll FLI driver for attached USB camera devices */
	myLastError = FLIList(FLIDOMAIN_USB | FLIDEVICE_CAMERA, &list);
	/** Reformat returned list into C strings with \0 terminators */
	for (i = 0; list[i] != NULL; i++) {
		for (j = 0; list[i][j] != '\0'; j++)
			if (list[i][j] == ';') {
			list[i][j] = '\0';
			break;
			}
	}
	/** Select first operational USB filterwheel for our use */
	for (i = 0; list[i] != NULL; i++) {
		myLastError = FLIOpen(&cameraDevice, list[i], FLIDOMAIN_USB | FLIDEVICE_CAMERA);
		if (myLastError != 0) {
			// There was an error
			string str("Cannot open FLIPROLINE camera");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			myCameraOK = false;
			FLIFreeList(list);
			return false;
		}
	}
	FLIFreeList(list);
	
	// Set up CCD into imaging mode
	myLastError = FLISetFrameType(cameraDevice, FLI_FRAME_TYPE_NORMAL);
	if(myLastError != 0)
	{
		string str("Error in FLIPROLINE");
		str += strerror((int)-myLastError);
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// Flush the CCD 2 times before exposing
	myLastError = FLISetNFlushes(cameraDevice, 2);
	
	// Set Initial exposure time
	myLastError = FLISetExposureTime(cameraDevice, (long) (myExposureTime*1000));
	
	// Get detector dimensions
	myLastError = FLIGetVisibleArea(cameraDevice, &myLeftOffset, &myTopOffset, &myRight, &myBottom);
	myWidth = myRight-myLeftOffset;
	myHeight = myBottom - myTopOffset;
	
	// Set the image dimensions to the total area
	myLastError =FLISetImageArea(cameraDevice, myLeft+myLeftOffset, myTop+myTopOffset, myLeft+myLeftOffset+myWidth/myBinX, myTop+myTopOffset+myHeight/myBinY);
	
	// Set the CCD temperature
	myLastError = FLISetTemperature(cameraDevice, myCCDSetTemp);
	
	// Read the initial temperature
	ReadTemperatureStatus();
	
	// We have sucessfully setup the camera
	myCameraOK = true;
	
	return myCameraOK;
}

/*!
    \fn CameraFLI::ReadTemperatureStatus()
 */
void CameraFLI::ReadTemperatureStatus()
{	
	myLastError = FLIGetTemperature(cameraDevice, &myCCDTemp);
	
	return;
}


/*!
    \fn CameraFLI::KillExposure()
 */
void CameraFLI::KillExposure()
{
	myLastError = FLICancelExposure(cameraDevice);
	if(myLastError != 0)
	{
		string str("Cannot kill exposure");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
	}
}


/*!
    \fn CameraFLI::ShutDownCamera
 */
int CameraFLI::ShutDownCamera()
{
	// Warm the CCD up to 0
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Warming up sensor to 0 C");
	FLISetTemperature(cameraDevice, 20.0);
	while(myCCDTemp < 0)
	{
		ReadTemperatureStatus();
		cout << myCCDTemp << endl;
		sleep(10);
	}
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Sensor warmed up successfully");
	
	FLIClose(cameraDevice);
	
	myCameraOK = false;
	
	return 0;
}


/*!
    \fn CameraFLI::SetCCDTemp()
 */
void CameraFLI::SetCCDTemp(double CCDTemp)
{
	// Place some error bounds on it
	if(CCDTemp > 20.)
		CCDTemp = 20.;
	
	if(myCameraOK)
	{
		// Set the temperature setting
		myLastError = FLISetTemperature(cameraDevice, CCDTemp);
		if(myLastError != 0)
		{
			string str("Error setting the CCD temperature setpoint of the FLIPROLINE camera");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			myCameraOK = false;
			return;
		}
		
		myCCDSetTemp = CCDTemp;
		isTempRegulating = true;
	}
}


/*!
    \fn CameraFLI::GetTemperatureRegulation()
 */
bool CameraFLI::GetTemperatureRegulation()
{
	return isTempRegulating;
}


/*!
    \fn CameraFLI::GetCCDTemp()
 */
double CameraFLI::GetCCDTemp()
{
	return myCCDTemp;
}


/*!
    \fn CameraFLI::GetExposureTimeLeft()
 */
int CameraFLI::GetExposureTimeLeft()
{
	return myExposureTimeLeft;
}


/*!
    \fn CameraFLI::Read(ifstream &is)
 */
void CameraFLI::Read(ifstream &is)
{
	// Reads the camera information from the configuration file
	// Search for the header from the beginning of the file
	is.seekg(0,ios::beg);
	
	string myLine;
	
	// Read in the next line of the file
	getline(is, myLine);
	
	// Find the header
	while(myLine != "#Caemra" && !is.eof())
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
	is >> myAutoTempRegulation;}


/*!
    \fn CameraFLI::Write(ofstream &os)
 */
void CameraFLI::Write(ofstream &os)
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
    \fn CameraFLI::WriteXML(ofstream &os)
 */
void CameraFLI::WriteXML(ofstream &os)
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


/*!
    \fn CameraFLI::SetImageSize(int width, int height)
 */
void CameraFLI::SetImageSize(int width, int height)
{
	myWidth = (long) width;
	myHeight = (long) height;
}


/*!
    \fn CameraFLI::SetBinnint(int X, int Y)
 */
void CameraFLI::SetBinnint(int X, int Y)
{
	myBinX = X;
	myBinY = Y;
}


/*!
    \fn CameraFLI::SetOffsets(int top, int left)
 */
void CameraFLI::SetOffsets(int top, int left)
{
	myTop = top;
	myLeft = left;
}


/*!
    \fn CameraFLI::SetExpTime(float exp)
 */
void CameraFLI::SetExpTime(float exp)
{
	myExposureTime = exp;
}


/*!
    \fn CameraFLI::SetBinX(int X)
 */
void CameraFLI::SetBinX(int X)
{
	myBinX = X;
}


/*!
    \fn CameraFLI::SetBinY(int Y)
 */
void CameraFLI::SetBinY(int Y)
{
	myBinY = Y;
}


/*!
    \fn CameraFLI::SetHeight(int height)
 */
void CameraFLI::SetHeight(int height)
{
	myHeight = height;
}


/*!
    \fn CameraFLI::SetWidth(int width)
 */
void CameraFLI::SetWidth(int width)
{
	myWidth = width;
}


/*!
    \fn CameraFLI::SetTop(int top)
 */
void CameraFLI::SetTop(int top)
{
	myTop = top;
}


/*!
    \fn CameraFLI::SetLeft(int left);
 */
void CameraFLI::SetLeft(int left)
{
	myLeft = left;
}


/*!
    \fn CameraFLI::GetWidth()
 */
int CameraFLI::GetWidth()
{
	return myWidth;
}


/*!
    \fn CameraFLI::GetHeight()
 */
int CameraFLI::GetHeight()
{
	return myHeight;
}


/*!
    \fn CameraFLI::GetTop()
 */
int CameraFLI::GetTop()
{
	return myTop;
}


/*!
    \fn CameraFLI::GetLeft()
 */
int CameraFLI::GetLeft()
{
	return myLeft;
}


/*!
    \fn CameraFLI::GetBinX();
 */
int CameraFLI::GetBinX()
{
	return myBinX;
}


/*!
    \fn CameraFLI::GetBinY()
 */
int CameraFLI::GetBinY()
{
	return myBinY;
}


/*!
    \fn CameraFLI::GetCCDSetTemp()
 */
double CameraFLI::GetCCDSetTemp()
{
	return myCCDSetTemp;
}


/*!
    \fn CameraFLI::GetExposureTime()
 */
float CameraFLI::GetExposureTime()
{
	return myExposureTime;
}


/*!
    \fn CameraFLI::SetDark(bool dark)
 */
void CameraFLI::SetDark(bool dark)
{
	myDark = dark;
}


/*!
    \fn CameraFLI::CaptureImage(Image *myImage, string filename, bool singleimage)
 */
bool CameraFLI::CaptureImage(Image *myImage, string filename, bool singleimage)
{
	// First check to make sure that the camera is initialized
	if(!myCameraOK)
	{
		// Camera isn't initialized, so try to initialize it
		if(!InitCamera())
		{
			// There was an error
			string str("Error initializing FLIPROLINE camera.");
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
		myLastError = FLISetFrameType(cameraDevice, FLI_FRAME_TYPE_DARK);
	}
	else
	{
		// We want to perform a normal image
		myLastError = FLISetFrameType(cameraDevice, FLI_FRAME_TYPE_NORMAL);
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
	
	// Save the header information
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
    \fn CameraFLI::DoExposureCollect()
 */
bool CameraFLI::DoExposureCollect()
{
	// Setup Image dimensions
	myLastError = FLISetImageArea(cameraDevice, myLeft+myLeftOffset, myTop+myTopOffset, myLeft+myLeftOffset+myWidth/myBinX, myTop+myTopOffset+myHeight/myBinY);
	if(myLastError != 0)
	{
		// There was an error
		string str("Error in image parameters.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// Set Binning
	myLastError = FLISetHBin(cameraDevice, myBinX);
	if(myLastError != 0)
	{
		// There was an error
		string str("Error in horizontal binning.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	myLastError = FLISetVBin(cameraDevice, myBinY);
	if(myLastError != 0)
	{
		// There was an error
		string str("Error in vertical binning.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// End and exposure that might have been in progress
	myLastError = FLICancelExposure(cameraDevice);
	if(myLastError != 0)
	{
		// There was an error
		string str("Error ending exposure.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// Set the exposure time
	myLastError = FLISetExposureTime(cameraDevice, (long) (myExposureTime*1000));
	if(myLastError != 0)
	{	// There was an error
		string str("Error setting exposure time.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
	
	// Start the exposure
	myLastError = FLIExposeFrame(cameraDevice);
	if(myLastError != 0)
	{	// There was an error
		string str("Error starting data acquisition.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
	
	// Calculate the start and stop times
//	time_t start, stop, position;
//	time(&start);
//	stop = (time_t) (start+myExposureTime);
	
	// Monitor the exposure
	long timeleft;
	do
	{
//		time(&position);
//		myExposureTimeLeft = (int) (stop-position);
		myLastError = FLIGetExposureStatus(cameraDevice, &timeleft);
		myExposureTimeLeft = timeleft*1000;
		sleep(1);
	}while(myExposureTimeLeft > 0);
	
	if(myLastError != 0)
	{
		// There was an error
		string str("Error in the image acquisition.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
	
	// Read the CCD
	u_int16_t *img;
	int row;
	double sum1, sum2;
	sum1 = sum2 = 0.;
	
	// allocate the image
	img = (u_int16_t*) malloc(myRows*myCols*sizeof(u_int16_t));
		
	for(row = 0; row < myRows; row++)
	{
		myLastError = FLIGrabRow(cameraDevice, &img[row*myCols], myCols);
		if(myLastError != 0)
		{
			// There was an error
			string str("Error reading data off of camera.");
			str += strerror((int)-myLastError);
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			myCameraOK = false;
			return false;
		}
		
		// Copy the data over into the image data pointer
		for(int col = 0; col < myCols; col++)
		{
			myCameraData[row*myCols + col] = (IMAGETYPE) img[row*myCols + col];
			sum1 += (double) myCameraData[row*myCols + col];
			sum2 += (double) (((double)myCameraData[row*myCols + col]) * myCameraData[row*myCols + col]);
		}
	}
	
	delete[] img;
	
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
    \fn CameraFLI::SetDarkData(IMAGETYPE *data, int width, int height)
 */
void CameraFLI::SetDarkData(IMAGETYPE *data, int width, int height)
{
	if(myDarkData != NULL)
		free(myDarkData);
	
	// Allocate the array data
	myDarkData = (IMAGETYPE *) new BYTE[sizeof(IMAGETYPE)*width*height];
	
	// Copy over the data
	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
			myDarkData[i*height+j] = data[i*height + j];}


/*!
    \fn CameraFLI::GetDarkDataPointer()
 */
IMAGETYPE * CameraFLI::GetDarkDataPointer()
{
	return myDarkData;
}


/*!
    \fn CameraFLI::SetGain(float gain)
 */
void CameraFLI::SetGain(float gain)
{
	// Note that this is not used in the FLIPROLIE camera, but we keep
	// it here so everyhting works.
	if(gain < 1)
		gain = 1;
	if(gain > 6)
		gain = 6;
	
	myUSBgain = gain;
}


/*!
    \fn CameraFLI::GetGain()
 */
float CameraFLI::GetGain()
{
	return myUSBgain;
}


/*!
    \fn CameraFLI::SetTemperatureRegulationOn()
 */
void CameraFLI::SetTemperatureRegulationOn()
{
	if(myCameraOK)
	{
		myLastError = FLISetTemperature(cameraDevice, myCCDSetTemp);
		if(myLastError != 0)
		{
			string str("Error turning CCD cooler on.");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			return;
		}
		
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation ON");
	isTempRegulating = true;
	}
}


/*!
    \fn CameraFLI::SetTemperatureRegulationOff()
 */
void CameraFLI::SetTemperatureRegulationOff()
{
	if(myCameraOK)
	{
		myLastError =FLISetTemperature(cameraDevice, 25.0);
		if(myLastError != 0)
		{
			string str("Error turning CCD cooler off (25 C).");
			mySystem.myConfiguration.myLogger.WriteErrorString(str);
			return;
		}
		
		mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation off");
		isTempRegulating = false;
	}
}


/*!
    \fn CameraFLI::GetAutoTemp()
 */
bool CameraFLI::GetAutoTemp()
{
	return myAutoTempRegulation;
}


/*!
    \fn CameraFLI::SetAutoTemp(bool auteReg)
 */
void CameraFLI::SetAutoTemp(bool autoReg)
{
	myAutoTempRegulation = autoReg;
}
