/***************************************************************************
                          camerapvcam  -  description
                             -------------------
    begin                : Fri Nov 21 2014
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
#include "camerapvcam.h"
#include "system.h"
#include "master.h"
#include "pvcam.h"
#include <sstream>

extern System mySystem;

CameraPVCAM::CameraPVCAM()
{
	// Initialize variables
	myCameraOK = false;
	myLibInitOK = false;
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
	myLastError = TRUE;
	myDark = false;
	
	if(!InitCamera())
	{
		// There was an error
		mySystem.myConfiguration.myCameraType = ERROR;
	}
	else
		mySystem.myConfiguration.myCameraType = PVCAM;
}

CameraPVCAM::~CameraPVCAM()
{
	if(myDarkData != NULL)
		free(myDarkData);
	
	// Destroy the CameraPVCAM object
	if(mySystem.myIsExposing)
		KillExposure();
	
	ShutDownCamera();
	mySystem.myConfiguration.myCameraType = ERROR;
}

/*!
    \fn CameraPVCAM::InitCamera()
 */
bool CameraPVCAM::InitCamera()
{
//	mySystem.myConfiguration.myLogger.WriteSystemLogString("DEBUG: InitCamera");
  
        if(myLibInitOK)
                // we've already initialized the camera libraries
                return true;

        char CameraName[CAM_NAME_LEN];
        short TotalCameras;

        if(!myCameraOK)
        {
                // Good, the camera isn't initialized so we should do it now
                pl_pvcam_init();
                pl_cam_get_total(&TotalCameras);
                pl_cam_get_name(0, CameraName); // We're only concerned with the camera in the first slot
                pl_cam_open(CameraName, &myCameraHandle, OPEN_EXCLUSIVE);
                int n = pl_error_code();

                if(n != 0)
                {
                        // We have an error, close the camera and try again
                        pl_cam_close(myCameraHandle);
                        pl_pvcam_uninit();
                        pl_pvcam_init();
                        pl_cam_get_total(&TotalCameras);
                        pl_cam_get_name(0, CameraName);
                        pl_cam_open(CameraName, &myCameraHandle, OPEN_EXCLUSIVE);

                        if(int n = pl_error_code() != 0)
                        {
                                // Still can't get the camera, must be an error
                                char msg[512];
                                pl_error_message(pl_error_code(), msg);
				  string str("ERROR: Cannot find Camera in CameraPVCAM");
				  mySystem.myConfiguration.myLogger.WriteErrorString(str);
                                pl_exp_abort(myCameraHandle, CCS_CLEAR_CLOSE_SHTR);
                                pl_cam_close(myCameraHandle);
                                pl_pvcam_uninit();
                                return false;
                        }
                }

                // OK, all seems to be well, set the desired temperature and other defaults
                pl_ccd_set_tmp_setpoint(myCameraHandle, (uns16) myCCDSetTemp*100);
                pl_ccd_set_clear_cycles(myCameraHandle, (uns16) 2);
                pl_ccd_set_clear_mode(myCameraHandle, CLEAR_PRE_EXPOSURE);

                myLibInitOK = true;
                myCameraOK = true;
        }

	
	return myCameraOK;
}


/*!
    \fn CameraPVCAM::KillExposure()
 */
void CameraPVCAM::KillExposure()
{
	myLastError = pl_exp_abort(myCameraHandle, CCS_HALT_CLOSE_SHTR);
	if(myLastError != TRUE)
	{
		string str("Cannot kill exposure.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
	}
}


/*!
    \fn CameraPVCAM::ShutDownCamera()
 */
int CameraPVCAM::ShutDownCamera()
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
	
	if(myCameraOK)
        {
                // the camera is open, we need to close it
                pl_cam_close(myCameraHandle);
                pl_pvcam_uninit();
        }

        myCameraOK = false;
        myLibInitOK = false;
	
	return 0;
}


/*!
    \fn CameraPVCAM::ReadTemperatureStatus()
 */
void CameraPVCAM::ReadTemperatureStatus()
{
	int16 TemperatureVal = 0;

        if(myCameraOK)
                pl_ccd_get_tmp(myCameraHandle, &TemperatureVal);
        else
        {
                // the camera isn't open, so open it
                if(!InitCamera())
                        
                pl_ccd_get_tmp(myCameraHandle, &TemperatureVal);
               
        }

		
	myCCDTemp = (double) TemperatureVal/100.;
	
	return;
}


/*!
    \fn CameraPVCAM::SetCCDTemp(int CCDTemp)
 */
void CameraPVCAM::SetCCDTemp(double CCDTemp)
{
    // Place some error bounds on it
	if(CCDTemp > 20.)
		CCDTemp = 20.;
	
	if(myCameraOK)
	{
		pl_ccd_set_tmp_setpoint(myCameraHandle, (uns16) myCCDSetTemp * 100);
		
		myCCDSetTemp = (int) CCDTemp;
		isTempRegulating = true;
	}
}


/*!
    \fn CameraPVCAM::GetTemperatureRegulation()
 */
bool CameraPVCAM::GetTemperatureRegulation()
{
	return isTempRegulating;
}


/*!
    \fn CameraPVCAM::GetCCDTemp()
 */
double CameraPVCAM::GetCCDTemp()
{
	return (double) myCCDTemp;
}


/*!
    \fn CameraPVCAM::GetExposureTimeLeft()
 */
int CameraPVCAM::GetExposureTimeLeft()
{
	return myExposureTimeLeft;
}


/*!
    \fn CameraPVCAM::Read(ifstream &is)
 */
void CameraPVCAM::Read(ifstream &is)
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
	is >> myAutoTempRegulation;
}


/*!
    \fn CameraPVCAM::Write(ofstream &os)
 */
void CameraPVCAM::Write(ofstream &os)
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
    \fn CameraPVCAM::WriteXML(ofstream &os)
 */
void CameraPVCAM::WriteXML(ofstream &os)
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
void CameraPVCAM::SetImageSize(int width, int height){
	myWidth = width;
	myHeight = height;
}

/** Set the binning in the x and y directions */
void CameraPVCAM::SetBinning(int X, int Y){
	myBinX = X;
	myBinY = Y;
}

/** Set the offsets from the top-left corner of the CCD */
void CameraPVCAM::SetOffsets(int top, int left){
	myTop = top;
	myLeft = left;
}

/** Sets the exposure time */
void CameraPVCAM::SetExpTime(float exp){
	myExposureTime = exp;
}

/** Sets the binning in the x-direction only */
void CameraPVCAM::SetBinX(int X){
	myBinX = X;
}

/** Sets binning in the y-direction only */
void CameraPVCAM::SetBinY(int Y){
	myBinY = Y;
}

/** Sets the width of the image */
void CameraPVCAM::SetWidth(int width){
	myWidth = width;
}

/** Sets the height of the image */
void CameraPVCAM::SetHeight(int height){
	myHeight = height;
//	mySystem.myConfiguration.myLogger.WriteSystemLogString("Set height in (SetHeight)");
}

/** Sets the left offset of the CCD */
void CameraPVCAM::SetLeft(int left){
	myLeft = left;
}

/** Sets the top offset of the CCD */
void CameraPVCAM::SetTop(int top){
	myTop = top;
}

/** Returns the width of the image */
int CameraPVCAM::GetWidth(){
	return myWidth;
}

/** Returns the height of the image */
int CameraPVCAM::GetHeight(){
	return myHeight;
}

/** Returns the left offset of the image */
int CameraPVCAM::GetLeft(){
	return myLeft;
}

/** Returns the top offset of the image */
int CameraPVCAM::GetTop(){
	return myTop;
}

/** Returns the binning in the x direction */
int CameraPVCAM::GetBinX(){
	return myBinX;
}

/** Returns the binning in the y direction */
int CameraPVCAM::GetBinY(){
	return myBinY;
}

/** Returns the CCD temperature setpoint */
double CameraPVCAM::GetCCDSetTemp(){
	return myCCDSetTemp;
}

/** Sets the dark image flag */
void CameraPVCAM::SetDark(bool dark){
	myDark = dark;
}

/** Returns the current exposure time */
float CameraPVCAM::GetExposureTime(){
	return myExposureTime;
}


/*!
    \fn CameraPVCAM::CaptureImage(Image *myImage, string filename, bool singleimage)
 */
bool CameraPVCAM::CaptureImage(Image *myImage, string filename, bool singleimage)
{
        unsigned long size;
	
        if(!myCameraOK)
        {
                // the camera isn't ready, this is a problem
                // try initializing the camera
                if(!InitCamera())
                {
			// There was an error
			string str("Error initializing PVCAM camera.");
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
                // We need to take a dark image so tell the shutter to stay closed
                pl_shtr_set_open_mode(myCameraHandle, OPEN_NEVER);
	}
	else
	{
                // We need to open the shutter
                pl_shtr_set_open_mode(myCameraHandle, OPEN_PRE_EXPOSURE);
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
    \fn CameraPVCAM::DoExposureCollect()
 */
bool CameraPVCAM::DoExposureCollect()
{
          // Set up the imaging region
          unsigned long size;
        rgn_type Region;
        Region.s1 = myLeft;
        Region.s2 = myLeft + myWidth - 1;
        Region.p1 = myTop;
        Region.p2 = myTop + myHeight - 1;
        Region.sbin = myBinX;
        Region.pbin = myBinY;
        size = myRows*myCols*2L;

        // initialize the sequence
        pl_exp_init_seq();
        myLastError = pl_exp_setup_seq(myCameraHandle, 1, 1, &Region, TIMED_MODE, (uns32) myExposureTime*1000, &size);
	
	if(myLastError != TRUE)
	{
		// There was an error
		string str("Error in image parameters.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// End and exposure that might have been in progress
	myLastError = pl_exp_abort(myCameraHandle, CCS_HALT_CLOSE_SHTR);
	if(myLastError != TRUE)
	{
		// There was an error
		string str("Error ending exposure.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		return false;
	}
	
	// Start the sequence that has been setup in InitCamera()
        myLastError = pl_exp_start_seq(myCameraHandle, myCameraData);
	if(myLastError != TRUE)
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
	short status;
	unsigned long bytes;
	do
	{
		time(&position);
		myExposureTimeLeft = (int) (stop-position);
		myLastError = pl_exp_check_status(myCameraHandle, &status, &bytes);
	}while(status!= READOUT_COMPLETE);
	
	if(myLastError != TRUE)
	{
		// There was an error
		string str("Error in the image acquisition.");
		mySystem.myConfiguration.myLogger.WriteErrorString(str);
		myCameraOK = false;
		return false;
	}
	  
	// We are done
	pl_exp_finish_seq(myCameraHandle, myCameraData, 0);
        pl_exp_uninit_seq();
	
	// Calculate image stats
	double sum1, sum2;
	sum1 = sum2 = 0.;
	for(int row = 0; row < myRows; row++)
	{
		for(int col = 0; col < myCols; col++)
		{
			sum1 += (double) myCameraData[row*myCols + col];
			sum2 += (double) (((double)myCameraData[row*myCols + col]) * myCameraData[row*myCols + col]);
		}
	}
		
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
    \fn CameraPVCAM::SetDarkData(IMAGETYPE * data, int width, int height)
 */
void CameraPVCAM::SetDarkData(IMAGETYPE * data, int width, int height)
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
    \fn CameraPVCAM::GetDarkDataPointer()
 */
IMAGETYPE * CameraPVCAM::GetDarkDataPointer()
{
	return myDarkData;
}

/*!
    \fn CameraPVCAM::SetGain(float gain)
 */
void CameraPVCAM::SetGain(float gain)
{
	// Note that this is not used in the ANDOR camera but is here to keep everything working
	if(gain < 1)
		gain = 1;
	if(gain > 6)
		gain = 6;
	
	myUSBgain = gain;
}


/*!
    \fn CameraPVCAM::GetGain()
 */
float CameraPVCAM::GetGain()
{
	return myUSBgain;
}


/*!
    \fn CameraPVCAM::SetTemperatureRegulationOn()
 */
void CameraPVCAM::SetTemperatureRegulationOn()
{
	if(myCameraOK)
	{
		myLastError = pl_ccd_set_tmp_setpoint(myCameraHandle,(uns16) myCCDSetTemp*100);
		if(myLastError != TRUE)
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
    \fn CameraPVCAM::SetTemperatureRegulationOff()
 */
void CameraPVCAM::SetTemperatureRegulationOff()
{
	if(myCameraOK)
	{
		myLastError = pl_ccd_set_tmp_setpoint(myCameraHandle,25*100.);
		if(myLastError != TRUE)
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
    \fn CameraPVCAM::SetAutoTemp(bool autoReg)
 */
void CameraPVCAM::SetAutoTemp(bool autoReg)
{
	myAutoTempRegulation = autoReg;
}


/*!
    \fn CameraPVCAM::GetAutoTemp()
 */
bool CameraPVCAM::GetAutoTemp()
{
	return myAutoTempRegulation;
}
