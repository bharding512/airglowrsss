/***************************************************************************
                          camerasbig.cpp  -  description
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

#include "camerasbig.h"
#include "system.h"

extern System mySystem;

CameraSBIG::CameraSBIG(){
  // Initialze variables
  myCameraOK = false;
  myCCDSetTemp = 0;
  myCameraHandle = (SbigUsbCam*)NULL;
  myCCDTemp = -999;
  myAmbientTemp = -999;
  myFanPower = -999;
  myBinX = 1;
  myBinY = 1;
  myTop = 0;
  myLeft = 0;
  myWidth = 1;
  myHeight = 1;
  myExposureTime = 0.;

  if(!InitCamera())
  {
    mySystem.myConfiguration.myCameraType = ERROR;
  }
  else mySystem.myConfiguration.myCameraType = SBIGUSB;    
}

CameraSBIG::~CameraSBIG(){
  if(mySystem.myIsExposing)
    KillExposure();

  ShutDownCamera();
  mySystem.myConfiguration.myCameraType = ERROR;
}

/** Captures an image given the input parameters:
	*myImage - a pointer to the CImage structure in which the image data will be saved
	ExpTime - number of seconds for the current exposure */
bool CameraSBIG::CaptureImage(Image * myImage, string filename){
  
  if(!myCameraOK)
  {
    // Camera isn't initialized, so try to initialize it
    if(!InitCamera())
    {
      // The was an error
      string str("Error initializing camera: ");
      str = str + myCameraHandle->GetStatusString();
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

  // Set the times in the image header
  mySystem.myConfiguration.myImage->SetTimesNow();
  // Update the temperature
  mySystem.myConfiguration.myCamera->ReadTemperatureStatus();
   
  // Get the number of columns and rows in our image
  myRows = myHeight / myBinY;
  myCols = myWidth / myBinX;

  // Reflect this in the image structure
  myImage->SetRows(myRows);
  myImage->SetCols(myCols);

  myCameraData = myImage->GetDataPointer();

  if(myDark)
  {
    // We want to perform a dark image
    myGIP.openShutter = 2;
  }
  else
  {
    // We want to perform a normal image
    myGIP.openShutter = 1;
  }

  // Setup the imaging region
  myGIP.top = myTop/myBinY;
  myGIP.left = myLeft/myBinX;
  myGIP.width = myWidth/myBinX;
  myGIP.height = myHeight/myBinY;

  // Set the exposuretime
  myGIP.exposureTime = (float) myExposureTime;

  // As a default use medium resolution mode
  myGIP.readoutMode = 1;
  
  if(myBinX == 1 and myBinY == 1)
  {
    // High resolution mode
    myGIP.readoutMode = 0;
  }
  if(myBinX == 2 and myBinY == 2)
  {
    // Medium resolution mode
    myGIP.readoutMode = 1;
  }
  if(myBinX == 3 and myBinY == 3)
  {
    // Low resolution mode
    myGIP.readoutMode = 2;
  }

  mySystem.myIsExposing = true;

  // Perform the exposure
  if(!DoExposureCollect())
  {
    // Error in the exposure, the error log should have been set in the function
    return false;
  }

  mySystem.myIsExposing = false;

  // Save header information into the image class
  mySystem.myConfiguration.myImage->SetGain((int) (GetGain() * 10));
  mySystem.myConfiguration.myImage->SetCCDTemp((float)myCCDTemp);
  mySystem.myConfiguration.myImage->SetBinning(mySystem.myConfiguration.myCamera->GetBinX(),
                                               mySystem.myConfiguration.myCamera->GetBinY());
  mySystem.myConfiguration.myImage->SetExposureTime(mySystem.myConfiguration.myCamera->GetExposureTime());
  mySystem.myConfiguration.myImage->SetSite(&mySystem.myConfiguration.mySite);
  mySystem.myConfiguration.myImage->AutoContrast();
  mySystem.myConfiguration.myImage->Normalize();

  // Save the file to disk
  ofstream myFile;
  string tifname;
  tifname = filename + ".tif";
  myFile.open(tifname.c_str(), ios::out | ios::binary);
  mySystem.myConfiguration.myImage->SaveTIF(myFile);
  myFile.close();

  return true;
}

/** Actually performs the exposure */
bool CameraSBIG::DoExposureCollect(){
  // Assumes the GrabImageParameters has been initialized properly

  // Check to make sure the camera is still there
  if(!myCameraHandle)
  {
    // The was an error
    string str("Cannot find SBIG USB camera: ");
    str = str + string(myCameraHandle->GetStatusString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check to make sure we're using a valid CCD
  if(myGIP.ccd > CCD_TRACKING)
  {
    // The was an error
    string str("Invalid CCD Parameter");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check the anti-blooming parameter
  if(myGIP.abgState > ABG_CLK_HI7)
  {
    // The was an error
    string str("Invalid anti-blooming Parameter");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check the shutter parameter
  if(myGIP.openShutter > SC_INITIALIZE_SHUTTER)
  {
    // The was an error
    string str("Invalid Shutter Parameter");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check the exposure time
  unsigned long exposureTime = (unsigned long) (myGIP.exposureTime * 100.0 + 0.5);
  if(exposureTime < MIN_ST7_EXPOSURE)
  {
    exposureTime = MIN_ST7_EXPOSURE;
    myExposureTime = ((float)exposureTime - 0.5)/100.;
  }

  // Get the CCD info
  myGCCDIP.request = (myGIP.ccd == CCD_IMAGING) ? 0:1;
  int status;
  if((status = myCameraHandle->GetCCDInfo(&myGCCDIP, (void *)&myGCCDIR0)) != CE_NO_ERROR)
  {
    // The was an error
    string str("Error getting CCD info: ");
    str = str + string(myCameraHandle->GetStatusString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Check readout mode
  short myReadoutMode = myGIP.readoutMode & 0xFF;
  short myYBin = myGIP.readoutMode >> 8;

  if(myGIP.ccd == CCD_TRACKING && myReadoutMode > 1)
  {
    // The was an error
    string str("Invalid readout mode for tacking CCD");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  if(myReadoutMode > 9)
  {
    // The was an error
    string str("Invalid readout mode");
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
      // The was an error
      string str("Invalid Y binning parameter");
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
//  if(myTop > (myCCDHeight - 1))
  if(myGIP.top > (myCCDHeight - 1))
  {
    // The was an error
    string str("Image top off of CCD");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

//  if(myLeft > (myCCDWidth - 1))
  if(myGIP.left > (myCCDWidth - 1))
  {
    // The was an error
    string str("Image left off of CCD");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

//  if(myWidth == 0)
  if(myGIP.width == 0)
  {
    // Use the entire CCD width
    myWidth = myCCDWidth;
    myLeft = 0;
  }
  else
  {
//    if((myLeft + myCols) > myCCDWidth)
    if((myGIP.left + myGIP.width) > myCCDWidth)
    {
      // The was an error
      string str("Image right off of CCD");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

//  if(myHeight == 0)
  if(myGIP.height == 0)
  {
    // Use the entire CCD height
    myHeight = myCCDHeight;
    myTop = 0;
  }
  else
  {
//    if((myTop + myRows) > myCCDHeight)
    if((myGIP.top + myGIP.height) > myCCDHeight)
    {
      // The was an error
      string str("Image bottom off of CCD");
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

  // Fill the StartExposureParams structure
  StartExposureParams mySEP;
  mySEP.ccd = myGIP.ccd;
  mySEP.abgState = myGIP.abgState;
  mySEP.openShutter = myGIP.openShutter;
  mySEP.exposureTime = exposureTime;

  // Start the exposure
  if((status = myCameraHandle->StartExposure(&mySEP)) != CE_NO_ERROR)
  {
    // The was an error
    string str("Error starting exposure: ");
    str = str + myCameraHandle->GetStatusString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);

    // Send an EndExposure
    EndExposureParams myEEP;
    myEEP.ccd = myGIP.ccd;
    myCameraHandle->EndExposure(&myEEP);

    myCameraOK = false;
    return false;
  }

  unsigned long start, end, position;
  // Calculate the start and stop times
  start = myCameraHandle->GetJiffies();
  end = start + exposureTime;

  // Monitor the exposure
  QueryCommandStatusResults myQCSR;
  QueryCommandStatusParams myQCSP;
  myQCSP.command = CC_START_EXPOSURE;
  if(myGIP.ccd == CCD_IMAGING)
  {
    do
    {
      position = myCameraHandle->GetJiffies();
      // Calculate the number of seconds remaining;
      myExposureTimeLeft = (int) ((end - position) * 100.);

      if((status = myCameraHandle->QueryCommandStatus(&myQCSP,&myQCSR)) != CE_NO_ERROR)
      {
        // The was an error
        string str("Error during exposure: ");
        str = str + myCameraHandle->GetStatusString();
        mySystem.myConfiguration.myLogger.WriteErrorString(str);

        // Send an EndExposure
        EndExposureParams myEEP;
        myEEP.ccd = myGIP.ccd;
        myCameraHandle->EndExposure(&myEEP);

        myCameraOK = false;
        return false;
      }
    } while((myQCSR.status & 3) != 3);
  }
  else
  {
    do
    {
      position = myCameraHandle->GetJiffies();
      // Calculate the number of seconds remaining;
      myExposureTimeLeft = (int) ((end - position) * 100.);

      if((status = myCameraHandle->QueryCommandStatus(&myQCSP,&myQCSR)) != CE_NO_ERROR)
      {
        // The was an error
        string str("Error during exposure: ");
        str = str + myCameraHandle->GetStatusString();
        mySystem.myConfiguration.myLogger.WriteErrorString(str);

        // Send an EndExposure
        EndExposureParams myEEP;
        myEEP.ccd = myGIP.ccd;
        myCameraHandle->EndExposure(&myEEP);

        myCameraOK = false;
        return false;
      }
    } while((myQCSR.status & 12) != 12);
  }

  // End exposure
  EndExposureParams myEEP;
  myEEP.ccd = myGIP.ccd;
  if((status = myCameraHandle->EndExposure(&myEEP)) != CE_NO_ERROR)
  {
    // The was an error
    string str("Error ending exposure: ");
    str = str + myCameraHandle->GetStatusString();
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Dump lines if neccessary
//  if(myTop > 0)
  if(myGIP.top > 0)
  {
    DumpLinesParams myDLP;
    myDLP.ccd = myGIP.ccd;
    myDLP.readoutMode = myGIP.readoutMode;
//    myDLP.lineLength = myTop;
    myDLP.lineLength = myGIP.top;
    status = myCameraHandle->DumpLines(&myDLP);
    if(status != CE_NO_ERROR)
    {
      // The was an error
      string str("Error dumping lines off of CCD");
      str = str + myCameraHandle->GetStatusString();
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      myCameraOK = false;
      return false;
    }
  }

  // Start readingout lines
  ReadoutLineParams myRLP;
  myRLP.ccd = myGIP.ccd;
  myRLP.readoutMode = myGIP.readoutMode;
//  myRLP.pixelStart = myLeft;
  myRLP.pixelStart = myGIP.left;
//  myRLP.pixelLength = myCols;
  myRLP.pixelLength = myGIP.width;

  unsigned short myLocalBuffer[MAX_DIG_WIDTH];
  int row, col;

  // Freeze the TE Cooler
  mySTRP.regulation = 3;
  myCameraHandle->SetTemperatureRegulation(&mySTRP);

  // variables for calculating the mean and stddev of the image
  double sum1;
  double sum2;
  sum1 = 0;
  sum2 = 0;
  
  // Begin readout
  if(myGIP.subtract)
  {
    // Performing on-chip subtration
    for(row = 0; row < myRows; row++)
    {
      // Copy the array to the local buffer
      for(col = 0; col < myCols; col ++)
      {
        myLocalBuffer[col] = myCameraData[row*myCols + col];
      }

      // readout CCD line & subtract local data
      if(myCameraHandle->ReadoutLine(&myRLP, myLocalBuffer, 1) != CE_NO_ERROR)
      {
        // The was an error
        string str("Error reading lines off of CCD");
        str = str + myCameraHandle->GetStatusString();
        mySystem.myConfiguration.myLogger.WriteErrorString(str);
        // End the readout but change the fan freeze mode
        mySTRP.regulation = 4;
        myCameraHandle->SetTemperatureRegulation(&mySTRP);
        myCameraOK = false;
        return false;
      }

      for(col = 0; col < myCols; col++)
      {
        myCameraData[row*myCols + col] = (IMAGETYPE) myLocalBuffer[col];
        sum1 += myLocalBuffer[col];
        sum2 += (double) (((double) myLocalBuffer[col]) * myLocalBuffer[col]);
      }
    }
  }
  else
  {
    for(row = 0; row < myRows; row++)
    {
      if((status = myCameraHandle->ReadoutLine(&myRLP, myLocalBuffer, 0)) != CE_NO_ERROR)
      {
        // The was an error
        string str("Error reading lines off of CCD: ");
        str = str + myCameraHandle->GetStatusString();
        mySystem.myConfiguration.myLogger.WriteErrorString(str);
        // End the readout but change the fan freeze mode
        mySTRP.regulation = 4;
        myCameraHandle->SetTemperatureRegulation(&mySTRP);
        myCameraOK = false;
        return false;
      }

      // Copy the data to the data array
      for(col = 0; col < myCols; col++)
      {
        myCameraData[row*myCols + col] = (IMAGETYPE) myLocalBuffer[col];
        sum1 += myLocalBuffer[col];
        sum2 += (double) (((double) myLocalBuffer[col]) * myLocalBuffer[col]);
      }
    }
  }

  // calculate the mean and stddev
  double mean;
  mean = sum1 / (myRows * myCols);
  double temp;
  temp = sum2 / (myRows * myCols);
  double stddev;
  if(temp < mean * mean)
    stddev = 0;
  else
    stddev = sqrt(temp - mean * mean);

  mySystem.myConfiguration.myImage->SetMean(mean);
  mySystem.myConfiguration.myImage->SetStdDev(stddev);
  
  // End the readout but change the fan freeze mode
  mySTRP.regulation = 4;
  myCameraHandle->SetTemperatureRegulation(&mySTRP);
  return true;
}

/** Returns the current CCD temperature */
double CameraSBIG::GetCCDTemp(){
  // Return the current CCD temperature
  return myCCDTemp;
}

/** Returns the estimated amount of time before the current exposure is done */
int CameraSBIG::GetExposureTimeLeft(){
  return myExposureTimeLeft;
}

/** Initializes the camera */
bool CameraSBIG::InitCamera(){
  // Detect the camera.
  // Currently assume the name is sbigusb0
  
  if(myCameraHandle) delete myCameraHandle;
  myCameraHandle = (SbigUsbCam*)NULL;
  myCameraHandle = new SbigUsbCam("sbigusb0");
  if(myCameraHandle->GetStatus() != CE_NO_ERROR)
  {
    // The was an error
    string str("Cannot open SBIG USB camera: ");
    str = str + string(myCameraHandle->GetStatusString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // No Error so try to establish a link
  myCameraHandle->EstablishLink();
  if(myCameraHandle->GetStatus() != CE_NO_ERROR)
  {
    // The was an error
    string str("Cannot establish link to SBIG USB camera: ");
    str = str + string(myCameraHandle->GetStatusString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Set safe defaults for the CCD imaging region
  myGIP.ccd = 0;  // We only have an imaging CCD
  myGIP.abgState = 0; // ABG shut off
  myGIP.openShutter = 2; // Default to keep shutter closed.  This is for safety reasons (e.g. don't want the thing to just open up
  myGIP.exposureTime = MIN_ST7_EXPOSURE;
  myGIP.readoutMode = 0;
  myGIP.top = 0; // top line
  myGIP.left = 0; // left column
  myGIP.width = 0;  // Lets CalcImageParams calculate this value
  myGIP.height = 0; // Lets CalcImageParams calculate this value
  myGIP.subtract = 0; // Don't perform subtraction

  // Get Driver Info
  myGDIP.request = 0;   // Perform a standard request
  if(myCameraHandle->GetDriverInfo(&myGDIP, (void *)&myGDIR0) != CE_NO_ERROR)
  {
    // The was an error
    string str("Cannot get driver info for SBIG USB camera: ");
    str = str + string(myCameraHandle->GetStatusString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // Set the Miscellaneous control parameters
  myMCP.fanEnable = (unsigned short) true;   // default have the fan running
  myMCP.shutterCommand = 2; // default have the shutter close
  myMCP.ledState = 0;       // keep the LED off
  myCameraHandle->MiscellaneousControl(&myMCP);

  myCameraHandle->QueryTemperatureStatus(&myQTSR);

  // Get CCD Info
  myGCCDIP.request = (myGIP.ccd == CCD_IMAGING) ? 0 : 1;
  myCameraHandle->GetCCDInfo(&myGCCDIP, (void *)&myGCCDIR0);

  if(myCameraHandle->GetStatus() != CE_NO_ERROR)
  {
    // The was an error
    string str("Cannot get CCD info for SBIG USB camera: ");
    str = str + string(myCameraHandle->GetStatusString());
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    myCameraOK = false;
    return false;
  }

  // We sucessfully initialized the camera
  myCameraOK = true;
  
  // Set the CCD temperature
  SetCCDTemp(myCCDSetTemp);

  // For good measure read the CCD temperature
  ReadTemperatureStatus();

  return myCameraOK;
}

/** Aborts the current exposure */
void CameraSBIG::KillExposure(){
  EndExposureParams myEEP;
  int status;

  myEEP.ccd = 0;
  
  if(myCameraOK)
  {
    if((status = myCameraHandle->EndExposure(&myEEP)) != CE_NO_ERROR)
    {
      // The was an error
      string str("Cannot kill exposure: ");
      str = str + string(myCameraHandle->GetStatusString());
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
    }
  }                
}

/** Sets the setpoint for the CCD */
void CameraSBIG::SetCCDTemp(double CCDTemp){
  myCCDSetTemp = CCDTemp;

  // place some error bounds on it
  if(myCCDSetTemp > 10.0)
    myCCDSetTemp = 10.0;

  if(myCameraOK)
  {
    mySTRP.regulation = 1;  // enable temperature regulation
    mySTRP.ccdSetpoint = myCameraHandle->CalcSetpoint(myCCDSetTemp);
    myCameraHandle->SetTemperatureRegulation(&mySTRP);

    if(myCameraHandle->GetStatus() != CE_NO_ERROR)
    {
      // The was an error
      string str("Error setting CCD temperature setpoint for SBIG USB camera: ");
      str = str + string(myCameraHandle->GetStatusString());
      mySystem.myConfiguration.myLogger.WriteErrorString(str);
      return;
    }
  }
}

/** Shuts down the camera */
int CameraSBIG::ShutDownCamera(){
  if(myCameraHandle)
  {
    delete myCameraHandle;
    myCameraOK = false;
  }

  return 0;
}

/** Returns the ambient temperature on the CCD chip */
double CameraSBIG::GetAmbientTemp(){
  // to save on bugging the CCD, only the GetCCDTemp actually poles these quantities
  return myAmbientTemp;
}

/** Gets the fan power on the CCD */
int CameraSBIG::GetFanPower(){
  // to save on bugging the CCD, only the GetCCDTemp actually poles these quantities
  return myFanPower;

}
/** Reads camera information from the configuration file */
void CameraSBIG::Read(ifstream & is){
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
}
/** Writes the camera information to the configuration file */
void CameraSBIG::Write(ofstream & os){
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
}
/** No descriptions */
void CameraSBIG::SetImagingSize(int width, int height){
  myWidth = width;
  myHeight = height;
}
/** No descriptions */
void CameraSBIG::SetBinning(int X, int Y){
  myBinX = X;
  myBinY = Y;
}
/** No descriptions */
void CameraSBIG::SetOffsets(int top, int left){
  myTop = top;
  myLeft = left;
}
/** No descriptions */
void CameraSBIG::SetExposureTime(float exp){
  myExposureTime = exp;
}
/** No descriptions */
void CameraSBIG::SetBinX(int X){
  myBinX = X;
}
/** No descriptions */
void CameraSBIG::SetBinY(int Y){
  myBinY = Y;
}
/** No descriptions */
void CameraSBIG::SetWidth(int width){
  myWidth = width;
}
/** No descriptions */
void CameraSBIG::SetHeight(int height){
  myHeight = height;
}
/** No descriptions */
void CameraSBIG::SetLeft(int left){
  myLeft = left;
}
/** No descriptions */
void CameraSBIG::SetTop(int top){
  myTop = top;
}
/** No descriptions */
int CameraSBIG::GetWidth(){
  return myWidth;
}
/** No descriptions */
int CameraSBIG::GetHeight(){
  return myHeight;
}
/** No descriptions */
int CameraSBIG::GetLeft(){
  return myLeft;
}
/** No descriptions */
int CameraSBIG::GetTop(){
  return myTop;
}
/** No descriptions */
int CameraSBIG::GetBinX(){
  return myBinX;
}
/** No descriptions */
int CameraSBIG::GetBinY(){
  return myBinY;
}
/** No descriptions */
double CameraSBIG::GetCCDSetTemp(){
  return myCCDSetTemp;
}
/** No descriptions */
void CameraSBIG::SetDark(bool dark){
  myDark = dark;
}
/** No descriptions */
float CameraSBIG::GetExposureTime(){
  return myExposureTime;
}
/** No descriptions */
void CameraSBIG::ReadTemperatureStatus(){
  if(myCameraOK)
  {
    myCameraHandle->QueryTemperatureStatus(&myQTSR);
    myFanPower = myCameraHandle->CalcPower(&myQTSR);
    myCCDTemp = myCameraHandle->CalcCcdTemperature(&myQTSR);
    myAmbientTemp = myCameraHandle->CalcAmbTemperature(&myQTSR);
  }
}
/** No descriptions */
void CameraSBIG::WriteXML(ofstream & os){
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

  os << "</camera>" << endl;
}
