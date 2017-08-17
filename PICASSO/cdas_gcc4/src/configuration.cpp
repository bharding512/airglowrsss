/***************************************************************************
                          configuration.cpp  -  description
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

#include "configuration.h"
#include "system.h"

extern System mySystem;

Configuration::Configuration(){
  myCamera = new Camera();
  myImage = new Image();
  myFilterWheel = new FilterWheel();
}

Configuration::~Configuration(){
  myFilterWheel->~FilterWheel();
  free(myFilterWheel);
  myCamera->~Camera();
  free(myCamera);
  myImage->~Image();
  free(myImage);
}

/** Copy constructor */
Configuration::Configuration(Configuration& c){
  Camera *myCamera = c.myCamera;
  Directories myDirectories = c.myDirectories;
  FilterWheel *myFilterWheel = c.myFilterWheel;
  Logger myLogger = c.myLogger;
  Schedule mySchedule = c.mySchedule;
  Site mySite = c.mySite;
  User myUser = c.myUser;
  Sequencer mySequencer = c.mySequencer;
}
/** Loads a configuration in from filename */
bool Configuration::LoadConfiguration(string filename){
  ifstream myFile;

  // open the file
  myFile.open(filename.c_str(), ifstream::in);

  if(!myFile)
  {
    // Error occured opening the file
    string str("Cannot open configuration file ");
    str = str + filename;
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    return false;
  }

  // Read in the structures
  myLogger.Read(myFile);
  myDirectories.Read(myFile);
  mySite.Read(myFile);
  myUser.Read(myFile);
  myCamera->Read(myFile);
  myFilterWheel->Read(myFile);
  mySchedule.Read(myFile);
  mySequencer.Read(myFile);
    
//  if(myCameraType == SBIGUSB)
  if(myCameraType == SBIGUNIV)
  {
    // Save the needed variables to temporary holders
    int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
    double tempTemp;
    float tempGain;
    bool tempRegulation;

    tempXBin = myCamera->GetBinX();
    tempYBin = myCamera->GetBinY();
    tempTop = myCamera->GetTop();
    tempLeft = myCamera->GetLeft();
    tempHeight = myCamera->GetHeight();
    tempWidth = myCamera->GetWidth();
    tempTemp = myCamera->GetCCDSetTemp();
		tempGain = myCamera->GetGain();
		tempRegulation = myCamera->GetAutoTemp();
    
    // close the old camera
    myCamera->~Camera();
    delete(myCamera);

    // open the new one
    //myCamera = new CameraSBIG();
    myCamera = new CameraSBIGUniv();
    myCamera->SetBinX(tempXBin);
    myCamera->SetBinY(tempYBin);
    myCamera->SetTop(tempTop);
    myCamera->SetLeft(tempLeft);
    myCamera->SetHeight(tempHeight);
    myCamera->SetWidth(tempWidth);
    myCamera->SetCCDTemp(tempTemp);
    myCamera->SetGain(tempGain);
    myCamera->SetAutoTemp(tempRegulation);
  }

  if(myCameraType == ANDOR)
  {
    // Save the needed variables to temporary holders
	  int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
	  double tempTemp;
	  float tempGain;
	  bool tempRegulation;

	  tempXBin = myCamera->GetBinX();
	  tempYBin = myCamera->GetBinY();
	  tempTop = myCamera->GetTop();
	  tempLeft = myCamera->GetLeft();
	  tempHeight = myCamera->GetHeight();
	  tempWidth = myCamera->GetWidth();
	  tempTemp = myCamera->GetCCDSetTemp();
	  tempGain = myCamera->GetGain();
	  tempRegulation = myCamera->GetAutoTemp();
    
    // close the old camera
	  myCamera->~Camera();
	  delete(myCamera);

    // open the new one
    //myCamera = new CameraSBIG();
	  myCamera = new CameraAndor();
	  myCamera->SetBinX(tempXBin);
	  myCamera->SetBinY(tempYBin);
	  myCamera->SetTop(tempTop);
	  myCamera->SetLeft(tempLeft);
	  myCamera->SetHeight(tempHeight);
	  myCamera->SetWidth(tempWidth);
	  myCamera->SetCCDTemp(tempTemp);
	  myCamera->SetGain(tempGain);
	  myCamera->SetAutoTemp(tempRegulation);
  }
  
  if(myCameraType == FLIPROLINE)
  {
    // Save the needed variables to temporary holders
	  int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
	  double tempTemp;
	  float tempGain;
	  bool tempRegulation;

	  tempXBin = myCamera->GetBinX();
	  tempYBin = myCamera->GetBinY();
	  tempTop = myCamera->GetTop();
	  tempLeft = myCamera->GetLeft();
	  tempHeight = myCamera->GetHeight();
	  tempWidth = myCamera->GetWidth();
	  tempTemp = myCamera->GetCCDSetTemp();
	  tempGain = myCamera->GetGain();
	  tempRegulation = myCamera->GetAutoTemp();
    
    // close the old camera
	  myCamera->~Camera();
	  delete(myCamera);

    // open the new one
    //myCamera = new CameraSBIG();
	  myCamera = new CameraFLI();
	  myCamera->SetBinX(tempXBin);
	  myCamera->SetBinY(tempYBin);
	  myCamera->SetTop(tempTop);
	  myCamera->SetLeft(tempLeft);
	  myCamera->SetHeight(tempHeight);
	  myCamera->SetWidth(tempWidth);
	  myCamera->SetCCDTemp(tempTemp);
	  myCamera->SetGain(tempGain);
	  myCamera->SetAutoTemp(tempRegulation);
  }
  
/*if(myCameraType == PVCAM)
  {
    // Save the needed variables to temporary holders
	  int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
	  double tempTemp;
	  float tempGain;
	  bool tempRegulation;

	  tempXBin = myCamera->GetBinX();
	  tempYBin = myCamera->GetBinY();
	  tempTop = myCamera->GetTop();
	  tempLeft = myCamera->GetLeft();
	  tempHeight = myCamera->GetHeight();
	  tempWidth = myCamera->GetWidth();
	  tempTemp = myCamera->GetCCDSetTemp();
	  tempGain = myCamera->GetGain();
	  tempRegulation = myCamera->GetAutoTemp();
    
    // close the old camera
	  myCamera->~Camera();
	  delete(myCamera);

    // open the new one
    //myCamera = new CameraSBIG();
	  myCamera = new CameraPVCAM();
	  myCamera->SetBinX(tempXBin);
	  myCamera->SetBinY(tempYBin);
	  myCamera->SetTop(tempTop);
	  myCamera->SetLeft(tempLeft);
	  myCamera->SetHeight(tempHeight);
	  myCamera->SetWidth(tempWidth);
	  myCamera->SetCCDTemp(tempTemp);
	  myCamera->SetGain(tempGain);
	  myCamera->SetAutoTemp(tempRegulation);
  }*/
/*  if(myCameraType == APOGEE)
  {
    // Save the needed variables to temporary holders
          int tempXBin, tempYBin, tempTop, tempLeft, tempHeight, tempWidth;
          double tempTemp;
          float tempGain;
          bool tempRegulation;

          tempXBin = myCamera->GetBinX();
          tempYBin = myCamera->GetBinY();
          tempTop = myCamera->GetTop();
          tempLeft = myCamera->GetLeft();
          tempHeight = myCamera->GetHeight();
          tempWidth = myCamera->GetWidth();
          tempTemp = myCamera->GetCCDSetTemp();
          tempGain = myCamera->GetGain();
          tempRegulation = myCamera->GetAutoTemp();
    
    // close the old camera
          myCamera->~Camera();
          delete(myCamera);

    // open the new one
    //myCamera = new CameraSBIG();
          myCamera = new CameraAPOGEE();
          myCamera->SetBinX(tempXBin);
          myCamera->SetBinY(tempYBin);
          myCamera->SetTop(tempTop);
          myCamera->SetLeft(tempLeft);
          myCamera->SetHeight(tempHeight);
          myCamera->SetWidth(tempWidth);
          myCamera->SetCCDTemp(tempTemp);
          myCamera->SetGain(tempGain);
          myCamera->SetAutoTemp(tempRegulation);
  }*/

  
  // Load the filterwheel)
  if(myFilterWheelType == SBIGCFWL)
  {
    // Save the needed variables to temporary holders
    int tempNumFilters, tempShutterPos;
    tempNumFilters = myFilterWheel->myNumFilters;
    tempShutterPos = myFilterWheel->myShutterPos;

    // close the old filterwheel
    myFilterWheel->~FilterWheel();
    delete(myFilterWheel);

    // open the new one
    myFilterWheel = new FilterWheelCFWL();
    myFilterWheel->myNumFilters = tempNumFilters;
    myFilterWheel->myShutterPos = tempShutterPos;
  }

  if(myFilterWheelType == FLIFWL)
  {
    // Save the needed variables to temporary holders

    int tempNumFilters, tempShutterPos;
    tempNumFilters = myFilterWheel->myNumFilters;
    tempShutterPos = myFilterWheel->myShutterPos;

    // close the old filterwheel
    myFilterWheel->~FilterWheel();
    delete(myFilterWheel);

    // open the new one
    myFilterWheel = new FilterWheelFLI();
    myFilterWheel->myNumFilters = tempNumFilters;
    myFilterWheel->myShutterPos = tempShutterPos;
  }

  if(myFilterWheelType == KEOFW)
  {
    // Save the needed variables to temporary holders
    int tempNumFilters, tempShutterPos;
    char tempFWDevice[256];
    char tempTempControlDevice[256];
    tempNumFilters = mySystem.myConfiguration.myFilterWheel->GetNumFilters();
    tempShutterPos = mySystem.myConfiguration.myFilterWheel->GetShutterPosition();
//    strcpy(tempFWDevice,myFilterWheel->myFWDevice);
//    strcpy(tempTempControlDevice,myFilterWheel->myTempControlDevice);
//    tempFWDevice = myFilterWheel->myFWDevice.c_str();
//    tempTempControlDevice = myFilterWheel->myTempControlDevice.c_str();

    // close the old filterwheel
    myFilterWheel->~FilterWheel();
    delete(myFilterWheel);

    // open the new one
    myFilterWheel = new FilterWheelKEO();
    mySystem.myConfiguration.myFilterWheel->SetNumFilters(tempNumFilters);
    mySystem.myConfiguration.myFilterWheel->SetShutterPosition(tempShutterPos);
//    myFilterWheel->SetFWDevice(tempFWDevice);
//    myFilterWheel->SetTempControlDevice(tempTempControlDevice);
    //strcpy(myFilterWheel->myFWDevice,tempFWDevice);
    //strcpy(myFilterWheel->myTempControlDevice,tempTempControlDevice);
//    myFilterWheel->myFWDevice = tempFWDevice.c_str();
//    myFilterWheel->myTempControlDevice = tempTempControlDevice.c_str();
//    myFilterWheel->SetFWDevice(tempFWDevice.c_str());
//    myFilterWheel->SetTempControlDevice(tempTempControlDevice.c_str());
    myFilterWheel->SetUp_FilterWheel();
    myFilterWheel->OpenShutter();
    sleep(1);
    myFilterWheel->CloseShutter();
  }

  if(myFilterWheelType == KEOSER)
  {
    // Save the needed variables to temporary holders
    int tempNumFilters, tempShutterPos;
    tempNumFilters = myFilterWheel->myNumFilters;
    tempShutterPos = myFilterWheel->myShutterPos;

    // close the old filterwheel
    myFilterWheel->~FilterWheel();
    delete(myFilterWheel);

    // open the new one
    myFilterWheel = new FilterWheelSERIAL();
    myFilterWheel->myNumFilters = tempNumFilters;
    myFilterWheel->myShutterPos = tempShutterPos;
    myFilterWheel->SetUp_FilterWheel();

    // turn to a shutter position
    myFilterWheel->OpenShutter();
    sleep(1);
    myFilterWheel->CloseShutter();
  }

  // Close the file and return
  myFile.close();

  // Get the rise and set times
  mySystem.myConfiguration.mySchedule.SetTodaysRiseSetTimes();
  return true;
}

/** Saves the configuration information */
bool Configuration::SaveConfiguration(string filename){
  ofstream myFile;

  // open the file
  myFile.open(filename.c_str(), ofstream::out);

  if(!myFile)
  {
    // Error occured opening the file
    string str("Cannot open configuration file ");
    str = str + filename;
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    return false;
  }

  // Write out the structures
  myLogger.Write(myFile);
  myDirectories.Write(myFile);
  mySite.Write(myFile);
  myUser.Write(myFile);
  myCamera->Write(myFile);
  myFilterWheel->Write(myFile);
  mySchedule.Write(myFile);
  mySequencer.Write(myFile);

  // Close the file and return
  myFile.close();
  return true;
}
/** Parses a line of XML.  Basically returns the string between <XXX>YYY</XXX>, returning YYY.  Does absolutely no error checking on if it is a valid XML string. */
string Configuration::ParseXML(string input){
  string temp;
  // this is cleaner since strtok needs to diddle with the pointer - ESM
  // keeps gcc happier, too
  char inputLine[255];
  strncpy(inputLine,input.c_str(),255);

  strtok(inputLine, ">");
  temp = strtok(NULL, "<");

  return temp;
}
/** Saves the configuration out to an XML file. */
bool Configuration::SaveXML(string filename){
  ofstream myFile;

  // open the file
  myFile.open(filename.c_str(), ofstream::out);

  if(!myFile)
  {
    // Error occured opening the file
    string str("Cannot open XML file ");
    str = str + filename;
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    return false;
  }

  // Write out the structures
  myFile << "<?xml version=\"1.0\"?>" << endl;
  myFile << "<configuration>" << endl;
  myLogger.WriteXML(myFile);
  myDirectories.WriteXML(myFile);
  mySite.WriteXML(myFile);
  myUser.WriteXML(myFile);
  myCamera->WriteXML(myFile);
  myFilterWheel->WriteXML(myFile);
  mySchedule.WriteXML(myFile);
  mySequencer.WriteXML(myFile);

  myFile << "<systemrunning>" << mySystem.myIsRunning << "</systemrunning>" << endl;

  myFile << "</configuration>" << endl;

  // Close the file and return
  myFile.close();
  return true;
}
