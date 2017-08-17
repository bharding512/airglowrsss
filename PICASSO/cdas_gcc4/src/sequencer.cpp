/***************************************************************************
                          sequencer.cpp  -  description
                             -------------------
    begin                : Mon May 24 2004
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

#include "sequencer.h"
#include <sstream>
#include "system.h"

extern System mySystem;

Sequencer::Sequencer(){
  // Initialize variables
  int i;
  for(i = 0; i < MAXFILTPOS; i++)
  {
    myNumberImages[i] = 0;
    myExposureTime[i] = 0.;
    myFilterName[i] = "";
    myEmissionHeight[i] = 0.;
  }

  for(i = 0; i < 10; i++)
  {
    // Initialize the four subsequences
    mySubSequences[0][i] = -1;
    mySubSequences[1][i] = -1;
    mySubSequences[2][i] = -1;
    mySubSequences[3][i] = -1;
  }

  myNumberInSubSequence[0] = 0;
  myNumberInSubSequence[1] = 0;
  myNumberInSubSequence[2] = 0;
  myNumberInSubSequence[3] = 0;
  
  for(i = 0; i < 50; i++)
  {
    // Initialize the sequence to -1 (invalid)
    mySequence[i] = -1;
  }

  myNumberInSequence = 0;
  myCurrentNumberInSequence = 0;

  for(i = 0; i < 5; i++)
  {
    mySubSequenceOrder[i] = -1;
  }

  myNumberInSubSequenceOrder = 0;

  myDoDark = false;
  myDarkFrequency = 0;
  myDarkExposureTime = 0.;
  myCurrentDarkCount = 0;
  myTotalImages = 0;

  myDoPeriodicity = false;
  myPeriodicityFrequency = 0.;
}

Sequencer::~Sequencer(){
}

/** Captures the next image in the sequence. */
void Sequencer::CaptureNextImage(){
//  int myCurrentFilter;
  string filename;
  char num[5];
  
  time_t start;
  start = time(NULL);
  char buff[20];
  strftime(buff, sizeof(buff), "%Y%m%d_%H%M%S", gmtime(&start));
  
  // Check to make sure that the current cycle number is out of bounds
  if(myCurrentNumberInSequence < 0 || myCurrentNumberInSequence >= myNumberInSequence)
    myCurrentNumberInSequence = 0;

  // See if we are doing a dark image or not
  if(myDoDark && (myTotalImages % myDarkFrequency == 0))
  {
    // Do a dark image
    // Move the filterwheel to the closed position
    mySystem.myConfiguration.myFilterWheel->CloseShutter();

    while(mySystem.myConfiguration.myFilterWheel->GetCurrentPosition() == -999)
      // Filterwheel still moving
      sleep(1);

    // Set the dark flag
    mySystem.myConfiguration.myCamera->SetDark(true);

    // Set the exposure time
    mySystem.myConfiguration.myCamera->SetExpTime(myDarkExposureTime);

    // Create the filename
    sprintf(num, "%.4d", myCurrentDarkCount);
//    filename = mySystem.myConfiguration.myDirectories.GetTodaysDir() + "DARK_" + string(num);
    filename = mySystem.myConfiguration.myDirectories.GetTodaysDir() + mySystem.myConfiguration.mySite.GetAbbreviation() + "_DARK_" + buff + "_" + num;

    // Capture the image
    mySystem.myConfiguration.myCamera->CaptureImage(mySystem.myConfiguration.myImage, mySystem.myConfiguration.myDirectories.GetDataPath() + filename, false);

    // Save a PNG
    if(mySystem.myConfiguration.myDirectories.GetDoPNG())
      mySystem.myConfiguration.myImage->SavePNG(mySystem.myConfiguration.myDirectories.GetPNGPath() + filename);

    // Increment the DarkCount
    myCurrentDarkCount++;
  }

  // Get the current time (for comparison later on if the Periodicity is set)
//  time_t start;
//  start = time(NULL);

  // Now do the next image in the cycle
  start = time(NULL);
  strftime(buff, sizeof(buff), "%Y%m%d_%H%M%S", gmtime(&start));

  // Make sure the shutter is open
  mySystem.myConfiguration.myFilterWheel->OpenShutter();

  // Move the filterwheel to the current position
  myCurrentFilter = mySequence[myCurrentNumberInSequence];
  mySystem.myConfiguration.myFilterWheel->MoveToPosition(myCurrentFilter);

  while(mySystem.myConfiguration.myFilterWheel->GetCurrentPosition() == -999)
    // Filterwheel still moving
    sleep(1);

  // Set the dark flag
  mySystem.myConfiguration.myCamera->SetDark(false);

  // Set the exposure time
  mySystem.myConfiguration.myCamera->SetExpTime(myExposureTime[myCurrentFilter-1]);

  // Create the filename
  sprintf(num, "%.4d", myNumberImages[myCurrentFilter-1]);
  
//  filename = mySystem.myConfiguration.myDirectories.GetTodaysDir() + myFilterName[myCurrentFilter-1] + "_" + num;
  filename = mySystem.myConfiguration.myDirectories.GetTodaysDir() + mySystem.myConfiguration.mySite.GetAbbreviation() + "_" + myFilterName[myCurrentFilter-1] + "_" + buff + "_" + num;

  // Capture the image
  mySystem.myConfiguration.myCamera->CaptureImage(mySystem.myConfiguration.myImage, mySystem.myConfiguration.myDirectories.GetDataPath() + filename, false);

  // Save a PNG
  if(mySystem.myConfiguration.myDirectories.GetDoPNG())
    mySystem.myConfiguration.myImage->SavePNG(mySystem.myConfiguration.myDirectories.GetPNGPath() + filename);

  // Increment the Image counts
  myNumberImages[myCurrentFilter-1]++;
  myCurrentNumberInSequence++;
  myTotalImages++;

  // Check periodicity
  if(myDoPeriodicity)
  {
    // Get the current time
    time_t now;
    now = time(NULL);

    // The MyPeriodicityFrequency is in seconds
    while(difftime(now, start) < myPeriodicityFrequency)
    {
      // We still need to wait some
      sleep(1);
      now = time(NULL);
    }
  }
}

/** Starts the sequence by reseting variables. */
void Sequencer::StartSequence(){
  // Reset various parameters
  myCurrentNumberInSequence = 0;
  myCurrentDarkCount = 0;
  myTotalImages = 0;

  // Reset the count for each of the filters
  for(int i = 0; i < mySystem.myConfiguration.myFilterWheel->GetNumFilters(); i++)
    myNumberImages[i] = 0;

  mySystem.myConfiguration.myDirectories.SetTodaysDir();
}
/** Sets the filter name for the indicated filter. */
void Sequencer::SetFilterName(int num, string name){
  // Make sure the requested filter number is valid
  if(num < 1 || num > mySystem.myConfiguration.myFilterWheel->GetNumFilters())
    return;
      
  myFilterName[num-1] = name;
}
/** Sets the emission height for the indicated filter. */
void Sequencer::SetEmissionHeight(int num, double height){
  // Make sure the requested filter number is valid
  if(num < 1 || num > mySystem.myConfiguration.myFilterWheel->GetNumFilters())
    return;

  myEmissionHeight[num-1] = height;  
}
/** Sets the exposure time for the indicated filter. */
void Sequencer::SetExposureTime(int num, double time){
  // Make sure the requested filter number is valid
  if(num < 1 || num > mySystem.myConfiguration.myFilterWheel->GetNumFilters())
    return;

  myExposureTime[num-1] = time;
}
/** Sets the shutter position (calls thru to FilterWheelCFWL, but does error checking first). */
void Sequencer::SetShutterPosition(int num){
  // Make sure the requested filter number is valid
  if(num < 1 || num > mySystem.myConfiguration.myFilterWheel->GetNumFilters())
    return;

  mySystem.myConfiguration.myFilterWheel->SetShutterPosition(num);
}

/** Write the Sequencer class to the config file. */
void Sequencer::Write(ofstream &os){
  // Put the header label
  os << "#Sequencer" << endl;

  // Write the settings to the file
  for(int i = 0; i < mySystem.myConfiguration.myFilterWheel->GetNumFilters(); i++)
  {
    os << myFilterName[i] << endl;
    os << myEmissionHeight[i] << endl;
    os << myExposureTime[i] << endl;
  }

  os << myDoDark << endl;
  os << myDarkFrequency << endl;
  os << myDarkExposureTime << endl;

  os << myDoPeriodicity << endl;
  os << myPeriodicityFrequency << endl;

  // Write out the sequence information
  os << myNumberInSubSequenceOrder << endl;
  for(int i = 0; i < myNumberInSubSequenceOrder; i++)
  {
    os << mySubSequenceOrder[i] << endl;
  }

  for(int j = 0; j < 4; j++)
  {
    os << myNumberInSubSequence[j] << endl;
    for(int k = 0; k < myNumberInSubSequence[j]; k++)
      os << mySubSequences[j][k] << endl;
  }
}

/** Reads the Sequencer class from the config file. */
void Sequencer::Read(ifstream &is){
  // Search for the header from the beginning of the file
  is.seekg(0,ios::beg);

  string myLine;

  // Read in the next line of the file
  getline(is, myLine);

  // Find the header
  while(myLine != "#Sequencer" && !is.eof())
    getline(is, myLine);

  if(is.eof())
  {
    mySystem.myConfiguration.myLogger.WriteErrorString("End of configuration file reached before #Sequencer found");
    return;
  }
  
  // Read the data
  for(int i = 0; i < mySystem.myConfiguration.myFilterWheel->GetNumFilters(); i++)
  {
    is >> myFilterName[i];
    is >> myEmissionHeight[i];
    is >> myExposureTime[i];
  }

  is >> myDoDark;
  is >> myDarkFrequency;
  is >> myDarkExposureTime;

  is >> myDoPeriodicity;
  is >> myPeriodicityFrequency;

  // Write out the sequence information
  is >> myNumberInSubSequenceOrder;
  for(int i = 0; i < myNumberInSubSequenceOrder; i++)
  {
    is >> mySubSequenceOrder[i];
  }

  for(int j = 0; j < 4; j++)
  {
    is >> myNumberInSubSequence[j];
    for(int k = 0; k < myNumberInSubSequence[j]; k++)
      is >> mySubSequences[j][k];
  }

  // Create the master sequence
  CreateSequence();
}

/** Writes the Sequencer class to XML. */
void Sequencer::WriteXML(ofstream &os){
}
/** Returns the name of the indicated filter. */
string Sequencer::GetFilterName(int num){
  if(num < 1 || num > mySystem.myConfiguration.myFilterWheel->GetNumFilters())
    return "ERROR";
    
  return myFilterName[num-1];
}
/** Returns the emission height of the indicated filter. */
double Sequencer::GetEmissionHeight(int num){
  if(num < 1 || num > mySystem.myConfiguration.myFilterWheel->GetNumFilters())
    return -999.;

  return myEmissionHeight[num-1];
}
/** Returns the exposure time of the indicated filter. */
double Sequencer::GetExposureTime(int num){
  if(num < 1 || num > mySystem.myConfiguration.myFilterWheel->GetNumFilters())
    return -999.;

  return myExposureTime[num-1];
}

/** Returns the shutter position. */
int Sequencer::GetShutterPosition(){
  return mySystem.myConfiguration.myFilterWheel->GetShutterPosition();
}

/** Creates the master sequence based on the sub-sequences and sequence order. */
void Sequencer::CreateSequence(){
  int myIndex;
  ostringstream sout;
  sout << "Filter Sequence";

  // Initialize the sequence counter
  myNumberInSequence = 0;

  for(int i = 0; i < myNumberInSubSequenceOrder; i++)
  {
    myIndex = mySubSequenceOrder[i];

    // We need to add another subsequence
    for(int j = 0; j < myNumberInSubSequence[myIndex]; j++)
    {
      mySequence[myNumberInSequence] = mySubSequences[myIndex][j];
      sout << ", " << mySubSequences[myIndex][j];
      myNumberInSequence++;
    }
  }

  // Write the sequence to the log
  mySystem.myConfiguration.myLogger.WriteSystemLogString(sout.str());
}
/** Sets the filter for the given position in a subsequence */
void Sequencer::SetSubSequence(int subsequence, int position, int filter){
  // Subsequences are 0 indexed
  subsequence--;
  
  if(subsequence < 0 || subsequence > 3)
    if(position < 0 || position >= myNumberInSubSequence[subsequence])
    return;
    
  mySubSequences[subsequence][position] = filter;
}

/** Sets the number of filled positions in the given subsequence. */
void Sequencer::SetNumberInSubSequence(int subsequence, int num){
  // Subsequences are 0 indexed
  subsequence--;
  
  if(subsequence < 0 || subsequence > 3 || num < 0 || num > 9)
    return;

  myNumberInSubSequence[subsequence] = num;
}

/** Returns the number of entries in the requested subsequence. */
int Sequencer::GetNumberInSubSequence(int subsequence){
  // Subsequences are 0 index
  subsequence--;
  
  if(subsequence < 0 || subsequence > 3)
    return 0;

  return myNumberInSubSequence[subsequence];
}

/** Returns the filter in the subsequence in the indicated position. */
int Sequencer::GetSubSequence(int subsequence, int position){
  // Subsequences are 0 index
  subsequence--;
  
  if(subsequence < 0 || subsequence > 3)
    if(position < 0 || position >= myNumberInSubSequence[subsequence])
      return 0;

  return mySubSequences[subsequence][position];
}
/** Sets the number of sub-sequences in the overall sequence. */
void Sequencer::SetNumberInSubSequenceOrder(int num){
  if(num < 0 || num > 5)
    return;

  myNumberInSubSequenceOrder = num;
}
/** Returns the number of sub-sequences in the overall sequence. */
int Sequencer::GetNumberInSubSequenceOrder(){
  return myNumberInSubSequenceOrder;
}
/** Sets the sub-sequence order in position. */
void Sequencer::SetSubSequenceOrder(int position, int sequence){
  // Subsequences are 0 index
  sequence--;

  if(sequence < 0 || sequence > 3 || position < 0 || position >= myNumberInSubSequenceOrder)
    return;

  mySubSequenceOrder[position] = sequence;
}
/** Returns the sequence in position of the sub-sequence order. */
int Sequencer::GetSubSequenceOrder(int position){
  if(position < 0 || position >= myNumberInSubSequenceOrder)
    return 0;

  // Subsequences are 0 index
  return mySubSequenceOrder[position] + 1;
}
/** Sets the doDark flag. */
void Sequencer::SetDoDark(bool flag){
  myDoDark = flag;
}
/** Sets the number of images between the dark images. */
void Sequencer::SetDarkFrequency(int freq){
  myDarkFrequency = freq;
}
/** Sets the exposure time for the dark images. */
void Sequencer::SetDarkExposureTime(double expTime){
  myDarkExposureTime = expTime;
}
/** Returns the doDark flag. */
bool Sequencer::GetDoDark(){
  return myDoDark;
}
/** Returns the number of images between dark images. */
int Sequencer::GetDarkFrequency(){
  return myDarkFrequency;
}
/** Returns the exposure time for the dark images. */
double Sequencer::GetDarkExposureTime(){
  return myDarkExposureTime;
}
/** Sets the myDoPeriodicity flag. */
void Sequencer::SetDoPeriodicity(bool flag){
  myDoPeriodicity = flag;
}
/** Sets the wait time between frames. */
void Sequencer::SetPeriodicityFrequency(double freq){
  myPeriodicityFrequency = freq;
}
/** Returns the Periodicity flag. */
bool Sequencer::GetDoPeriodicity(){
  return myDoPeriodicity;
}
/** Returns the time between frames. */
double Sequencer::GetPeriodicityFrequency(){
  return myPeriodicityFrequency;
}
/** Returns the current filter position. */
int Sequencer::GetCurrentFilter(){
  return myCurrentFilter;
}
