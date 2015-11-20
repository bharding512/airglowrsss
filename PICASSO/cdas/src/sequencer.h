/***************************************************************************
                          sequencer.h  -  description
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

#ifndef SEQUENCER_H
#define SEQUENCER_H


/**  The class to handle the sequencing of the camera.
     This actually performs the bulk of the work and is called from
     the control thread.  It will move the filter position to the correct
     location and take the exposure.
  *@author Jonathan Makela
  */

#include <stdlib.h>
#include <iostream>
   
#ifdef MAXFILTPOS
#else
  #define MAXFILTPOS 6
#endif

using namespace std;
  
class Sequencer {
public: 
	Sequencer();
	~Sequencer(); 
  /** Starts the sequence by reseting variables. */
  void StartSequence();
  /** Captures the next image in the sequence. */
  void CaptureNextImage();
  /** Sets the emission height for the indicated filter. */
  void SetEmissionHeight(int num, double height);
  /** Sets the filter name for the indicated filter. */
  void SetFilterName(int num, string name);
  /** Sets the exposure time for the indicated filter. */
  void SetExposureTime(int num, double time);
  /** Sets the shutter position (calls thru to FilterWheelCFWL, but does error checking first). */
  void SetShutterPosition(int num);
  /** Write the Sequencer class to the config file. */
  void Write(ofstream &os);
  /** Reads the Sequencer class from the config file. */
  void Read(ifstream &is);
  /** Writes the Sequencer class to XML. */
  void WriteXML(ofstream &os);
  /** Returns the emission height of the indicated filter. */
  double GetEmissionHeight(int num);
  /** Returns the name of the indicated filter. */
  string GetFilterName(int num);
  /** Returns the shutter position. */
  int GetShutterPosition();
  /** Returns the exposure time of the indicated filter. */
  double GetExposureTime(int num);
  /** Creates the master sequence based on the sub-sequences and sequence order. */
  void CreateSequence();
  /** Sets the number of filled positions in the given subsequence. */
  void SetNumberInSubSequence(int subsequence, int num);
  /** Sets the filter for the given position in a subsequence */
  void SetSubSequence(int subsequence, int position, int filter);
  /** Returns the filter in the subsequence in the indicated position. */
  int GetSubSequence(int subsequence, int position);
  /** Returns the number of entries in the requested subsequence. */
  int GetNumberInSubSequence(int subsequence);
  /** Sets the sub-sequence order in position. */
  void SetSubSequenceOrder(int position, int sequence);
  /** Returns the number of sub-sequences in the overall sequence. */
  int GetNumberInSubSequenceOrder();
  /** Sets the number of sub-sequences in the overall sequence. */
  void SetNumberInSubSequenceOrder(int num);
  /** Returns the sequence in position of the sub-sequence order. */
  int GetSubSequenceOrder(int position);
  /** Sets the exposure time for the dark images. */
  void SetDarkExposureTime(double expTime);
  /** Sets the number of images between the dark images. */
  void SetDarkFrequency(int freq);
  /** Sets the doDark flag. */
  void SetDoDark(bool flag);
  /** Returns the exposure time for the dark images. */
  double GetDarkExposureTime();
  /** Returns the number of images between dark images. */
  int GetDarkFrequency();
  /** Returns the doDark flag. */
  bool GetDoDark();
  /** Sets the wait time between frames. */
  void SetPeriodicityFrequency(double freq);
  /** Sets the myDoPeriodicity flag. */
  void SetDoPeriodicity(bool flag);
  /** Returns the time between frames. */
  double GetPeriodicityFrequency();
  /** Returns the Periodicity flag. */
  bool GetDoPeriodicity();
  /** Returns the current filter position. */
  int GetCurrentFilter();
private: // Public attributes
  /** Contains the current number of images during the current data run. */
  int myNumberImages[MAXFILTPOS];
  /** The label for the filters. */
  string myFilterName[MAXFILTPOS];
  /** The exposure time in seconds for each filter. */
  double myExposureTime[MAXFILTPOS];
  /** The sequencing of filters to be used. */
  int mySequence[50];
  /** The height in kilometers of the emission. */
  double myEmissionHeight[MAXFILTPOS];
  /** The total number in the sequence. */
  int myNumberInSequence;
  /** The current location in the sequence. */
  int myCurrentNumberInSequence;
  /** The subsequences that can make up the final sequence. */
  int mySubSequences[4][10];
  /** The total number in each subsequence */
  int myNumberInSubSequence[4];
  /** Flag on whether to do dark images or not. */
  bool myDoDark;
  /** The frequency (in nunber of images) to take dark images within the sequence. */
  int myDarkFrequency;
  /** The exposure time in seconds for the dark images. */
  double myDarkExposureTime;
  /** Flag to determine if periodicity will be used. */
  bool myDoPeriodicity;
  /** The number of minutes to wait in between images. */
  double myPeriodicityFrequency;
  /** The order in which to perform the subsequences. */
  int mySubSequenceOrder[5];
  /** The current number of images since the last dark image was taken. */
  int myCurrentDarkCount;
  /** The number of sub sequences in the overall sequence */
  int myNumberInSubSequenceOrder;
  /** The total number of images in the current sequence */
  int myTotalImages;
  /** The current filter position. */
  int myCurrentFilter;
};

#endif
