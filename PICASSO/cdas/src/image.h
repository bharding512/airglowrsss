/***************************************************************************
                          image.h  -  description
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

#ifndef IMAGE_H
#define IMAGE_H

#include "site.h"
#include "calibration.h"
#include <png.h>

/**
  *@author Jonathan Makela
  */

#define IMAGETYPE unsigned short
#define BYTE char
#define DWORD unsigned long
#define WORD unsigned int
  
class Image {
public: 
	Image();
	~Image();
  /** Writes the data to an extended tif file. */
//  void SaveTIF(ofstream &of);
  void SaveTIF(string filename);
  /** Reads a TIF file into the Image structure */
  void LoadTIF(ifstream &is);
  /** No descriptions */
  IMAGETYPE* GetDataPointer();
  /** No descriptions */
  void SetCols(int col);
  /** No descriptions */
  void SetRows(int row);
  /** No descriptions */
  void SetRowsCols(int rows, int cols);
  /** No descriptions */
  void SetCCDTemp(float temp);
  /** No descriptions */
  void SetExposureTime(float expTime);
  /** No descriptions */
  void SetBinning(int x, int y);
  /** No descriptions */
  void SetSite(Site *site);
  /** No descriptions */
  void SetTimesNow();
  /** No descriptions */
  void CalcStdDev();
  /** No descriptions */
  void CalcMean();
  /** No descriptions */
  void SetStdDev(double dev);
  /** No descriptions */
  void SetMean(double mean);
  /** Sets the myMin and myMax values for the purpose of autocontrasting an image.  The values of myMin and myMax are set to be +/- one stddev from the mean.

Autocontrast assumes that myMean and myStdDev have been set by calls to CalcMean() and CalcStdDev() respectively. */
  void AutoContrast();
  /** Saves a PNG version of the image.  The image is scaled to an 8-bit image based on myMin and myMax. */
  void SavePNG(string filename);
  /** Normalizes the data between myMin and myMax to 0-255 */
  void Normalize();
  /** Sets the filter name. */
  void SetFilterName(string name);
  /** Normalizes the data between myMin and myMax to 0-255 and subtracts the dark image (if dark == true) */
  void Normalize(bool dark);
  /** Sets the CCD gain. */
  void SetGain(int gain);
private: // Private attributes
  /**  */
  IMAGETYPE *myData;
  int myRows;
  int myCols;
  IMAGETYPE myMin;
  IMAGETYPE myMax;
  struct tm myLocalTime;
  struct tm myUTTime;
  string myFilterName;
  float myExposureTime;
  int myGain;  // Contains the gain times 10
  int myXBinning;
  int myYBinning;
  float myCCDTemp;
  float myFWTemp;
  float myXSpacing;
  float myYSpacing;
  Site mySite;
  string myComment;
  float myProjectionAltitude;
  /**  */
  int myDataOffset;
  /**  */
  int myCustomHeaderOffset;
  /**  */
  Calibration myCalibrationCoefs;
  /**  */
  float myEmissionHeight;
  /**  */
  double myMean;
  double myStdDev;
  float myProjectionLon;
private: // Private methods
  /** No descriptions */
  void SaveCustomHeader(ofstream &of);
public: // Public attributes
public: // Public attributes
  /**  */
  png_byte *myNormalized;
};

#endif
