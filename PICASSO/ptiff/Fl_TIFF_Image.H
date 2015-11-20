//
// "$Id: Fl_TIFF_Image.H,v 1.1 2002/10/11 03:13:52 easysw Exp $"
//
// PhotoCD image header file for flphoto.
//
// Copyright 2002 by Michael Sweet.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//

#ifndef Fl_TIFF_Image_H
#define Fl_TIFF_Image_H
#  include <stdio.h>
#  include <FL/Fl_Image.H>
#include "time.h"
#include "math.h"

#define IMAGETYPE unsigned short

typedef unsigned short ushort;

class FL_EXPORT Fl_TIFF_Image : public Fl_RGB_Image {
  void load_image(FILE *fp);
  void dark_image(FILE *fp);

  public:
  Fl_TIFF_Image(const char* filename);
  Fl_TIFF_Image(FILE *fp, FILE *fpd);
 
  static Fl_Image *check(const char *name, uchar *header, int headerlen);

private:
  IMAGETYPE *myData;
  unsigned short myRows;
  unsigned short myCols;
  IMAGETYPE myMin;
  IMAGETYPE myMax;
  struct tm myLocalTime;
  struct tm myUTTime;
  char myFilterName[12];
  float myExposureTime;
  int myGain;  // Contains the gain times 10
  int myXBinning;
  int myYBinning;
  float myCCDTemp;
  float myFWTemp;
  float myXSpacing;
  float myYSpacing;
  char myComment[100];
  float myProjectionAltitude;
  int myDataOffset;
  int myCustomHeaderOffset;
  float myEmissionHeight;
  double myMean;
  double myStdDev;
  float myProjectionLon;

  // Site data
  float myHeight;
  float myLatitude;
  float myLongitude;
  char myName[30];
  char myAbbreviation[4];
};

#endif

//
// End of "$Id: Fl_TIFF_Image.H,v 1.1 2002/10/11 03:13:52 easysw Exp $".
//
