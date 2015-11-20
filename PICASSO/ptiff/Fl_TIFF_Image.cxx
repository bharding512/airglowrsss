//
// "$Id: Fl_TIFF_Image.cxx,v 1.11 2006/06/29 15:33:53 ethan Exp $".
//
// TIFF image routines for flphoto.
//
// Copyright 2006 by Ethan Miller.
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
// Contents:
//
//   Fl_TIFF_Image::Fl_TIFF_Image() - Load a TIFF image...
//   Fl_TIFF_Image::Fl_TIFF_Image() - Load a TIFF image...
//   Fl_TIFF_Image::check()        - Try loading the named file.
//   Fl_TIFF_Image::load_image()   - Load the image from the file.
//

//
// Include necessary header files...
//

#include <FL/Fl.H>
#include "Fl_TIFF_Image.H"
#include <FL/Fl_JPEG_Image.H>
#include <stdlib.h>
#include "string.h"
#include <iostream>
#include <sys/stat.h>
#include <time.h>

//
// 'Fl_TIFF_Image::Fl_TIFF_Image()' - Load a TIFF image...
//

Fl_TIFF_Image::Fl_TIFF_Image(const char *name)	// I - Filename to read from
  : Fl_RGB_Image(0,0,0) {
  FILE          *fpd = NULL;
  FILE		*fp;				// File pointer
  struct stat   statbufI;
  struct stat   statbufD;
  double        mindif = 1e6;
  int           mindifi = -1;
  int           i = 0, r = 0;
  char          fname[512]; 
  char          path[512];
  char          *tempStr;

  strncpy(path,name,512);
  tempStr = rindex(path,'/');
  tempStr[0] = '\0'; 
  stat(name,&statbufI);
  for (i = 0; i < 100; i++) {
    sprintf(fname,"%s/DARK_%04i.tif",path,i);
    r = stat(fname,&statbufD);
    if ((r == 0) & (fabs(difftime(statbufI.st_mtime,statbufD.st_mtime)) < mindif)) {
      mindif = fabs(difftime(statbufI.st_mtime,statbufD.st_mtime));
      mindifi = i;
    }
  }

  if (mindifi != -1) {
    sprintf(fname,"%s/DARK_%04i.tif",path,mindifi);
    fpd = fopen(fname,"rb");
  } else {
    fpd = NULL;
  }

  // Open the image file...
  fp = fopen(name, "rb");
  if (fp == NULL)
    return;

  // Load the image...
  load_image(fp);
  fclose(fp);
  if (fpd != NULL) {
    dark_image(fpd);
    fclose(fpd);
  }
}


//
// 'Fl_TIFF_Image::Fl_TIFF_Image()' - Load a TIFF image...
//

Fl_TIFF_Image::Fl_TIFF_Image(FILE *fp, FILE *fpd)	// I - File stream to read from
  : Fl_RGB_Image(0,0,0) {
  // Load the image...
  load_image(fp);
  if (fpd != NULL) {
    dark_image(fpd);
  }
}

//
// 'Fl_TIFF_Image::check()' - Try loading the named file.
//

Fl_Image *					// O - New image or NULL
Fl_TIFF_Image::check(const char *name,
		    uchar      *header,		
		     int        headerlen)	// I - Number of header bytes

{
  FILE		*fp;				// Image file
  FILE          *fpd;
  char          header0[4];
  Fl_Image	*img;				// New image
  struct stat   statbufI;
  struct stat   statbufD;
  double        mindif = 1e6;
  long int      r = 0;
  int           mindifi = -1;
  int           i = 0;
  char          fname[512];
  char          path[512];
  char          *tempStr; 

  // Open the file to see if this is a TIFF image...
  fp = fopen(name, "rb");
  if (fp == NULL) {
    return (0);
  }

  strncpy(path,name,512);
  tempStr = rindex(path,'/');
  tempStr[0] = '\0'; 
  stat(name,&statbufI);
  for (i = 0; i < 100; i++) {
    sprintf(fname,"%s/DARK_%04i.tif",path,i);
    r = stat(fname,&statbufD);
    if ((r == 0) & (fabs(difftime(statbufI.st_mtime,statbufD.st_mtime)) < mindif)) {
      mindif = fabs(difftime(statbufI.st_mtime,statbufD.st_mtime));
      mindifi = i;
    }
  }

  if (mindifi != -1) {
    sprintf(fname,"%s/DARK_%04i.tif",path,mindifi);
    fpd = fopen(fname,"rb");
  } else {
    fpd = NULL;
  }

  // Look for TIFF header
  fread(header0,1,4,fp);
  if (memcmp(header0,"II*\0",4)==0) {
    img = new Fl_TIFF_Image(fp,fpd);
  } else {
    img = 0;
  }
  // Close the file and return the image we loaded, if any...
  if (fpd != NULL) {
    fclose(fpd);
  }
  fclose(fp);
  return (img);
}


//
// 'Fl_TIFF_Image::load_image()' - Load the image from the file.
//

void
Fl_TIFF_Image::load_image(FILE *fp)		// I - File to load
{
  int		left,				// Left position of image
		top,				// Top position of image
		rotation;			// Rotation of image
  long int      tempLong = 0;
  char		header1[3];			// More header data
  char          header2[88];
  bool          isGMS;
  uchar         *arrayx;
  long int      offset16 = 0;

  fseek(fp, 200L, SEEK_SET);
  fread(header1,1,3,fp);
  if(!memcmp(header1,"CDF",3) == 0)
    isGMS = true;
  else
    isGMS = false;
  fseek(fp, 203L, SEEK_SET);
  fread(&myRows,2,1,fp);
  fread(&myCols,2,1,fp);
  fread(&myMin,2,1,fp);
  fread(&myMax,2,1,fp);
  fread(&myLocalTime,1,36,fp);
  fread(&myUTTime,1,36,fp);
  fread(&myFilterName,1,12,fp);
  fread(&myEmissionHeight,4,1,fp);
  fread(&myExposureTime,4,1,fp);
  fread(&myGain,1,1,fp);
  fread(&myXBinning,1,1,fp);
  fread(&myYBinning,1,1,fp);
  fread(&myCCDTemp,4,1,fp);
  fread(&myFWTemp,4,1,fp);
  fread(&myXSpacing,4,1,fp);
  fread(&myYSpacing,4,1,fp);
  fread(header2,1,48,fp);   //Calibration coefficients
  fread(&myLatitude,4,1,fp);
  fread(&myLongitude,4,1,fp);
  if (isGMS) {
    myLongitude = 360 - myLongitude;
  }
  fread(&myHeight,4,1,fp);
  fread(&myName,1,30,fp);
  fread(&myComment,1,100,fp);
  fread(header2,1,88,fp);

  // Get the image size and orientation...
  rotation = 0;
  left     = 0;
  top      = 0;

  // Allocate and initialize...
  w(myCols);
  h(myRows);
  d(3);

  alloc_array = 1;
  offset16 = ((long int) w())*((long int) h());
  arrayx = new uchar[offset16*5];
  // Read the image...
  fread(&arrayx[0]+offset16,1,offset16*2,fp);
  for (int i = 0; i < offset16; i++) {
    tempLong = ((arrayx[i*2+offset16]+256*arrayx[i*2+offset16+1]) - myMin)*256/(myMax - myMin);
    if (tempLong > 255) {
      arrayx[i] = 255;
    } else if (tempLong < 0) {
      arrayx[i] = 0;
    } else {
      arrayx[i] = (char) tempLong;
    }
    arrayx[i+offset16*3] = 0;
    arrayx[i+offset16*4] = 0;
  }
  array = &arrayx[0];
  d(1);
}

void
Fl_TIFF_Image::dark_image(FILE *fp)		// I - File to load
{
  int		left,				// Left position of image
		top,				// Top position of image
		rotation;			// Rotation of image
  long int      tempLong = 0;
  char		header1[3];			// More header data
  char          header2[88];
  bool          isGMS;
  long int      offset16 = 0;
  uchar         *arrayx;

  fseek(fp, 200L, SEEK_SET);
  fread(header1,1,3,fp);
  if(!memcmp(header1,"CDF",3) == 0)
    isGMS = true;
  else
    isGMS = false;
  fseek(fp, 203L, SEEK_SET);
  fread(&myRows,2,1,fp);
  fread(&myCols,2,1,fp);
  fread(&myMin,2,1,fp);
  fread(&myMax,2,1,fp);
  fread(&myLocalTime,1,36,fp);
  fread(&myUTTime,1,36,fp);
  fread(&myFilterName,1,12,fp);
  fread(&myEmissionHeight,4,1,fp);
  fread(&myExposureTime,4,1,fp);
  fread(&myGain,1,1,fp);
  fread(&myXBinning,1,1,fp);
  fread(&myYBinning,1,1,fp);
  fread(&myCCDTemp,4,1,fp);
  fread(&myFWTemp,4,1,fp);
  fread(&myXSpacing,4,1,fp);
  fread(&myYSpacing,4,1,fp);
  fread(header2,1,48,fp);   //Calibration coefficients
  fread(&myLatitude,4,1,fp);
  fread(&myLongitude,4,1,fp);
  if (isGMS) {
    myLongitude = 360 - myLongitude;
  }
  fread(&myHeight,4,1,fp);
  fread(&myName,1,30,fp);
  fread(&myComment,1,100,fp);
  fread(header2,1,88,fp);

  // Get the image size and orientation...
  rotation = 0;
  left     = 0;
  top      = 0;

  offset16 = ((long int) myRows)*((long int) myCols);
  arrayx = (uchar *) &array[0]+offset16*3;
  fread(arrayx,1,offset16*2,fp);
}

//
// End of "$Id: Fl_TIFF_Image.cxx,v 1.11 2006/06/29 15:33:53 ethan Exp $".
//
