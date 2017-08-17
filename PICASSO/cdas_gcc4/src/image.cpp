/***************************************************************************
                          image.cpp  -  description
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

#include <stdio.h>
#include <stdlib.h>

#include "system.h"
#include "image.h"

extern System mySystem;

Image::Image(){
  myCCDTemp = 0;
  myCols = 0;
  myComment = "";
  myCustomHeaderOffset = 200;
  myData = NULL;
  myNormalized = NULL;
  myDataOffset = 600;
  myEmissionHeight = 0;
  myExposureTime = 0;
  myFWTemp = 0;
  myFilterName = "";
  myGain = 10;
//  myLocalTime = NULL;
  myMax = 0;
  myMin = 0;
  myProjectionAltitude = 0;
  myProjectionLon = 0;
  myRows = 0;
//  mySite = "";
//  myUTTime = NULL;
  myXBinning = 1;
  myXSpacing = 0;
  myYBinning = 1;
  myYSpacing = 0;
}
Image::~Image(){
  if(myData != NULL)
    free(myData);
  if(myNormalized != NULL)
    free(myNormalized);
}
/** Writes the data to an extended tif file. */
//void Image::SaveTIF(ofstream &of){
void Image::SaveTIF(string filename){
  // Open the file for writing
  ofstream of;
  of.open(filename.c_str(), ios::out | ios::binary);
  
  // First write the default TIF header info
  BYTE space[4];
  space[0] = 'I';
  space[1] = 'I';
  space[2] = 42;
  space[3] = 0;
  of.write(space, 4);

  DWORD myDWORD;
  WORD myWORD;

  // The IFD offset
  myDWORD = 8;
  of.write((char *)&myDWORD, 4);

  // The IFD
  myWORD = 13;
  of.write((char *)&myWORD, 2);

  // Each entry in the IFD is 12 bytes long
  // Bytes 0-1  The tag that identifies the field
  // Bytes 2-3  The field type
  // Bytes 4-7  The number of values. Count of the indicated type
  // Bytes 8-11 The value offset, the file offset of the value for the field
  //
  // NOTES:
  // Types:
  //  1 = BYTE (8-bit unsigned int)
  //  2 = ASCII (8-bit byte that contains a 7-bit ASCII code
  //  3 = SHORT (16-bit (2-byte) unsigned int)
  //  4 = LONG (32-bit (4-byte) unsigned int)
  //  5 = RATIONAL (two LONGS)
  //
  // Values/Offset: If it is shorted than 4 bytes, the value can be put in the field

  // Byte offset 10
  // NewSubfileType
  myWORD = 0xFE;
  of.write((char *)&myWORD, 2);
  myWORD = 4;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 0;
  of.write((char *)&myDWORD, 4);

  // Byte offset 22
  // Image columns
  myWORD = 0x100;
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = myCols;
  of.write((char *)&myDWORD, 4);

  // Byte offset 34
  // Image rows
  myWORD = 0x101;
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = myRows;
  of.write((char *)&myDWORD, 4);

  // Byte offset 46
  // Image bits per pixel
  myWORD = 0x102;              
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 16;
  of.write((char *)&myDWORD, 4);

  // Byte offset 58
  // Compression code
  myWORD = 0x103;
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);

  // Byte offset 70
  // Photometric interpretation
  myWORD = 0x106;
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);

  // Byte offset 82
  // Strip offset (the Byte offset for the start of the data)
  myWORD = 0x111;
  of.write((char *)&myWORD, 2);
  myWORD = 4;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = myDataOffset;
  of.write((char *)&myDWORD, 4);

  // Byte offset 94
  // Samples per pixel
  myWORD = 0x115;
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);

  // Byte offset 106
  // Rows per strip
  myWORD = 0x116;
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = myRows;
  of.write((char *)&myDWORD, 4);

  // Byte offset 118
  // Strip byte count
  myWORD = 0x117;
  of.write((char *)&myWORD, 2);
  myWORD = 4;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = myRows * myCols * 2;
  of.write((char *)&myDWORD, 4);

  // Byte offset 130
  // X-resolution
  myWORD = 0x11A;
  of.write((char *)&myWORD, 2);
  myWORD = 5;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 100;
  of.write((char *)&myDWORD, 4);

  // Byte offset 142
  // Y-resolution
  myWORD = 0x11B;
  of.write((char *)&myWORD, 2);
  myWORD = 5;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 100;
  of.write((char *)&myDWORD, 4);

  // Byte offset 154
  // Resolution Unit
  myWORD = 0x128;
  of.write((char *)&myWORD, 2);
  myWORD = 3;
  of.write((char *)&myWORD,2);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);
  myDWORD = 1;
  of.write((char *)&myDWORD, 4);

  // Byte offset 166
  // That was the last record
  myDWORD = 0;
  of.write((char *)&myDWORD, 4);

  // Now get to the beginning of where the custom header is going to go
  of.seekp(myCustomHeaderOffset,ios::beg);

  // Save the custom header information
  SaveCustomHeader(of);

  // Go to the data offset
  of.seekp(myDataOffset, ios::beg);

  // Write the data
  of.write((char *)myData, myRows * myCols * 2);

  // close the file
  of.close();
}
/** No descriptions */
void Image::SaveCustomHeader(ofstream &of){
  // Saves the custom TIF header
  float myFloat;
  BYTE myTemp[3];
  myTemp[0] = 'C';
  myTemp[1] = 'D';
  myTemp[2] = 'F';

  of.write((char *)myTemp, 3);
  of.write((char *)&myRows, 2);
  of.write((char *)&myCols, 2);
  of.write((char *)&myMin, 2);
  of.write((char *)&myMax, 2);
  of.write((char *)&myLocalTime, 36);
  of.write((char *)&myUTTime, 36);
  of.write((char *)myFilterName.c_str(), 12);
  of.write((char *)&myEmissionHeight,4);
  of.write((char *)&myExposureTime,4);
  of.write((char *)&myGain,1);
  of.write((char *)&myXBinning,1);
  of.write((char *)&myYBinning,1);
  of.write((char *)&myCCDTemp,4);
  of.write((char *)&myFWTemp,4);
  of.write((char *)&myXSpacing,4);
  of.write((char *)&myYSpacing,4);
  of.write((char *)&myCalibrationCoefs.a0, 8);
  of.write((char *)&myCalibrationCoefs.a1, 8);
  of.write((char *)&myCalibrationCoefs.a2, 8);
  of.write((char *)&myCalibrationCoefs.b0, 8);
  of.write((char *)&myCalibrationCoefs.b1, 8);
  of.write((char *)&myCalibrationCoefs.b2, 8);
  myFloat = mySite.GetLatitude();
  of.write((char *)&myFloat, 4);
  myFloat = mySite.GetLongitude();
  of.write((char *)&myFloat, 4);
  myFloat = mySite.GetHeight();
  of.write((char *)&myFloat, 4);
  of.write((char *)mySite.GetName().c_str(), 30);
  of.write((char *)myComment.c_str(), 100);
  of.write((char *)&myCalibrationCoefs.myIsAllSky, 1);
  of.write((char *)&myCalibrationCoefs.myCenterAz, 8);
  of.write((char *)&myCalibrationCoefs.myCenterEl, 8);
  of.write((char *)&myProjectionAltitude, 4);
  of.write((char *)&myProjectionLon, 4);
}
/** Reads a TIF file into the Image structure */
void Image::LoadTIF(ifstream &is){
  char myTempChar[200];
  char myTemp[200];
  char myType[3];
  char myFilter[12];
  char myTempSiteName[30];
  char myTempComment[100];
  float myTempFloat;
  bool isGMS;
  string myString;

  is.read(myTemp, 200);
  is.read(myType, 3);

  if(!strncmp(myType,"CDF",3) == 0)
    isGMS = true;
  else
    isGMS = false;
    
//  is.read(myInt, 2);
  is.read((char *)&myRows, 2);
  is.read((char *)&myCols, 2);
  is.read((char *)&myMin, 2);
  is.read((char *)&myMax, 2);
  is.read((char *)&myLocalTime, 36);
  is.read((char *)&myUTTime, 36);
  is.read((char *)&myFilter, 12);//myTempChar,12);
  myFilterName = myFilter;
  is.read((char *)&myEmissionHeight, 4);
  is.read((char *)&myExposureTime, 4);
  is.read((char *)&myGain, 1);
  is.read((char *)&myXBinning, 1);
  is.read((char *)&myYBinning, 1);
  is.read((char *)&myCCDTemp, 4);
  is.read((char *)&myFWTemp, 4);
  is.read((char *)&myXSpacing, 4);
  is.read((char *)&myYSpacing, 4);
  is.read((char *)&myCalibrationCoefs.a0, 8);
  is.read((char *)&myCalibrationCoefs.a1, 8);
  is.read((char *)&myCalibrationCoefs.a2, 8);
  is.read((char *)&myCalibrationCoefs.b0, 8);
  is.read((char *)&myCalibrationCoefs.b1, 8);
  is.read((char *)&myCalibrationCoefs.b2, 8);
  is.read((char *)&myTempFloat, 4);
  mySite.SetLatitude(myTempFloat);
  is.read((char *)&myTempFloat, 4);
  if(isGMS)
  {
    // GMS files are referenced from west not east
    myTempFloat = 360 - myTempFloat;
  }
  mySite.SetLongitude(myTempFloat);
  is.read((char *)&myTempFloat, 4);
  mySite.SetHeight(myTempFloat);
  is.read(myTempSiteName,30);
  mySite.SetSiteName(myTempSiteName);
  is.read(myTempComment,100);
  myComment = myTempChar;

  if(!isGMS)
  {
    // Read in new parameters added from original format
    is.read((char *)&myCalibrationCoefs.myIsAllSky, 1);
    is.read((char *)&myCalibrationCoefs.myCenterAz, 8);
    is.read((char *)&myCalibrationCoefs.myCenterEl, 8);
    is.read((char *)&myProjectionAltitude, 4);
    is.read((char *)&myProjectionLon, 4);
    is.read(myTempChar, 63);
  }
  else
    is.read(myTempChar, 88);

  // Set the calibration structure information
  myCalibrationCoefs.SetCoefficients();

  // Clear the data pointers if needed
  if(myData != NULL)
    free(myData);
  myData = (IMAGETYPE*) new BYTE[myRows*myCols*sizeof(IMAGETYPE)];
  is.read((char *)myData, myRows*myCols*sizeof(IMAGETYPE));  
}
/** No descriptions */
IMAGETYPE* Image::GetDataPointer(){
  return myData;
}
/** No descriptions */
void Image::SetRows(int row){
  SetRowsCols(row, myCols);
}
/** No descriptions */
void Image::SetCols(int col){
  SetRowsCols(myRows, col);
  
}
/** No descriptions */
void Image::SetRowsCols(int rows, int cols){
  myRows = rows;
  myCols = cols;

  if(myData != NULL)
    free(myData);

  // allocate the array data
  myData = (IMAGETYPE *) new BYTE[sizeof(IMAGETYPE)*myRows*myCols];
}
/** No descriptions */
void Image::SetCCDTemp(float temp){
  myCCDTemp = temp;
}
/** No descriptions */
void Image::SetBinning(int x, int y){
  myXBinning = x;
  myYBinning = y;
}
/** No descriptions */
void Image::SetExposureTime(float expTime){
  myExposureTime = expTime;
}
/** No descriptions */
void Image::SetSite(Site *site){
  mySite.SetHeight(site->GetHeight());
  mySite.SetLatitude(site->GetLatitude());
  mySite.SetLongitude(site->GetLongitude());
  mySite.SetSiteName(site->GetName());
  mySite.SetAbbreviation(site->GetAbbreviation());
}
/** No descriptions */
void Image::SetTimesNow(){
  time_t now;
  struct tm* temp_tm;

  // Get the current time
  time(&now);

  temp_tm = gmtime(&now);
  myUTTime = *temp_tm;
  temp_tm = localtime(&now);
  myLocalTime = *temp_tm;  
}
/** Calculates the mean of the intensity data in the myData array
    The mean is stored in the variable myMean */
void Image::CalcMean(){
  long pix;
  double sum;
  double size;
  IMAGETYPE *ptr;
  IMAGETYPE curPix;

  if(myData == NULL)
    return;

  // initialize the counters
  sum = 0.;
  size = myRows * myCols;
  ptr = myData;

  for(pix = 0; pix < myRows * myCols; pix++)
  {
    curPix = (*ptr++);
    if(curPix != 0)
      sum += (double) curPix;
  }

  myMean = sum/size;
}
/** Calculates the standard deviation of the intensity data in the myData
    array.  The result is stored in the variable myStdDev */
void Image::CalcStdDev(){
  long pix;
  double sum;
  double sum2;
  long size;
  IMAGETYPE *ptr;
  IMAGETYPE curPix;

  if(myData == NULL)
    return;

  // Initialize the counters
  sum = 0.;
  sum2 = 0.;
  size = myRows * myCols;
  ptr = myData;

  for(pix = 0; pix < myRows * myCols; pix++)
  {
    curPix = *ptr;
    if(curPix != 0)
    {
      sum += (double) (((double) curPix) * (curPix));
      sum2 += (double) curPix;
    }
    ptr++;
  }

  sum /= size;
  sum2 /= size;

  if(sum < sum2 * sum2)
    myStdDev = 0.;
  else
    myStdDev = sqrt(sum - sum2 * sum2); 
}
/** Sets the myMin and myMax values for the purpose of autocontrasting an image.  The values of myMin and myMax are set to be +/- one stddev from the mean.

Autocontrast assumes that myMean and myStdDev have been set by calls to CalcMean() and CalcStdDev() respectively. */
void Image::AutoContrast(){
  if(myMean == 0 and myStdDev == 0)
  {
    // Seems like the mean and stddev havn't been caluculated
    CalcMean();
    CalcStdDev();
  }
  
  myMin = (IMAGETYPE) (myMean - myStdDev);
  myMax = (IMAGETYPE) (myMean + myStdDev);

  if(myMin < 0)
    myMin = 0;
  if(myMax > 65530)
    myMax = 65530;
}
/** No descriptions */
void Image::SetMean(double mean){
  myMean = mean;
}
/** No descriptions */
void Image::SetStdDev(double dev){
  myStdDev = dev;
}

/** Saves a PNG version of the image.  The image is scaled to an 8-bit image based on myMin and myMax. */
void Image::SavePNG(string filename){

  filename += ".png";

  // grab pointers to the rows
  png_byte **row_pointers;

  row_pointers = (png_byte **) new BYTE[sizeof(png_byte *)*myRows];
  for(int i = 0; i < myRows; i++)
  {
    row_pointers[i] = (png_byte *)&myNormalized[i*myCols];
  }

  FILE *fp = fopen(filename.c_str(), "wb");
  if(!fp)
  {
    // The was an error
    string str("ERROR: opening png file: ");
    str = str + filename;
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    free(row_pointers);
    return;
  }

  png_structp png_ptr;

  // Initialize the png
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if(!png_ptr)
  {
    // The was an error
    string str("ERROR: Initializing PNG pointer");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    free(row_pointers);
    return;
  }

  png_infop info_ptr;

  // Initialize the info pointer
  info_ptr = png_create_info_struct(png_ptr);

  if(!info_ptr)
  {
    // The was an error
    string str("ERROR: Initializing info pointer");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    free(row_pointers);
    return;
  }

  // Set the rows
  png_set_rows(png_ptr, info_ptr, row_pointers);

  if(setjmp(png_jmpbuf(png_ptr)))
  {
    // The was an error
    string str("ERROR: Error during init_io");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    free(row_pointers);
    return;
  }

  png_init_io(png_ptr, fp);

  // write header
  if(setjmp(png_jmpbuf(png_ptr)))
  {
    // The was an error
    string str("ERROR: Writing PNG header");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    free(row_pointers);
    return;
  }

  png_byte color_type;
  png_byte bit_depth;

  color_type = PNG_COLOR_TYPE_GRAY;
  bit_depth = 8;

  png_set_IHDR(png_ptr, info_ptr, myCols, myRows,
    bit_depth, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
    PNG_FILTER_TYPE_BASE);

  png_write_info(png_ptr, info_ptr);

  // write bytes
  if(setjmp(png_jmpbuf(png_ptr)))
  {
    // The was an error
    string str("ERROR: Writing PNG bytes");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    free(row_pointers);
    return;
  }

//  png_write_image(png_ptr, row_pointers);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

  // end write
  if(setjmp(png_jmpbuf(png_ptr)))
  {
    // The was an error
    string str("ERROR: End of write");
    mySystem.myConfiguration.myLogger.WriteErrorString(str);
    free(row_pointers);
    return;
  }

  png_write_end(png_ptr, info_ptr);

  png_destroy_write_struct(&png_ptr, &info_ptr);

  fclose(fp);
  
  free(row_pointers);
}
/** Normalizes the data between myMin and myMax to 0-255 */
void Image::Normalize(){
  IMAGETYPE pix;
  
  if(myNormalized != NULL)
    free(myNormalized);

  // allocate the memory
  myNormalized = (png_byte *) new BYTE[myRows * myCols * sizeof(png_byte)];

  int i, j;

  if(myData == NULL)
  {
    // no data, so fill with 0s
    for(i = 0; i < myRows; i++)
      for(j = 0; j < myCols; j++)
        myNormalized[i*myCols + j] = 0;
  }
  else
  {
    // valid data, so normalize it to the myMin-myMax range
    for(i = 0; i < myRows; i++)
    {
      for(j = 0; j < myCols; j++)
      {
        pix = myData[i*myCols + j];
        if(pix >= myMean+3*myStdDev)
          myNormalized[i*myCols + j] = (png_byte) 0;  // For large outliers
        else if(pix >= myMax)
          myNormalized[i*myCols + j] = (png_byte) 255;
        else if(pix <= myMin)
          myNormalized[i*myCols + j] = (png_byte) 0;
        else                                     
          myNormalized[i*myCols + j] = (png_byte)((double)(pix - myMin)/(myMax-myMin)*255);
      }
    }
  }
}
/** Sets the filter name. */
void Image::SetFilterName(string name){
  myFilterName = name;
}
/** Normalizes the data between myMin and myMax to 0-255 and subtracts the dark image (if dark == true) */
void Image::Normalize(bool dark){
  IMAGETYPE pix, darkpix;

  if(myNormalized != NULL)
    free(myNormalized);

  // allocate the memory
  myNormalized = (png_byte *) new BYTE[myRows * myCols * sizeof(png_byte)];

  int i, j;

  IMAGETYPE *DarkData;
  DarkData = mySystem.myConfiguration.myCamera->GetDarkDataPointer();

  if(myData == NULL || DarkData == NULL)
  {
    // no data, so fill with 0s
    for(i = 0; i < myRows; i++)
      for(j = 0; j < myCols; j++)
        myNormalized[i*myCols + j] = 0;
  }
  else
  {
    // valid data, so normalize it to the myMin-myMax range
    for(i = 0; i < myRows; i++)
    {
      for(j = 0; j < myCols; j++)
      {
        pix = myData[i*myCols + j];
        if(dark)
          darkpix = DarkData[i*myCols + j];
        else
          darkpix = 0;

        // Check if the dark image pixel is greater than the image pixel value
        if(pix < darkpix)
        	pix = 0;
        else
        	pix = pix - darkpix;
          
        if(pix >= myMax)
          myNormalized[i*myCols + j] = (png_byte) 255;
        else if(pix <= myMin)
          myNormalized[i*myCols + j] = (png_byte) 0;
        else
          myNormalized[i*myCols + j] = (png_byte)((double)(pix - myMin)/(myMax-myMin)*255);
      }
    }
  }
}
/** Sets the CCD gain. */
void Image::SetGain(int gain){
	myGain = gain;
}
