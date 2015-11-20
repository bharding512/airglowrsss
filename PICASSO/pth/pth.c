/* pth.c - a quick and dirty PICASSO TIFF Header reader especially for use
           with shell and Perl scripts.

   Ethan S. Miller (esmiller@uiuc.edu)   30 May 2007

   To build:  gcc -o pth pth.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define IMAGETYPE unsigned short
typedef unsigned short ushort;
typedef unsigned char uchar;

int main (int argc, char * argv[]) {
  char          header1[3];                     // More header data
  char          header2[88];
  uchar         *arrayx;
  long int      offset16 = 0;
  FILE          *fp;
  char          filename[512];
  char          csite[10];
  int           i = 0;

  IMAGETYPE *myData;
  unsigned short myRows;
  unsigned short myCols;
  IMAGETYPE myMin;
  IMAGETYPE myMax;
  struct tm myLocalTime;
  struct tm myUTTime;
  char myFilterName[12];
  float myExposureTime;
  int myGain;   /* Contains the gain times 10 */
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

  /* Site data */
  float myHeight;
  float myLatitude;
  float myLongitude;
  char myName[30];
  char myAbbreviation[4];

  if (argc == 1) {
    printf("PICASSO TIFF header information - version 1.0 - esmiller@uiuc.edu\n\n");
    printf("USAGE:  pth [filename1,[filename2,[filename3,...]]]\n\n");
    printf("All shell-expandable wildcards may be used.\n\n");
    return 0;
  }

  for (i = 1; i < argc; i++) {

    strncpy(filename,argv[i],512);
    fp = fopen(filename,"r");

    fseek(fp, 200L, SEEK_SET);
    fread(header1,1,3,fp);

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

    fread(&myHeight,4,1,fp);
    fread(&myName,1,30,fp);
    fread(&myAbbreviation,1,4,fp);
    fread(&myComment,1,100,fp);
    fread(header2,1,88,fp);

    (void)strncpy(csite,myName,8);
    csite[9] = 0;
    fclose(fp);

    printf("%13s (%1s) %03d.%04d %02d%02dUT (%s;%dx%d;%dx%d;%2.2f;%2.0f)\n",filename,csite,myUTTime.tm_yday+1,myUTTime.tm_year+1900,myUTTime.tm_hour,myUTTime.tm_min,myFilterName,myRows,myCols,myXBinning,myYBinning,myExposureTime,myCCDTemp);
  }

  return 0;
}
