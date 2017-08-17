/***************************************************************************
                          cameraapogee  -  description
                             -------------------
    begin                : Fri Aug 17 2007
    copyright            : (C) 2007 by Jonathan Makela
    email                : jmakela@uiuc.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "cameraapogee.h"
#include "system.h"
#include <sstream>

extern System mySystem;

unsigned short hextoi(char* instr);

#define MAX_CCD_BUFFERS  1000
#define MAXTOTALCOLUMNS 16383
#define MAXTOTALROWS 	16383

int verbose=0;

/* Error codes returned from config_load */
const long CCD_OPEN_NOERR   = 0;          // No error detected
const long CCD_OPEN_CFGNAME = 1;          // No config file specified
const long CCD_OPEN_CFGDATA = 2;          // Config missing or missing required data
const long CCD_OPEN_LOOPTST = 3;          // Loopback test failed, no camera found
const long CCD_OPEN_ALLOC   = 4;          // Memory alloc failed - system error

void trimstr(char* s);

typedef struct {
	unsigned short *pixels;
	int            size;
	short          xdim;
	short          ydim;
	short          zdim;
	short          xbin;
	short          ybin;
	short          type;
	char           name[64];
	int            shmid;
	size_t         shmsize;
	char           *shmem;
} CCD_FRAME;

CCD_FRAME CCD_Frame[MAX_CCD_BUFFERS];

bool CfgGet ( FILE* inifp,
			 char* inisect,
			 char* iniparm,
			 char* retbuff,
			 short bufflen,
			 short* parmlen);
		
CameraAPOGEE::CameraAPOGEE()
{
	// Initialize variables
	myCameraOK = false;
	myCCDSetTemp = 0;
	myCCDTemp = 25;
	myBinX = 1;
	myBinY = 1;
	myTop = 1;
	myLeft = 1;
	myLeftOffset = 0;
	myTopOffset = 0;
	myWidth = 1;
	myHeight = 1;
	myExposureTime = 0.1;
	myDarkData = NULL;
	myAutoTempRegulation = false;
	myUSBgain = 1.0;
	myDark = false;
	myCount = 0;
	
	if(!InitCamera())
		// There was an error
		mySystem.myConfiguration.myCameraType = ERROR;
	else
		// Alls well
		mySystem.myConfiguration.myCameraType = APOGEE;
}


CameraAPOGEE::~CameraAPOGEE()
{
	if(myDarkData != NULL)
		free(myDarkData);
	
	// Stop the exposure, if it is occuring
	if(mySystem.myIsExposing)
		KillExposure();
	
	ShutDownCamera();
	
	mySystem.myConfiguration.myCameraType = ERROR;
}

/*!
    \fn CameraAPOGEE::InitCamera
 */
bool CameraAPOGEE::InitCamera()
{
	/*	Parse camera configuration  file */
	char *hdir;
	char cfgname[256];
	hdir = getenv("HOME");
	sprintf(cfgname,"%s/.apccd.ini",hdir);
	long ret = config_load( cfgname, -1, -1 );
	if ( ret != 0 ) {
		switch ( ret )
		{
			case CCD_OPEN_CFGNAME:
				mySystem.myConfiguration.myLogger.WriteErrorString("cameraapogee: No config file specified.");
				break;
			case CCD_OPEN_CFGDATA:
				mySystem.myConfiguration.myLogger.WriteErrorString("cameraapogee: Config file missing or missing required data.");
				break;
			case CCD_OPEN_LOOPTST:
				mySystem.myConfiguration.myLogger.WriteErrorString("cameraapogee: Loopback test failed, no camera found");
				break;
			case CCD_OPEN_ALLOC:
				mySystem.myConfiguration.myLogger.WriteErrorString("cameraapogee: Memory allocation failed - system error");
				break;
		}
		printf("Cannot load configuration\n");
		myCameraOK = false;
		return myCameraOK;
	}

	/*	Assume only 1 camera, and it is number 0 */
	cameraDevice->InitDriver(0);
	
	/*	Do a system reset to ensure known state, flushing enabled etc */
	cameraDevice->Reset();
	
	/*	Parse config file for camera parameters */
	cameraDevice->Flush();
	
	// Set to slow readouts	
	cameraDevice->write_FastReadout(false);

	// Set to long cable
	cameraDevice->write_LongCable(true);
	
	// Read the initial temperature
	ReadTemperatureStatus();
	
	// We have sucessfully setup the camera
	myCameraOK = true;
	
	return myCameraOK;
}

/*!
    \fn CameraAPOGEE::ReadTemperatureStatus()
 */
void CameraAPOGEE::ReadTemperatureStatus()
{	
	myCCDTemp = cameraDevice->read_Temperature();
	
	return;
}

/*!
    \fn CameraAPOGEE::KillExposure()
 */
void CameraAPOGEE::KillExposure()
{
  cameraDevice->Reset();
}

/*!
    \fn CameraAPOGEE::ShutDownCamera
 */
int CameraAPOGEE::ShutDownCamera()
{
	// Warm the CCD up to 0
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Warming up sensor to 0 C");
	cameraDevice->write_CoolerMode(1);
	cameraDevice->write_CoolerSetPoint(25.0);
	while(myCCDTemp < 0)
	{
		ReadTemperatureStatus();
		cout << myCCDTemp << endl;
		sleep(10);
	}
	cameraDevice->write_CoolerMode(0);
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Sensor warmed up successfully");
		
	myCameraOK = false;
	
	return 0;
}

/*!
    \fn CameraAPOGEE::SetCCDTemp()
 */
void CameraAPOGEE::SetCCDTemp(double CCDTemp)
{
	// Place some error bounds on it
	if(CCDTemp > 20.)
		CCDTemp = 20.;
	
	if(myCameraOK)
	{
		// Set the temperature setting
		fflush(NULL);
		printf("Turning cooler on\n");
		cameraDevice->write_CoolerMode(1);
		cameraDevice->write_CoolerSetPoint(CCDTemp);
		printf("Cooler setpoint %f\n", cameraDevice->read_CoolerSetPoint());
		myCCDSetTemp = cameraDevice->read_CoolerSetPoint();
		isTempRegulating = true;
	}
}

/*!
    \fn CameraAPOGEE::GetTemperatureRegulation()
 */
bool CameraAPOGEE::GetTemperatureRegulation()
{
	return isTempRegulating;
}


/*!
    \fn CameraAPOGEE::GetCCDTemp()
 */
double CameraAPOGEE::GetCCDTemp()
{
	return myCCDTemp;
}


/*!
    \fn CameraAPOGEE::GetExposureTimeLeft()
 */
int CameraAPOGEE::GetExposureTimeLeft()
{
	return myExposureTimeLeft;
}

/*!
    \fn CameraAPOGEE::Read(ifstream &is)
 */
void CameraAPOGEE::Read(ifstream &is)
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
	is >> myAutoTempRegulation;}


/*!
    \fn CameraAPOGEE::Write(ofstream &os)
 */
void CameraAPOGEE::Write(ofstream &os)
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
    \fn CameraAPOGEE::WriteXML(ofstream &os)
 */
void CameraAPOGEE::WriteXML(ofstream &os)
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


/*!
    \fn CameraAPOGEE::SetImageSize(int width, int height)
 */
void CameraAPOGEE::SetImageSize(int width, int height)
{
	myWidth = (long) width;
	myHeight = (long) height;
}


/*!
    \fn CameraAPOGEE::SetBinnint(int X, int Y)
 */
void CameraAPOGEE::SetBinnint(int X, int Y)
{
	myBinX = X;
	myBinY = Y;
}


/*!
    \fn CameraAPOGEE::SetOffsets(int top, int left)
 */
void CameraAPOGEE::SetOffsets(int top, int left)
{
	myTop = top;
	myLeft = left;
}


/*!
    \fn CameraAPOGEE::SetExpTime(float exp)
 */
void CameraAPOGEE::SetExpTime(float exp)
{
	myExposureTime = exp;
}


/*!
    \fn CameraAPOGEE::SetBinX(int X)
 */
void CameraAPOGEE::SetBinX(int X)
{
	myBinX = X;
}


/*!
    \fn CameraAPOGEE::SetBinY(int Y)
 */
void CameraAPOGEE::SetBinY(int Y)
{
	myBinY = Y;
}


/*!
    \fn CameraAPOGEE::SetHeight(int height)
 */
void CameraAPOGEE::SetHeight(int height)
{
	myHeight = height;
}


/*!
    \fn CameraAPOGEE::SetWidth(int width)
 */
void CameraAPOGEE::SetWidth(int width)
{
	myWidth = width;
}


/*!
    \fn CameraAPOGEE::SetTop(int top)
 */
void CameraAPOGEE::SetTop(int top)
{
	myTop = top;
}


/*!
    \fn CameraAPOGEE::SetLeft(int left);
 */
void CameraAPOGEE::SetLeft(int left)
{
	myLeft = left;
}


/*!
    \fn CameraAPOGEE::GetWidth()
 */
int CameraAPOGEE::GetWidth()
{
	return myWidth;
}


/*!
    \fn CameraAPOGEE::GetHeight()
 */
int CameraAPOGEE::GetHeight()
{
	return myHeight;
}


/*!
    \fn CameraAPOGEE::GetTop()
 */
int CameraAPOGEE::GetTop()
{
	return myTop;
}


/*!
    \fn CameraAPOGEE::GetLeft()
 */
int CameraAPOGEE::GetLeft()
{
	return myLeft;
}


/*!
    \fn CameraAPOGEE::GetBinX();
 */
int CameraAPOGEE::GetBinX()
{
	return myBinX;
}


/*!
    \fn CameraAPOGEE::GetBinY()
 */
int CameraAPOGEE::GetBinY()
{
	return myBinY;
}


/*!
    \fn CameraAPOGEE::GetCCDSetTemp()
 */
double CameraAPOGEE::GetCCDSetTemp()
{
	return myCCDSetTemp;
}


/*!
    \fn CameraAPOGEE::GetExposureTime()
 */
float CameraAPOGEE::GetExposureTime()
{
	return myExposureTime;
}


/*!
    \fn CameraAPOGEE::SetDark(bool dark)
 */
void CameraAPOGEE::SetDark(bool dark)
{
	myDark = dark;
}

////// *********

/*!
    \fn CameraAPOGEE::CaptureImage(Image *myImage, string filename, bool singleimage)
 */
bool CameraAPOGEE::CaptureImage(Image *myImage, string filename, bool singleimage)
{
	// First check to make sure that the camera is initialized
	if(!myCameraOK)
	{
		// Camera isn't initialized, so try to initialize it
		if(!InitCamera())
		{
			// There was an error
			string str("Error initializing APOGEE camera.");
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
	
	// Save the header information
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
    \fn CameraAPOGEE::DoExposureCollect()
 */
bool CameraAPOGEE::DoExposureCollect()
{
	Camera_Status myStatus;
	
	// Setup Image dimensions
	cameraDevice->m_StartX = myLeft+myLeftOffset;
	cameraDevice->m_StartY = myTop+myTopOffset;
	cameraDevice->m_NumX = myWidth;
	cameraDevice->m_NumY = myHeight;
		
	// Set Binning
//	cameraDevice->m_BinX = myBinX;
//	cameraDevice->m_BinY = myBinY;
	
	if(myDark) {
	  cameraDevice->Expose(myExposureTime, 0);
	}
	else {
	  cameraDevice->Expose(myExposureTime, 1);
	}
	
	do
	{
		sleep(1);
		myStatus = cameraDevice->read_Status();
	}while(myStatus != 0);
		
	// Read the CCD
	unsigned short *image;
	int bnum;
	short nx, ny;
	int row;
	double sum1, sum2;
	sum1 = sum2 = 0.;
	char im_name[256];
	
	sprintf(im_name,"im%03d", myCount); 
	myCount++;
	
	// Readout the image and save in a named buffer (tempobs);
	myStatus = cameraDevice->BufferImage(im_name);
	
	// Use the libccd routine to find the corresponding buffer index
	bnum = CCD_locate_buffernum(im_name);
	
	// Obtain the memory address of the actual image data, and x,y dimensions
	image = CCD_Frame[bnum].pixels;
	nx = CCD_Frame[bnum].xdim;
	ny = CCD_Frame[bnum].ydim;
	
	for(row = 0; row < myRows; row++)
	{		
		// Copy the data over into the image data pointer
		for(int col = 0; col < myCols; col++)
		{
			myCameraData[row*myCols + col] = (IMAGETYPE) image[row*myCols + col];
			sum1 += (double) myCameraData[row*myCols + col];
			sum2 += (double) (((double)myCameraData[row*myCols + col]) * myCameraData[row*myCols + col]);
		}
	}
	
	delete image;
	
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
    \fn CameraAPOGEE::SetDarkData(IMAGETYPE *data, int width, int height)
 */
void CameraAPOGEE::SetDarkData(IMAGETYPE *data, int width, int height)
{
	if(myDarkData != NULL)
		free(myDarkData);
	
	// Allocate the array data
	myDarkData = (IMAGETYPE *) new BYTE[sizeof(IMAGETYPE)*width*height];
	
	// Copy over the data
	for(int i = 0; i < width; i++)
		for(int j = 0; j < height; j++)
			myDarkData[i*height+j] = data[i*height + j];}


/*!
    \fn CameraAPOGEE::GetDarkDataPointer()
 */
IMAGETYPE * CameraAPOGEE::GetDarkDataPointer()
{
	return myDarkData;
}


/*!
    \fn CameraAPOGEE::SetGain(float gain)
 */
void CameraAPOGEE::SetGain(float gain)
{
	if(gain < 1)
		gain = 1;
	if(gain > 12)
		gain = 12;
	
	cameraDevice->m_Gain = gain;
	
	myUSBgain = gain;
}


/*!
    \fn CameraAPOGEE::GetGain()
 */
float CameraAPOGEE::GetGain()
{
	return myUSBgain;
}

/*!
    \fn CameraAPOGEE::SetTemperatureRegulationOn()
 */
void CameraAPOGEE::SetTemperatureRegulationOn()
{
	if(myCameraOK)
	{
		// Set the temperature setting
		fflush(NULL);
		cameraDevice->write_CoolerMode(0);
		cameraDevice->write_CoolerMode(1);
		cameraDevice->write_CoolerSetPoint(myCCDTemp);
		
	mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation ON");
	isTempRegulating = true;
	}
}

/*!
    \fn CameraAPOGEE::SetTemperatureRegulationOff()
 */
void CameraAPOGEE::SetTemperatureRegulationOff()
{
	if(myCameraOK)
	{
		// Set the temperature setting
		fflush(NULL);
		cameraDevice->write_CoolerSetPoint(25.0);
		cameraDevice->write_CoolerMode(0);
		
		mySystem.myConfiguration.myLogger.WriteSystemLogString("Temperature regulation off");
		isTempRegulating = false;
	}
}


/*!
    \fn CameraAPOGEE::GetAutoTemp()
 */
bool CameraAPOGEE::GetAutoTemp()
{
	return myAutoTempRegulation;
}


/*!
    \fn CameraAPOGEE::SetAutoTemp(bool auteReg)
 */
void CameraAPOGEE::SetAutoTemp(bool autoReg)
{
	myAutoTempRegulation = autoReg;
}

// Convert a string to a decimal or hexadecimal integer
unsigned short hextoi(char *instr)
{
    unsigned short val, tot = 0;
	bool IsHEX = false;
	
	long n = strlen( instr );
	if ( n > 1 )
	{	// Look for hex format e.g. 8Fh, A3H or 0x5D
		if ( instr[ n - 1 ] == 'h' || instr[ n - 1 ] == 'H' )
			IsHEX = true;
		else if ( *instr == '0' && *(instr+1) == 'x' )
		{
			IsHEX = true;
			instr += 2;
		}
	}
	
	if ( IsHEX )
	{
		while (instr && *instr && isxdigit(*instr))
		{
			val = *instr++ - '0';
			if (9 < val)
				val -= 7;
			tot <<= 4;
			tot |= (val & 0x0f);
		}
	}
	else
		tot = atoi( instr );
	
	return tot;
}

// Trim trailing spaces from a string
void trimstr(char *s)
{
    char *p;

    p = s + (strlen(s) - 1);
    while (isspace(*p))
        p--;
    *(++p) = '\0';
}

//-------------------------------------------------------------
// CfgGet
//
// Retrieve a parameter from an INI file. Returns a status code
// and the paramter string in retbuff.
//-------------------------------------------------------------
bool CfgGet ( FILE* inifp,
               char  *inisect,
               char  *iniparm,
               char  *retbuff,
               short bufflen,
               short *parmlen)
{
    short gotsect;
    char  tbuf[256];
    char  *ss, *eq, *ps, *vs, *ptr;

	rewind( inifp );

    // find the target section

    gotsect = 0;
    while (fgets(tbuf,256,inifp) != NULL) {
        if ((ss = strchr(tbuf,'[')) != NULL) {
            if (strncasecmp(ss+1,inisect,strlen(inisect)) == 0) {
                gotsect = 1;
                break;
                }
            }
        }

    if (!gotsect) {                             // section not found
        return false;
        }

    while (fgets(tbuf,256,inifp) != NULL) {     // find parameter in sect

        if ((ptr = strrchr(tbuf,'\n')) != NULL) // remove newline if there
            *ptr = '\0';

        ps = tbuf+strspn(tbuf," \t");           // find the first non-blank

        if (*ps == ';')                         // Skip line if comment
            continue;

        if (*ps == '[') {                       // Start of next section
            return false;
            }

        eq = strchr(ps,'=');                    // Find '=' sign in string

        if (eq)
            vs = eq + 1 + strspn(eq+1," \t");   // Find start of value str
        else
            continue;

        // found the target parameter

        if (strncasecmp(ps,iniparm,strlen(iniparm)) == 0) {

            if ((ptr = strchr(vs,';')) != NULL) // cut off an EOL comment
                *ptr = '\0';

            if (short(strlen(vs)) > bufflen - 1) {// not enough buffer space
                strncpy(retbuff,vs,bufflen - 1);
                retbuff[bufflen - 1] = '\0';    // put EOL in string
                *parmlen = bufflen;
                if (verbose) printf("Configuration %s.%s = %s\n",inisect,iniparm,retbuff);
                return true;
                }
            else {
                strcpy(retbuff,vs);             // got it
                trimstr(retbuff);               // trim any trailing blanks
                *parmlen = strlen(retbuff);
                if (verbose) printf("Configuration %s.%s = %s\n",inisect,iniparm,retbuff);
                return true;
                }
            }
        }

    return false;                         // parameter not found
}

// Initializes internal variables to their default value and reads the parameters in the
// specified INI file. Then initializes the camera using current settings. If BaseAddress
// or RegOffset parameters are specified (not equal to -1) then one or both of these
// values are used for the m_BaseAddress and m_RegisterOffset properties overriding those
// settings in the INI file.
long CameraAPOGEE::config_load( char* cfgname, short BaseAddress, short RegOffset )
{
    short plen;
    char retbuf[256];
    int m_BaseAddress;
	
    if ((strlen(cfgname) == 0) || (cfgname[0] == '\0')) return CCD_OPEN_CFGNAME;
	
    // attempt to open INI file
    FILE* inifp = NULL;
	
	if ((inifp = fopen(cfgname,"r")) == NULL) return CCD_OPEN_CFGDATA;
	
	
    // System
    if (CfgGet (inifp, "system", "interface", retbuf, sizeof(retbuf), &plen))
	{
		cameraDevice = new CCameraIO();
		
		if ( cameraDevice == NULL )
		{
			fclose( inifp );
			return CCD_OPEN_ALLOC;
		}
	}
    else
	{
		fclose( inifp );
		return CCD_OPEN_CFGDATA;
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// Settings which are stored in a class member (not in firmware) are already set
	// to a default value in the constructor. Settings accessed by get/put functions
	// must be set to a default value in this routine, after the base address and
	// communication protocal is setup.
	
	/////////////////////////////////////////////////////////////////////////////////
	// These settings must done first since they affect communication with the camera
	// In the Linux drivers this is taken care of by the /dev/apppi0 device  file
	// so we just set a dummy m_BaseAddress instead
	if ( BaseAddress == -1 )
	{
		if (CfgGet (inifp, "system", "base", retbuf, sizeof(retbuf), &plen))
			m_BaseAddress = hextoi(retbuf) & 0xFFF;
		else
		{
			fclose( inifp );
			delete cameraDevice;
			cameraDevice = NULL;
			return CCD_OPEN_CFGDATA;           // base address MUST be defined
		}
	}
	else
		m_BaseAddress = BaseAddress & 0xFFF;
	
	if ( RegOffset == -1 )
	{
		if (CfgGet (inifp, "system", "reg_offset", retbuf, sizeof(retbuf), &plen))
		{
			unsigned short val = hextoi(retbuf);
			cameraDevice->m_RegisterOffset = val & 0xF0;
		}
	}
	else
	{
		if ( RegOffset >= 0x0 && RegOffset <= 0xF0 ) cameraDevice->m_RegisterOffset = RegOffset & 0xF0;
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// Necessary geometry settings
	
	if (CfgGet (inifp, "geometry", "rows", retbuf, sizeof(retbuf), &plen))
	{
		short val = hextoi(retbuf);
		if ( val >= 1 && val <= MAXTOTALROWS ) cameraDevice->m_Rows = val;
	}
	else
	{
		fclose( inifp );
		delete cameraDevice;
		cameraDevice = NULL;
		return CCD_OPEN_CFGDATA;           // rows MUST be defined
	}
	
	if (CfgGet (inifp, "geometry", "columns", retbuf, sizeof(retbuf), &plen))
	{
		short val = hextoi(retbuf);
		if ( val >= 1 && val <= MAXTOTALCOLUMNS ) cameraDevice->m_Columns = val;
	}
	else
	{
		fclose( inifp );
		delete cameraDevice;
		cameraDevice = NULL;
		return CCD_OPEN_CFGDATA;           // columns MUST be defined
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	
	if (CfgGet (inifp, "system", "pp_repeat", retbuf, sizeof(retbuf), &plen))
	{
		short val = hextoi( retbuf );
		if ( val > 0 && val <= 1000 ) cameraDevice->m_PPRepeat = val;
	}
	
	/////////////////////////////////////////////////////////////////////////////////
	// First actual communication with camera if in PPI mode
	if ( !cameraDevice->InitDriver(0) )
	{
		delete cameraDevice;
		cameraDevice = NULL;
		fclose( inifp );
		return CCD_OPEN_LOOPTST;
	}
	/////////////////////////////////////////////////////////////////////////////////
	// First actual communication with camera if in ISA mode
	cameraDevice->Reset();	// Read in command register to set shadow register known state
	/////////////////////////////////////////////////////////////////////////////////
	
	if (CfgGet (inifp, "system", "cable", retbuf, sizeof(retbuf), &plen))
	{
		if (!strcmp("LONG",retbuf))
			cameraDevice->write_LongCable( true );
		else if ( !strcmp("SHORT",retbuf) )
			cameraDevice->write_LongCable( false );
	}
	else
		cameraDevice->write_LongCable( false );	// default
	if ( !cameraDevice->read_Present() )
	{
		delete cameraDevice;
		cameraDevice = NULL;
		fclose( inifp );
		
		return CCD_OPEN_LOOPTST;
	}
	/////////////////////////////////////////////////////////////////////////////////
	// Set default setting and read other settings from ini file
	
	cameraDevice->write_UseTrigger( false );
	cameraDevice->write_ForceShutterOpen( false );
	
	if (CfgGet (inifp, "system", "high_priority", retbuf, sizeof(retbuf), &plen))
	{
		if (!strcmp("ON",retbuf) || !strcmp("TRUE",retbuf) || !strcmp("1",retbuf))
		{
			cameraDevice->m_HighPriority = true;
		}
		else if (!strcmp("OFF",retbuf) || !strcmp("FALSE",retbuf) || !strcmp("0",retbuf))
		{
			cameraDevice->m_HighPriority = false;
		}
	}
	
	if (CfgGet (inifp, "system", "data_bits", retbuf, sizeof(retbuf), &plen))
	{
		short val = hextoi( retbuf );
		if ( val >= 8 && val <= 18 ) cameraDevice->m_DataBits = val;
	}
	
	if (CfgGet (inifp, "system", "sensor", retbuf, sizeof(retbuf), &plen))
	{
		if ( strcmp( "ccd", retbuf ) == 0 )
		{
			cameraDevice->m_SensorType = Camera_SensorType_CCD;
		}
		else if ( strcmp( "cmos", retbuf ) == 0 )
		{
			cameraDevice->m_SensorType = Camera_SensorType_CMOS;
		}
	}
	
    if (CfgGet (inifp, "system", "mode", retbuf, sizeof(retbuf), &plen))
	{
        unsigned short val = hextoi(retbuf) & 0xF;
        cameraDevice->write_Mode( val );
    }
	else
		cameraDevice->write_Mode( 0 );			// default
	
    if (CfgGet (inifp, "system", "test", retbuf, sizeof(retbuf), &plen))
	{
        unsigned short val = hextoi(retbuf) & 0xF;
        cameraDevice->write_TestBits( val );
    }
	else
		cameraDevice->write_TestBits( 0 );		//default
	
    if (CfgGet (inifp, "system", "test2", retbuf, sizeof(retbuf), &plen))
	{
        unsigned short val = hextoi(retbuf) & 0xF;
        cameraDevice->write_Test2Bits( val );
    }
	else
		cameraDevice->write_Test2Bits( 0 );	// default
	
	cameraDevice->write_FastReadout( false );	//default
	
    if (CfgGet (inifp, "system", "shutter_speed", retbuf, sizeof(retbuf), &plen))
	{
        if (!strcmp("normal",retbuf))
		{
			cameraDevice->m_FastShutter = false;
			cameraDevice->m_MaxExposure = 10485.75;
			cameraDevice->m_MinExposure = 0.01;
		}
		else if (!strcmp("fast",retbuf))
		{
			cameraDevice->m_FastShutter = true;
			cameraDevice->m_MaxExposure = 1048.575;
			cameraDevice->m_MinExposure = 0.001;
		}
		else if ( !strcmp("dual",retbuf))
		{
			cameraDevice->m_FastShutter = true;
			cameraDevice->m_MaxExposure = 10485.75;
			cameraDevice->m_MinExposure = 0.001;
		}
    }
	
    if (CfgGet (inifp, "system", "shutter_bits", retbuf, sizeof(retbuf), &plen))
	{
		unsigned short val = hextoi(retbuf);
		cameraDevice->m_FastShutterBits_Mode = val & 0x0F;
		cameraDevice->m_FastShutterBits_Test = ( val & 0xF0 ) >> 4;
	}
	
    if (CfgGet (inifp, "system", "maxbinx", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXHBIN ) cameraDevice->m_MaxBinX = val;
	}
	
    if (CfgGet (inifp, "system", "maxbiny", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXVBIN ) cameraDevice->m_MaxBinY = val;
	}
	
    if (CfgGet (inifp, "system", "guider_relays", retbuf, sizeof(retbuf), &plen))
	{
		if (!strcmp("ON",retbuf) || !strcmp("TRUE",retbuf) || !strcmp("1",retbuf))
		{
			cameraDevice->m_GuiderRelays = true;
		}
		else if (!strcmp("OFF",retbuf) || !strcmp("FALSE",retbuf) || !strcmp("0",retbuf))
		{
			cameraDevice->m_GuiderRelays = false;
		}
	}
	
    if (CfgGet (inifp, "system", "timeout", retbuf, sizeof(retbuf), &plen))
	{
        double val = atof(retbuf);
		if ( val >= 0.0 && val <= 10000.0 ) cameraDevice->m_Timeout = val;
    }
	
	/////////////////////////////////////////////////////////////////////////////////
    // Geometry
	
    if (CfgGet (inifp, "geometry", "bic", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXCOLUMNS ) cameraDevice->m_BIC = val;
	}
	
    if (CfgGet (inifp, "geometry", "bir", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXROWS ) cameraDevice->m_BIR = val;
	}
	
	if (CfgGet (inifp, "geometry", "skipc", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 0 && val <= MAXCOLUMNS ) cameraDevice->m_SkipC = val;
	}
	
	if (CfgGet (inifp, "geometry", "skipr", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 0 && val <= MAXROWS ) cameraDevice->m_SkipR = val;
	}
	
    if (CfgGet (inifp, "geometry", "imgcols", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXTOTALCOLUMNS ) cameraDevice->m_ImgColumns = val;
	}
	else
		cameraDevice->m_ImgColumns = cameraDevice->m_Columns - cameraDevice->m_BIC - cameraDevice->m_SkipC;
	
    if (CfgGet (inifp, "geometry", "imgrows", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXTOTALROWS ) cameraDevice->m_ImgRows = val;
	}
	else
		cameraDevice->m_ImgRows = cameraDevice->m_Rows - cameraDevice->m_BIR - cameraDevice->m_SkipR;
	
    if (CfgGet (inifp, "geometry", "hflush", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXHBIN ) cameraDevice->m_HFlush = val;
	}
	
    if (CfgGet (inifp, "geometry", "vflush", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
        if ( val >= 1 && val <= MAXVBIN ) cameraDevice->m_VFlush = val;
	}
	
	// Default to full frame
	cameraDevice->m_NumX = cameraDevice->m_ImgColumns;
	cameraDevice->m_NumY = cameraDevice->m_ImgRows;
	
	/////////////////////////////////////////////////////////////////////////////////
	// Temperature
	
    if (CfgGet (inifp, "temp", "control", retbuf, sizeof(retbuf), &plen))
	{
		if (!strcmp("ON",retbuf) || !strcmp("TRUE",retbuf) || !strcmp("1",retbuf))
		{
			cameraDevice->m_TempControl = true;
		}
		else if (!strcmp("OFF",retbuf) || !strcmp("FALSE",retbuf) || !strcmp("0",retbuf))
		{
			cameraDevice->m_TempControl = false;
		}
    }
	
    if (CfgGet (inifp, "temp", "cal", retbuf, sizeof(retbuf), &plen))
	{
        short val = hextoi(retbuf);
		if ( val >= 1 && val <= 255 ) cameraDevice->m_TempCalibration = val;
    }
	
    if (CfgGet (inifp, "temp", "scale", retbuf, sizeof(retbuf), &plen))
	{
        double val = atof(retbuf);
		if ( val >= 1.0 && val <= 10.0 ) cameraDevice->m_TempScale = val;
    }
	
    if (CfgGet (inifp, "temp", "target", retbuf, sizeof(retbuf), &plen))
	{
        double val = atof(retbuf);
        if ( val >= -60.0 && val <= 40.0 )
			cameraDevice->write_CoolerSetPoint( val );
		else
			cameraDevice->write_CoolerSetPoint( -10.0 );
    }
	else
		cameraDevice->write_CoolerSetPoint( -10.0 );	//default
	
	/////////////////////////////////////////////////////////////////////////////////
	// CCD
	
	if (CfgGet (inifp, "ccd", "sensor", retbuf, sizeof(retbuf), &plen))
	{
		if ( plen > 256 ) plen = 256;
		memcpy( cameraDevice->m_Sensor, retbuf, plen );
    }
	
	if (CfgGet (inifp, "ccd", "color", retbuf, sizeof(retbuf), &plen))
	{
		if (!strcmp("ON",retbuf) || !strcmp("TRUE",retbuf) || !strcmp("1",retbuf))
		{
			cameraDevice->m_Color = true;
		}
		else if (!strcmp("OFF",retbuf) || !strcmp("FALSE",retbuf) || !strcmp("0",retbuf))
		{
			cameraDevice->m_Color = false;
		}
    }
	
    if (CfgGet (inifp, "ccd", "noise", retbuf, sizeof(retbuf), &plen))
	{
		cameraDevice->m_Noise = atof( retbuf );
    }
	
	if (CfgGet (inifp, "ccd", "gain", retbuf, sizeof(retbuf), &plen))
	{
		cameraDevice->m_Gain = atof( retbuf );
    }
	
    if (CfgGet (inifp, "ccd", "pixelxsize", retbuf, sizeof(retbuf), &plen))
	{
		cameraDevice->m_PixelXSize = atof( retbuf );
    }
	
    if (CfgGet (inifp, "ccd", "pixelysize", retbuf, sizeof(retbuf), &plen))
	{
		cameraDevice->m_PixelYSize = atof( retbuf );
    }
	
	fclose( inifp );
    return CCD_OPEN_NOERR;
}
